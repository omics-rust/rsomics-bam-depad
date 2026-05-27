//! `samtools depad` port: convert padded BAM alignments to unpadded BAM.
//!
//! Padded BAM arises from multiple-sequence assembly tools (e.g. CAF, ACE
//! format pipelines). The padded reference has gap columns (`*`/`-`); the
//! BAM header's `@SQ LN` values reflect the padded (longer) lengths; read
//! CIGARs and POS values are in padded coordinates. This tool converts to
//! unpadded coordinates so the result can be used with a standard reference.
//!
//! ## Algorithm (mirrors `padding.c:bam_pad2unpad`)
//!
//! 1. For each reference, build a **position map** (`posmap`) that converts a
//!    padded 0-based position to the corresponding unpadded 0-based position.
//!    The source of truth is either:
//!    - An **embedded reference record** in the BAM stream: a read whose `QNAME`
//!      equals the `@SQ` name it is aligned to, at `POS == 0`. Its CIGAR (M/D
//!      ops; D = gap/pad site) plus SEQ define the padded reference.
//!    - A **padded FASTA** supplied via `-T`: gap characters (`*` or `-`) mark
//!      pad positions.
//! 2. For each query read, reconstruct its CIGAR by comparing the query bases
//!    (from M/D positions) against the padded reference base by base:
//!    - query base ≠ 0, ref base ≠ 0  → `M`
//!    - query base ≠ 0, ref base  = 0  → `I`
//!    - query base  = 0, ref base ≠ 0  → `D`
//!    - query base  = 0, ref base  = 0  → `P`
//! 3. Remove redundant `P` operators that occur between `M`/`D` operators
//!    (e.g. `5M2P10M` → `15M`). Leading `P` ops are preserved.
//! 4. Update `POS` (and `PNEXT` if on the same reference) via `posmap`.
//! 5. Recalculate `BIN` from the new POS + end-pos pair.
//!
//! Input reads must not contain `P` or `I` CIGAR operators.

use std::io::Write;
use std::num::NonZero;
use std::path::Path;

use noodles::bam;
use noodles::bgzf;
use noodles::sam;
use noodles::sam::Header;
use noodles::sam::alignment::io::Write as AlnWrite;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::Program;
use noodles::sam::header::record::value::map::program::tag as program_tag;
use rsomics_bamio::raw::{self, RawRecord};
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

// BAM CIGAR op codes (low nibble of the packed u32): M=0 I=1 D=2 N=3 S=4 H=5 P=6 ==7 X=8.
const CIGAR_M: u8 = 0;
const CIGAR_I: u8 = 1;
const CIGAR_D: u8 = 2;
const CIGAR_N: u8 = 3;
const CIGAR_S: u8 = 4;
const CIGAR_H: u8 = 5;
const CIGAR_P: u8 = 6;

// FLAG bits (SAMv1 §1.4).
const FLAG_UNMAPPED: u16 = 0x4;

// BAM payload layout constants (offsets from start of payload, after block_size).
const POS: usize = 4;
const L_READ_NAME: usize = 8;
const N_CIGAR: usize = 12;
const FLAG: usize = 14;
const L_SEQ: usize = 16;
const NEXT_REF_ID: usize = 20;
const NEXT_POS: usize = 24;
const FIXED_HEAD: usize = 32;

/// htslib `bam_reg2bin(b, e)` = `hts_reg2bin(b, e, 14, 5)`.
fn reg2bin(beg: i64, end: i64) -> u16 {
    let end = if end > 0 { end - 1 } else { 0 };
    if beg >> 14 == end >> 14 {
        return (((1 << 15) - 1) / 7 + (beg >> 14)) as u16;
    }
    if beg >> 17 == end >> 17 {
        return (((1 << 12) - 1) / 7 + (beg >> 17)) as u16;
    }
    if beg >> 20 == end >> 20 {
        return (((1 << 9) - 1) / 7 + (beg >> 20)) as u16;
    }
    if beg >> 23 == end >> 23 {
        return (((1 << 6) - 1) / 7 + (beg >> 23)) as u16;
    }
    if beg >> 26 == end >> 26 {
        return (1 + (beg >> 26)) as u16;
    }
    0
}

/// Build a position map: `posmap[padded_pos] = unpadded_pos`.
/// `ref_seq[i] == 0` indicates a pad column; non-zero is a real base.
fn build_posmap(ref_seq: &[u8]) -> Vec<i32> {
    let mut map = Vec::with_capacity(ref_seq.len());
    let mut k: i32 = 0;
    for &base in ref_seq {
        map.push(k);
        if base != 0 {
            k += 1;
        }
    }
    map
}

/// Decode the padded sequence of a read into `out`: each output element is
/// either 0 (a gap from a D op) or a non-zero nibble code (from M/=/X ops).
/// S ops consume query bases without emitting. H/P ops emit nothing and consume
/// nothing. N ops are treated as D with a warning (mirrors samtools).
/// I ops are rejected (must not appear in input).
fn decode_seq(record: &[u8], out: &mut Vec<u8>, qname: &[u8]) -> Result<()> {
    let name_len = usize::from(record[L_READ_NAME]);
    let n_cigar = usize::from(u16::from_le_bytes([record[N_CIGAR], record[N_CIGAR + 1]]));

    let cigar_start = FIXED_HEAD + name_len;
    let seq_start = cigar_start + n_cigar * 4;

    out.clear();
    let mut q_idx = 0usize;

    for ci in 0..n_cigar {
        let raw_op = u32::from_le_bytes(
            record[cigar_start + ci * 4..cigar_start + ci * 4 + 4]
                .try_into()
                .unwrap(),
        );
        let op = (raw_op & 0xf) as u8;
        let ol = (raw_op >> 4) as usize;

        match op {
            CIGAR_M | 7 | 8 => {
                for _ in 0..ol {
                    let nibble = if q_idx.is_multiple_of(2) {
                        record[seq_start + q_idx / 2] >> 4
                    } else {
                        record[seq_start + q_idx / 2] & 0x0f
                    };
                    // 0 in seq_nt16 means '='; treat as a non-zero sentinel (has a base).
                    out.push(if nibble == 0 { 0xff } else { nibble });
                    q_idx += 1;
                }
            }
            CIGAR_D => {
                for _ in 0..ol {
                    out.push(0);
                }
            }
            CIGAR_N => {
                eprintln!(
                    "[depad] WARNING: CIGAR op N treated as op D in read {}",
                    String::from_utf8_lossy(qname)
                );
                for _ in 0..ol {
                    out.push(0);
                }
            }
            CIGAR_S => {
                q_idx += ol;
            }
            CIGAR_H | CIGAR_P => {}
            CIGAR_I => {
                return Err(RsomicsError::InvalidInput(format!(
                    "[depad] ERROR: Didn't expect CIGAR op I in read {}",
                    String::from_utf8_lossy(qname)
                )));
            }
            _ => {
                return Err(RsomicsError::InvalidInput(format!(
                    "[depad] ERROR: Unknown CIGAR op {op} in read {}",
                    String::from_utf8_lossy(qname)
                )));
            }
        }
    }

    Ok(())
}

/// Return (hard_clip, soft_clip) lengths at the 5' end of the CIGAR.
fn leading_clips(cigar: &[(u8, u32)]) -> (u32, u32) {
    let mut h = 0u32;
    let mut s = 0u32;
    let mut i = 0usize;
    if !cigar.is_empty() && cigar[0].0 == CIGAR_H {
        h = cigar[0].1;
        i = 1;
    }
    if i < cigar.len() && cigar[i].0 == CIGAR_S {
        s = cigar[i].1;
    }
    (h, s)
}

/// Return (soft_clip, hard_clip) lengths at the 3' end of the CIGAR.
fn trailing_clips(cigar: &[(u8, u32)]) -> (u32, u32) {
    let n = cigar.len();
    if n == 0 {
        return (0, 0);
    }
    let mut h = 0u32;
    let mut s = 0u32;
    let mut i = n;
    if cigar[n - 1].0 == CIGAR_H {
        h = cigar[n - 1].1;
        i -= 1;
    }
    if i > 0 && cigar[i - 1].0 == CIGAR_S {
        s = cigar[i - 1].1;
    }
    (s, h)
}

/// Pack (length, op) into the BAM packed u32: length in [31:4], op in [3:0].
#[inline]
fn cigar_gen(len: u32, op: u8) -> u32 {
    (len << 4) | u32::from(op)
}

/// Decode the original CIGAR ops from the BAM payload.
fn decode_cigar(record: &[u8]) -> Vec<(u8, u32)> {
    let name_len = usize::from(record[L_READ_NAME]);
    let n_cigar = usize::from(u16::from_le_bytes([record[N_CIGAR], record[N_CIGAR + 1]]));
    let cigar_start = FIXED_HEAD + name_len;
    (0..n_cigar)
        .map(|ci| {
            let raw = u32::from_le_bytes(
                record[cigar_start + ci * 4..cigar_start + ci * 4 + 4]
                    .try_into()
                    .unwrap(),
            );
            ((raw & 0xf) as u8, raw >> 4)
        })
        .collect()
}

/// Rebuild the CIGAR for a query read after depadding. Mirrors
/// `bam_pad2unpad`'s CIGAR-reconstruction loop (padding.c).
fn rebuild_cigar(record: &[u8], q: &[u8], ref_seq: &[u8], pos: usize) -> Vec<u32> {
    let orig_cigar = decode_cigar(record);
    let mut out: Vec<u32> = Vec::with_capacity(q.len() + 4);

    // Passthrough leading H/S clips.
    let (lead_h, lead_s) = leading_clips(&orig_cigar);
    if lead_h > 0 {
        out.push(cigar_gen(lead_h, CIGAR_H));
        if lead_s > 0 {
            out.push(cigar_gen(lead_s, CIGAR_S));
        }
    } else if lead_s > 0 {
        out.push(cigar_gen(lead_s, CIGAR_S));
    }

    if q.is_empty() {
        let (trail_s, trail_h) = trailing_clips(&orig_cigar);
        if trail_s > 0 {
            out.push(cigar_gen(trail_s, CIGAR_S));
        }
        if trail_h > 0 {
            out.push(cigar_gen(trail_h, CIGAR_H));
        }
        return out;
    }

    // Determine the per-base CIGAR operator for each aligned position.
    let mut ops: Vec<u8> = Vec::with_capacity(q.len());
    for (i, &qb) in q.iter().enumerate() {
        let k = pos + i;
        let rb = ref_seq[k];
        ops.push(match (qb != 0, rb != 0) {
            (true, true) => CIGAR_M,
            (true, false) => CIGAR_I,
            (false, true) => CIGAR_D,
            (false, false) => CIGAR_P,
        });
    }

    // Count leading pad columns immediately before `pos` to handle the case
    // where a read starts at or just after a pad column. This mirrors the
    // `include any pads if starts with an insert` block in padding.c.
    let lead_pad_k: u32 = if matches!(ops[0], CIGAR_I | CIGAR_P) {
        let mut scan = pos;
        let mut k = 0u32;
        while scan > 0 && ref_seq[scan - 1] == 0 {
            k += 1;
            scan -= 1;
        }
        k
    } else {
        0
    };

    if lead_pad_k > 0 && ops[0] == CIGAR_I {
        out.push(cigar_gen(lead_pad_k, CIGAR_P));
    }

    // Compress consecutive identical ops into run-length CIGAR entries.
    let start_len: u32 = if ops[0] == CIGAR_P { lead_pad_k + 1 } else { 1 };
    let mut cur_op = ops[0];
    let mut cur_len = start_len;

    for &op in &ops[1..] {
        if op == cur_op {
            cur_len += 1;
        } else {
            out.push(cigar_gen(cur_len, cur_op));
            cur_op = op;
            cur_len = 1;
        }
    }
    out.push(cigar_gen(cur_len, cur_op));

    // Remove redundant P operators between M/D operators
    // (e.g. 5M 2P 10M → 15M). Mirrors samtools' post-processing pass.
    remove_redundant_p(&mut out);

    // Passthrough trailing H/S clips.
    let (trail_s, trail_h) = trailing_clips(&orig_cigar);
    if trail_h > 0 {
        if trail_s > 0 {
            out.push(cigar_gen(trail_s, CIGAR_S));
        }
        out.push(cigar_gen(trail_h, CIGAR_H));
    } else if trail_s > 0 {
        out.push(cigar_gen(trail_s, CIGAR_S));
    }

    out
}

/// Remove redundant P operators from the CIGAR in place. A P between two M or D
/// operators is redundant (padding.c post-processing pass): mark as `0` and
/// compact. If the flanking ops are the same, merge their lengths.
fn remove_redundant_p(cigar: &mut Vec<u32>) {
    let n = cigar.len();
    if n < 3 {
        return;
    }
    let mut i = 2;
    while i < cigar.len() {
        let pre_op = (cigar[i - 2] & 0xf) as u8;
        let mid_op = (cigar[i - 1] & 0xf) as u8;
        let post_op = (cigar[i] & 0xf) as u8;
        if mid_op == CIGAR_P
            && (pre_op == CIGAR_M || pre_op == CIGAR_D)
            && (post_op == CIGAR_M || post_op == CIGAR_D)
        {
            cigar[i - 1] = 0; // mark P as 0 (= 0M, filtered below)
            if pre_op == post_op {
                let combined = (cigar[i - 2] >> 4) + (cigar[i] >> 4);
                cigar[i] = cigar_gen(combined, post_op);
                cigar[i - 2] = 0;
            }
        }
        i += 1;
    }
    cigar.retain(|&op| op != 0);
}

/// Compute the reference span of a CIGAR (M/D/N/=/X ops).
fn ref_span(cigar: &[u32]) -> i64 {
    cigar
        .iter()
        .filter(|&&op| matches!((op & 0xf) as u8, CIGAR_M | CIGAR_D | CIGAR_N | 7 | 8))
        .map(|&op| i64::from(op >> 4))
        .sum()
}

/// Replace the CIGAR in a raw BAM payload, returning a new byte vector. The
/// name/seq/qual/aux regions are copied verbatim; only CIGAR and `n_cigar`
/// change.
fn replace_cigar(record: &[u8], new_cigar: &[u32]) -> Vec<u8> {
    let name_len = usize::from(record[L_READ_NAME]);
    let old_n_cigar = usize::from(u16::from_le_bytes([record[N_CIGAR], record[N_CIGAR + 1]]));
    let cigar_start = FIXED_HEAD + name_len;
    let after_old = cigar_start + old_n_cigar * 4;

    let mut out = Vec::with_capacity(
        FIXED_HEAD + name_len + new_cigar.len() * 4 + (record.len() - after_old),
    );
    out.extend_from_slice(&record[..FIXED_HEAD + name_len]);
    for &op in new_cigar {
        out.extend_from_slice(&op.to_le_bytes());
    }
    out.extend_from_slice(&record[after_old..]);
    let nc = new_cigar.len() as u16;
    out[N_CIGAR..N_CIGAR + 2].copy_from_slice(&nc.to_le_bytes());
    out
}

/// Read a padded FASTA (gap chars = `*` or `-`) and return the sequence as a
/// nibble array: 0 = gap, non-zero = real base (seq_nt16 code).
fn load_fasta_ref(fasta_path: &Path, ref_name: &str, padded_len: usize) -> Result<Vec<u8>> {
    use std::io::{BufRead, BufReader};
    let f = std::fs::File::open(fasta_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", fasta_path.display())))?;
    let mut reader = BufReader::new(f);
    let mut line = String::new();
    let mut in_target = false;
    let mut seq: Vec<u8> = Vec::with_capacity(padded_len);

    loop {
        line.clear();
        let n = reader.read_line(&mut line).map_err(RsomicsError::Io)?;
        if n == 0 {
            break;
        }
        let trimmed = line.trim_end();
        if let Some(name) = trimmed.strip_prefix('>') {
            let this_name = name.split_whitespace().next().unwrap_or("");
            if in_target {
                break; // started a new record after ours
            }
            in_target = this_name == ref_name;
        } else if in_target {
            for b in trimmed.bytes() {
                let code = if b == b'-' || b == b'*' {
                    0u8
                } else {
                    let c = SEQ_NT16_TABLE[b as usize];
                    if c == 16 {
                        return Err(RsomicsError::InvalidInput(format!(
                            "[depad] invalid base '{}' (ASCII {b}) in reference {ref_name}",
                            b as char
                        )));
                    }
                    if c == 0 { 0xff } else { c } // '=' → non-zero sentinel
                };
                seq.push(code);
            }
        }
    }

    if seq.len() != padded_len {
        return Err(RsomicsError::InvalidInput(format!(
            "[depad] reference {ref_name}: FASTA length {} != expected {padded_len}",
            seq.len()
        )));
    }
    Ok(seq)
}

/// htslib `seq_nt16_table`: ASCII byte → 4-bit nucleotide code.
#[rustfmt::skip]
const SEQ_NT16_TABLE: [u8; 256] = [
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
     1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8, 8, 7, 9, 15,10,15,15, 15,15,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8, 8, 7, 9, 15,10,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
];

/// Depad options.
#[derive(Debug, Clone, Default)]
pub struct DepadOpts {
    /// `-s`: emit SAM (default: BAM).
    pub sam_output: bool,
    /// `-u`: uncompressed BAM.
    pub uncompressed: bool,
    /// `-1`: level-1 BAM compression.
    pub fast_compression: bool,
    /// `-T FILE`: padded FASTA reference.
    pub reference: Option<std::path::PathBuf>,
    /// `--no-PG`: suppress @PG line.
    pub no_pg: bool,
}

#[derive(Debug, Default, Clone, Serialize)]
pub struct DepadStats {
    pub records: u64,
    pub embedded_refs: u64,
}

/// Write a single raw-payload record (with a 4-byte block_size prefix) to `out`.
fn write_bytes<W: Write>(out: &mut W, bytes: &[u8]) -> Result<()> {
    let block_size = u32::try_from(bytes.len())
        .map_err(|e| RsomicsError::InvalidInput(format!("record too large: {e}")))?;
    out.write_all(&block_size.to_le_bytes())
        .map_err(RsomicsError::Io)?;
    out.write_all(bytes).map_err(RsomicsError::Io)?;
    Ok(())
}

/// Entry point called from `cli.rs`.
pub fn depad(input: &Path, output_path: Option<&Path>, opts: &DepadOpts) -> Result<DepadStats> {
    let workers = NonZero::<usize>::MIN;
    let mut reader = rsomics_bamio::open_with_workers(input, workers)?;
    let mut header: Header = reader.read_header().map_err(RsomicsError::Io)?;

    if !opts.no_pg {
        let program = Map::<Program>::builder()
            .insert(program_tag::NAME, "rsomics-bam-depad")
            .insert(program_tag::VERSION, env!("CARGO_PKG_VERSION"))
            .build()
            .map_err(|e| RsomicsError::InvalidInput(format!("building @PG: {e}")))?;
        header
            .programs_mut()
            .add("rsomics-bam-depad", program)
            .map_err(RsomicsError::Io)?;
    }

    let ref_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|k| String::from_utf8_lossy(k.as_ref()).into_owned())
        .collect();

    let ref_lengths: Vec<usize> = header
        .reference_sequences()
        .values()
        .map(|sq| usize::from(sq.length()))
        .collect();

    if opts.sam_output {
        // SAM output: write to a BufWriter over stdout/file, using noodles SAM writer
        // for the header, then decode each record to SAM text.
        match output_path {
            Some(path) => {
                let file = std::fs::File::create(path).map_err(|e| {
                    RsomicsError::InvalidInput(format!("creating {}: {e}", path.display()))
                })?;
                let mut out = std::io::BufWriter::new(file);
                run_sam(
                    &mut reader,
                    &mut out,
                    header,
                    &ref_names,
                    &ref_lengths,
                    opts,
                )
            }
            None => {
                let stdout = std::io::stdout();
                let mut out = std::io::BufWriter::new(stdout.lock());
                run_sam(
                    &mut reader,
                    &mut out,
                    header,
                    &ref_names,
                    &ref_lengths,
                    opts,
                )
            }
        }
    } else {
        // BAM output: wrap in BGZF. Use the multithreaded writer for
        // compression levels; for uncompressed use NONE, for fast use level 1.
        let compression = if opts.uncompressed {
            bgzf::io::writer::CompressionLevel::NONE
        } else if opts.fast_compression {
            bgzf::io::writer::CompressionLevel::FAST
        } else {
            bgzf::io::writer::CompressionLevel::default()
        };

        match output_path {
            Some(path) => {
                let file = std::fs::File::create(path).map_err(|e| {
                    RsomicsError::InvalidInput(format!("creating {}: {e}", path.display()))
                })?;
                let bgzf_writer = bgzf::io::multithreaded_writer::Builder::default()
                    .set_compression_level(compression)
                    .set_worker_count(workers)
                    .build_from_writer(file);
                let mut bam_writer = bam::io::Writer::from(bgzf_writer);
                run_bam(
                    &mut reader,
                    &mut bam_writer,
                    header,
                    &ref_names,
                    &ref_lengths,
                    opts,
                )
            }
            None => {
                // stdout().lock() is not Send so cannot use MultithreadedWriter;
                // fall back to single-threaded BGZF writer (same as ampliconclip).
                let mut bam_writer = bam::io::Writer::new(std::io::stdout().lock());
                run_bam(
                    &mut reader,
                    &mut bam_writer,
                    header,
                    &ref_names,
                    &ref_lengths,
                    opts,
                )
            }
        }
    }
}

fn run_bam<W: Write>(
    reader: &mut rsomics_bamio::ParallelBamReader,
    writer: &mut bam::io::Writer<W>,
    header: Header,
    ref_names: &[String],
    ref_lengths: &[usize],
    opts: &DepadOpts,
) -> Result<DepadStats> {
    writer.write_header(&header).map_err(RsomicsError::Io)?;
    let inner = writer.get_mut();
    run_core(reader, inner, &header, ref_names, ref_lengths, opts)
}

fn run_sam<W: Write>(
    reader: &mut rsomics_bamio::ParallelBamReader,
    writer: &mut W,
    header: Header,
    ref_names: &[String],
    ref_lengths: &[usize],
    opts: &DepadOpts,
) -> Result<DepadStats> {
    // Write SAM header as text; scope the writer so the borrow on writer ends
    // before run_core_sam takes it.
    {
        let mut sam_writer = noodles::sam::io::Writer::new(&mut *writer);
        sam_writer.write_header(&header).map_err(RsomicsError::Io)?;
    }
    // Write records as SAM text via a scratch BAM buffer decoded via noodles.
    // This is correct but slower — SAM output is the rare path.
    run_core_sam(reader, writer, &header, ref_names, ref_lengths, opts)
}

/// Core depad loop writing raw BAM bytes (block_size + payload) to `out`.
fn run_core<W: Write>(
    reader: &mut rsomics_bamio::ParallelBamReader,
    out: &mut W,
    _header: &Header,
    ref_names: &[String],
    ref_lengths: &[usize],
    opts: &DepadOpts,
) -> Result<DepadStats> {
    let mut ref_seq: Vec<u8> = Vec::new();
    let mut posmap: Vec<i32> = Vec::new();
    let mut cur_tid: i32 = -2;

    let mut stats = DepadStats::default();
    let mut rec = RawRecord::default();
    let inner = reader.get_mut();

    loop {
        if raw::read_record(inner, &mut rec)? == 0 {
            break;
        }
        stats.records += 1;

        let bytes = rec.as_bytes();
        let flags = u16::from_le_bytes([bytes[FLAG], bytes[FLAG + 1]]);
        let tid = i32::from_le_bytes(bytes[0..4].try_into().unwrap());
        let pos = i32::from_le_bytes(bytes[POS..POS + 4].try_into().unwrap());
        let n_cigar = usize::from(u16::from_le_bytes([bytes[N_CIGAR], bytes[N_CIGAR + 1]]));
        let name_len = usize::from(bytes[L_READ_NAME]);
        let qname = &bytes[FIXED_HEAD..FIXED_HEAD + name_len.saturating_sub(1)];

        if flags & FLAG_UNMAPPED != 0 {
            let mut out_bytes = bytes.to_vec();
            // Remap POS and PNEXT through posmap if available.
            if pos >= 0 && !posmap.is_empty() && (pos as usize) < posmap.len() {
                let mapped = posmap[pos as usize];
                out_bytes[POS..POS + 4].copy_from_slice(&mapped.to_le_bytes());
            }
            let mtid =
                i32::from_le_bytes(out_bytes[NEXT_REF_ID..NEXT_REF_ID + 4].try_into().unwrap());
            let mpos = i32::from_le_bytes(out_bytes[NEXT_POS..NEXT_POS + 4].try_into().unwrap());
            if mtid == cur_tid && mpos >= 0 && !posmap.is_empty() && (mpos as usize) < posmap.len()
            {
                let mapped = posmap[mpos as usize];
                out_bytes[NEXT_POS..NEXT_POS + 4].copy_from_slice(&mapped.to_le_bytes());
            }
            write_bytes(out, &out_bytes)?;
            continue;
        }

        // Detect embedded reference: qname == ref_name AND pos == 0.
        let is_embedded_ref = pos == 0
            && tid >= 0
            && (tid as usize) < ref_names.len()
            && qname == ref_names[tid as usize].as_bytes();

        if is_embedded_ref {
            decode_seq(bytes, &mut ref_seq, qname)?;
            let expected_len = ref_lengths.get(tid as usize).copied().unwrap_or(0);
            if ref_seq.len() != expected_len {
                return Err(RsomicsError::InvalidInput(format!(
                    "[depad] embedded reference {} length {} != header LN {expected_len}",
                    String::from_utf8_lossy(qname),
                    ref_seq.len()
                )));
            }
            posmap = build_posmap(&ref_seq);
            cur_tid = tid;
            stats.embedded_refs += 1;

            let seq_len = u32::from_le_bytes(bytes[L_SEQ..L_SEQ + 4].try_into().unwrap());
            let new_cigar = [cigar_gen(seq_len, CIGAR_M)];
            let new_bytes = replace_cigar(bytes, &new_cigar);
            write_bytes(out, &new_bytes)?;
            continue;
        }

        if n_cigar == 0 {
            // No CIGAR: pass through with POS remapped.
            let mut out_bytes = bytes.to_vec();
            if pos >= 0 && !posmap.is_empty() && (pos as usize) < posmap.len() {
                let mapped = posmap[pos as usize];
                out_bytes[POS..POS + 4].copy_from_slice(&mapped.to_le_bytes());
            }
            write_bytes(out, &out_bytes)?;
            continue;
        }

        // Regular mapped read with CIGAR.
        if tid < 0 {
            return Err(RsomicsError::InvalidInput(format!(
                "[depad] read '{}' has CIGAR but no RNAME",
                String::from_utf8_lossy(qname)
            )));
        }

        if tid != cur_tid {
            let Some(ref fasta_path) = opts.reference else {
                return Err(RsomicsError::InvalidInput(format!(
                    "[depad] Missing {} embedded reference sequence (and no FASTA file)",
                    ref_names
                        .get(tid as usize)
                        .map(|s| s.as_str())
                        .unwrap_or("?")
                )));
            };
            let ref_name = ref_names.get(tid as usize).ok_or_else(|| {
                RsomicsError::InvalidInput(format!("[depad] tid {tid} out of range"))
            })?;
            let padded_len = ref_lengths.get(tid as usize).copied().unwrap_or(0);
            ref_seq = load_fasta_ref(fasta_path, ref_name, padded_len)?;
            posmap = build_posmap(&ref_seq);
            cur_tid = tid;
        }

        let mut q: Vec<u8> = Vec::new();
        decode_seq(bytes, &mut q, qname)?;

        let new_cigar = rebuild_cigar(bytes, &q, &ref_seq, pos as usize);
        let mut new_bytes = replace_cigar(bytes, &new_cigar);

        let new_pos = posmap[pos as usize];
        new_bytes[POS..POS + 4].copy_from_slice(&new_pos.to_le_bytes());

        let end_pos = new_pos as i64 + ref_span(&new_cigar);
        let bin = reg2bin(new_pos as i64, end_pos);
        new_bytes[10..12].copy_from_slice(&bin.to_le_bytes());

        let mtid = i32::from_le_bytes(new_bytes[NEXT_REF_ID..NEXT_REF_ID + 4].try_into().unwrap());
        let mpos = i32::from_le_bytes(new_bytes[NEXT_POS..NEXT_POS + 4].try_into().unwrap());
        if mtid == tid && mpos >= 0 && (mpos as usize) < posmap.len() {
            let new_mpos = posmap[mpos as usize];
            new_bytes[NEXT_POS..NEXT_POS + 4].copy_from_slice(&new_mpos.to_le_bytes());
        }

        write_bytes(out, &new_bytes)?;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(stats)
}

/// SAM-output path: same core logic but decodes final records to SAM text via
/// noodles. Used only for `-s` (rare path, correct > fast).
fn run_core_sam<W: Write>(
    reader: &mut rsomics_bamio::ParallelBamReader,
    out: &mut W,
    header: &Header,
    ref_names: &[String],
    ref_lengths: &[usize],
    opts: &DepadOpts,
) -> Result<DepadStats> {
    // Collect the depadded records into a temporary in-memory buffer using the
    // BAM binary format (without the BGZF wrapper — raw block_size + payload),
    // then decode each record via noodles and emit as SAM text.
    let mut raw_buf: Vec<u8> = Vec::new();
    let stats = run_core(reader, &mut raw_buf, header, ref_names, ref_lengths, opts)?;

    // Wrap the raw record bytes in a minimal BAM container so noodles can decode them.
    // We need to prepend a BAM header (magic + SAM text + ref dict) then the records.
    let mut bam_container: Vec<u8> = Vec::new();
    {
        let bgzf_writer = bgzf::io::Writer::new(&mut bam_container);
        let mut bam_writer = bam::io::Writer::from(bgzf_writer);
        bam_writer.write_header(header).map_err(RsomicsError::Io)?;
        bam_writer
            .get_mut()
            .write_all(&raw_buf)
            .map_err(RsomicsError::Io)?;
        // BGZF EOF block appended on drop.
    }

    let mut sam_writer = sam::io::Writer::new(&mut *out);
    let cursor = std::io::Cursor::new(bam_container);
    let mut bam_reader = bam::io::Reader::new(cursor);
    let _ = bam_reader.read_header().map_err(RsomicsError::Io)?;

    let mut record = bam::Record::default();
    loop {
        match bam_reader
            .read_record(&mut record)
            .map_err(RsomicsError::Io)?
        {
            0 => break,
            _ => {
                sam_writer
                    .write_alignment_record(header, &record)
                    .map_err(RsomicsError::Io)?;
            }
        }
    }
    out.flush().map_err(RsomicsError::Io)?;
    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reg2bin_basic() {
        // hts_reg2bin(0, 1, 14, 5) = 4681 (verified against samtools).
        assert_eq!(reg2bin(0, 1), 4681);
    }

    #[test]
    fn posmap_basic() {
        // ref_seq: [A, C, 0(gap), T] → posmap = [0, 1, 2, 2]
        let ref_seq = [1u8, 2, 0, 8];
        let posmap = build_posmap(&ref_seq);
        assert_eq!(posmap, [0, 1, 2, 2]);
    }

    #[test]
    fn remove_redundant_p_merges() {
        // 5M 2P 10M → 15M
        let mut cigar = vec![
            cigar_gen(5, CIGAR_M),
            cigar_gen(2, CIGAR_P),
            cigar_gen(10, CIGAR_M),
        ];
        remove_redundant_p(&mut cigar);
        assert_eq!(cigar, vec![cigar_gen(15, CIGAR_M)]);
    }

    #[test]
    fn remove_redundant_p_md_removes_p() {
        // 5M 2P 3D: P between M and D is removed but ops differ so no merge.
        let mut cigar = vec![
            cigar_gen(5, CIGAR_M),
            cigar_gen(2, CIGAR_P),
            cigar_gen(3, CIGAR_D),
        ];
        remove_redundant_p(&mut cigar);
        assert_eq!(cigar, vec![cigar_gen(5, CIGAR_M), cigar_gen(3, CIGAR_D)]);
    }

    #[test]
    fn cigar_gen_roundtrip() {
        for op in 0u8..9 {
            for len in [1u32, 100, 0xffff] {
                let enc = cigar_gen(len, op);
                assert_eq!((enc & 0xf) as u8, op);
                assert_eq!(enc >> 4, len);
            }
        }
    }
}
