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

use crate::bam_layout::{
    CIGAR_M, FIXED_HEAD, FLAG, FLAG_UNMAPPED, L_READ_NAME, L_SEQ, N_CIGAR, NEXT_POS, NEXT_REF_ID,
    POS, cigar_gen, ref_span, reg2bin, write_bytes,
};
use crate::cigar::{rebuild_cigar, replace_cigar};
use crate::posmap::{build_posmap, decode_seq, load_fasta_ref};

/// Depad options passed from `cli.rs`.
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
    {
        let mut sam_writer = noodles::sam::io::Writer::new(&mut *writer);
        sam_writer.write_header(&header).map_err(RsomicsError::Io)?;
    }
    // SAM output: correct but slower — collect via BAM buffer, then decode.
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

/// SAM-output path: depad via `run_core`, then decode raw BAM bytes to SAM text.
///
/// Used only for `-s`; correct over fast.
fn run_core_sam<W: Write>(
    reader: &mut rsomics_bamio::ParallelBamReader,
    out: &mut W,
    header: &Header,
    ref_names: &[String],
    ref_lengths: &[usize],
    opts: &DepadOpts,
) -> Result<DepadStats> {
    let mut raw_buf: Vec<u8> = Vec::new();
    let stats = run_core(reader, &mut raw_buf, header, ref_names, ref_lengths, opts)?;

    // Wrap raw record bytes in a minimal BAM container so noodles can decode them.
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
