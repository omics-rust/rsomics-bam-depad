use std::path::Path;

use rsomics_common::{Result, RsomicsError};

use crate::bam_layout::{
    CIGAR_D, CIGAR_H, CIGAR_I, CIGAR_M, CIGAR_N, CIGAR_P, CIGAR_S, FIXED_HEAD, L_READ_NAME, N_CIGAR,
};

/// Build the padded→unpadded position map. `ref_seq[i] == 0` is a pad column.
pub(crate) fn build_posmap(ref_seq: &[u8]) -> Vec<i32> {
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

/// htslib `seq_nt16_table`: ASCII byte → 4-bit nucleotide code.
#[rustfmt::skip]
pub(crate) const SEQ_NT16_TABLE: [u8; 256] = [
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

/// Decode the padded reference sequence of a read into `out`.
///
/// Each element is 0 (D/gap) or a non-zero nibble (M/=/X base).
/// S ops advance the query index without output; H/P emit nothing.
/// N is treated as D (mirrors samtools). I is an error.
pub(crate) fn decode_seq(record: &[u8], out: &mut Vec<u8>, qname: &[u8]) -> Result<()> {
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
                    // seq_nt16 code 0 means '='; treat as non-zero (has a base).
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

/// Read a padded FASTA (gap chars `*` or `-`) for `ref_name` into a nibble array.
///
/// 0 = gap column, non-zero = real base (seq_nt16 code).
pub(crate) fn load_fasta_ref(
    fasta_path: &Path,
    ref_name: &str,
    padded_len: usize,
) -> Result<Vec<u8>> {
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
                break;
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
                    if c == 0 { 0xff } else { c }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn posmap_basic() {
        // [A, C, 0(gap), T] → posmap = [0, 1, 2, 2]
        let ref_seq = [1u8, 2, 0, 8];
        let posmap = build_posmap(&ref_seq);
        assert_eq!(posmap, [0, 1, 2, 2]);
    }
}
