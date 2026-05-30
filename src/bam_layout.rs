use std::io::Write;

use rsomics_common::{Result, RsomicsError};

// BAM CIGAR op codes (low nibble of the packed u32): M=0 I=1 D=2 N=3 S=4 H=5 P=6 ==7 X=8.
pub(crate) const CIGAR_M: u8 = 0;
pub(crate) const CIGAR_I: u8 = 1;
pub(crate) const CIGAR_D: u8 = 2;
pub(crate) const CIGAR_N: u8 = 3;
pub(crate) const CIGAR_S: u8 = 4;
pub(crate) const CIGAR_H: u8 = 5;
pub(crate) const CIGAR_P: u8 = 6;

// FLAG bits (SAMv1 §1.4).
pub(crate) const FLAG_UNMAPPED: u16 = 0x4;

// BAM payload layout constants (offsets from start of payload, after block_size).
pub(crate) const POS: usize = 4;
pub(crate) const L_READ_NAME: usize = 8;
pub(crate) const N_CIGAR: usize = 12;
pub(crate) const FLAG: usize = 14;
pub(crate) const L_SEQ: usize = 16;
pub(crate) const NEXT_REF_ID: usize = 20;
pub(crate) const NEXT_POS: usize = 24;
pub(crate) const FIXED_HEAD: usize = 32;

/// `hts_reg2bin(b, e, 14, 5)` — same as htslib `bam_reg2bin`.
pub(crate) fn reg2bin(beg: i64, end: i64) -> u16 {
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

/// Reference span of a CIGAR (M/D/N/=/X ops consume reference).
pub(crate) fn ref_span(cigar: &[u32]) -> i64 {
    cigar
        .iter()
        .filter(|&&op| matches!((op & 0xf) as u8, CIGAR_M | CIGAR_D | CIGAR_N | 7 | 8))
        .map(|&op| i64::from(op >> 4))
        .sum()
}

/// Pack (length, op) into the BAM u32: length in [31:4], op in [3:0].
#[inline]
pub(crate) fn cigar_gen(len: u32, op: u8) -> u32 {
    (len << 4) | u32::from(op)
}

/// Write a raw BAM payload with its 4-byte block_size prefix.
pub(crate) fn write_bytes<W: Write>(out: &mut W, bytes: &[u8]) -> Result<()> {
    let block_size = u32::try_from(bytes.len())
        .map_err(|e| RsomicsError::InvalidInput(format!("record too large: {e}")))?;
    out.write_all(&block_size.to_le_bytes())
        .map_err(RsomicsError::Io)?;
    out.write_all(bytes).map_err(RsomicsError::Io)?;
    Ok(())
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
