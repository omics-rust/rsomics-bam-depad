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

mod bam_layout;
mod cigar;
mod driver;
mod posmap;

pub use driver::{DepadOpts, DepadStats, depad};
