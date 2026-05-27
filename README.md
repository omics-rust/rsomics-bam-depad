# rsomics-bam-depad

Convert padded BAM alignments to unpadded BAM — Rust port of `samtools depad`.

```sh
rsomics-bam-depad in.padded.bam -o out.bam          # embedded reference in BAM
rsomics-bam-depad in.padded.bam -T padded_ref.fa    # FASTA reference (corrects @SQ LN)
rsomics-bam-depad in.padded.bam -s                   # SAM output
rsomics-bam-depad in.padded.bam -u -o out.bam        # uncompressed BAM
```

## What it does

Padded BAM arises from multiple-sequence assembly pipelines (ACE/CAF output,
padded-reference aligners). The padded reference has gap columns (`*`/`-`); the
BAM header `@SQ LN` fields and read `POS` values are in padded coordinates.
`rsomics-bam-depad` converts to unpadded coordinates so the output can be used
with a standard linear reference.

The algorithm:

1. For each reference, build a **position map** from an embedded reference record
   (a read whose `QNAME == @SQ SN` at `POS = 0`) or a padded FASTA (`-T`).
2. Reconstruct each read's CIGAR by comparing query vs. padded reference
   base-by-base (M / I / D / P depending on whether each position has a real
   base in query and/or reference).
3. Remove redundant `P` operators that fall between `M`/`D` operators
   (e.g. `5M2P10M` → `15M`). Leading `P` ops are preserved.
4. Update `POS` and `PNEXT` through the position map; recalculate BAM `BIN`.

Input reads must not contain `P` or `I` CIGAR operators (same constraint as
`samtools depad`).

## Options

| Flag | Meaning |
|---|---|
| `-s, --sam` | Write SAM instead of BAM (default: BAM). |
| `-u, --uncompressed` | Write uncompressed BAM. |
| `-1, --fast-compression` | Write level-1 compressed BAM. |
| `-T, --reference FILE` | Padded reference FASTA. Corrects `@SQ LN` fields and serves as fallback when no embedded reference is present. |
| `-o, --output FILE` | Output file path (default stdout). |
| `--no-PG` | Do not add a `@PG` line. |

## Origin

This crate is a Rust reimplementation of `samtools depad`, informed by the
upstream MIT-licensed source (`padding.c:bam_pad2unpad`, `samtools` 1.23.1).
The implementation follows the source algorithm directly: embedded-reference
detection by `QNAME == @SQ SN` at `POS == 0`, the per-base CIGAR reconstruction
loop (query × reference nibble comparison → M/I/D/P), the redundant-P removal
pass, position-map construction and application, and BIN recalculation via
`hts_reg2bin`. The FASTA-reference path (`-T`) reads padded gap characters
(`*` or `-`) and maps them to zero the same way `load_unpadded_ref` does.

License: MIT OR Apache-2.0.
Upstream credit: [samtools](https://github.com/samtools/samtools) (MIT/Expat).
