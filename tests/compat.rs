//! Field-level compat against `samtools depad`.
//!
//! The test converts a padded BAM (with an embedded reference) through both
//! samtools and our binary, then diffs the SAM output record-by-record
//! (FLAG / POS / CIGAR / MAPQ / SEQ / QUAL / aux all match; the @PG lines
//! differ and are excluded from the comparison, as is TLEN which samtools
//! doesn't recompute).
//!
//! Skips gracefully when samtools is absent or below 1.3 (when `depad` was
//! introduced).

use std::path::{Path, PathBuf};
use std::process::Command;

fn ours() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam-depad"))
}

fn golden(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn samtools_version() -> Option<(u32, u32)> {
    let out = Command::new("samtools").arg("--version").output().ok()?;
    if !out.status.success() {
        return None;
    }
    let stdout = String::from_utf8_lossy(&out.stdout);
    let ver_str = stdout
        .lines()
        .next()?
        .split_whitespace()
        .nth(1)
        .unwrap_or("");
    let mut it = ver_str.split('.');
    let major: u32 = it.next()?.parse().ok()?;
    let minor: u32 = it.next().and_then(|s| s.parse().ok()).unwrap_or(0);
    Some((major, minor))
}

fn samtools_ready() -> bool {
    match samtools_version() {
        Some((major, minor)) if major > 1 || (major == 1 && minor >= 3) => true,
        Some((major, minor)) => {
            eprintln!("SKIP depad compat: samtools {major}.{minor} (need >= 1.3 for `depad`)");
            false
        }
        None => {
            eprintln!("SKIP depad compat: samtools not found");
            false
        }
    }
}

/// Run `samtools depad -s <bam>` and return the stdout, minus @PG lines.
fn samtools_depad_sam(bam: &Path) -> Vec<u8> {
    let out = Command::new("samtools")
        .args(["depad", "-s", "--no-PG"])
        .arg(bam)
        .output()
        .unwrap();
    assert!(
        out.status.success() || {
            // samtools exits 1 even on success when printing warnings to stderr.
            let err = String::from_utf8_lossy(&out.stderr);
            err.contains("Warning") && !err.contains("ERROR")
        },
        "samtools depad failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    filter_pg_lines(&out.stdout)
}

/// Run `rsomics-bam-depad -s <bam>` and return the stdout, minus @PG lines.
fn our_depad_sam(bam: &Path) -> Vec<u8> {
    let out = ours().args(["-s", "--no-PG"]).arg(bam).output().unwrap();
    assert!(
        out.status.success(),
        "rsomics-bam-depad failed:\nstderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    filter_pg_lines(&out.stdout)
}

fn filter_pg_lines(sam: &[u8]) -> Vec<u8> {
    sam.split(|&b| b == b'\n')
        .filter(|line| !line.starts_with(b"@PG"))
        .flat_map(|line| line.iter().copied().chain(std::iter::once(b'\n')))
        .collect()
}

#[test]
fn depad_matches_samtools_golden() {
    if !samtools_ready() {
        return;
    }
    let bam = golden("padded.bam");
    let expected = samtools_depad_sam(&bam);
    let got = our_depad_sam(&bam);

    if expected != got {
        let exp_str = String::from_utf8_lossy(&expected);
        let got_str = String::from_utf8_lossy(&got);
        panic!(
            "depad output differs from samtools on golden fixture.\n\
             === samtools ===\n{exp_str}\n=== ours ===\n{got_str}"
        );
    }
}

/// Verify BAM-mode output also matches by decoding both to SAM.
#[test]
fn depad_bam_mode_matches_samtools() {
    if !samtools_ready() {
        return;
    }
    let bam = golden("padded.bam");

    // samtools depad → BAM → SAM
    let samtools_bam = {
        let dir = std::env::temp_dir().join("rsomics-bam-depad-compat");
        let _ = std::fs::create_dir_all(&dir);
        let out_bam = dir.join("samtools_out.bam");
        let status = Command::new("samtools")
            .args(["depad", "--no-PG", "-o"])
            .arg(&out_bam)
            .arg(&bam)
            .status()
            .unwrap();
        assert!(status.success(), "samtools depad (BAM) failed");
        let view = Command::new("samtools")
            .args(["view", "-h"])
            .arg(&out_bam)
            .output()
            .unwrap();
        filter_pg_lines(&view.stdout)
    };

    // ours depad → BAM → SAM
    let our_bam = {
        let dir = std::env::temp_dir().join("rsomics-bam-depad-compat");
        let out_bam = dir.join("ours_out.bam");
        let out = ours()
            .args(["--no-PG", "-o"])
            .arg(&out_bam)
            .arg(&bam)
            .output()
            .unwrap();
        assert!(
            out.status.success(),
            "rsomics-bam-depad (BAM) failed: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let view = Command::new("samtools")
            .args(["view", "-h"])
            .arg(&out_bam)
            .output()
            .unwrap();
        filter_pg_lines(&view.stdout)
    };

    if samtools_bam != our_bam {
        let exp_str = String::from_utf8_lossy(&samtools_bam);
        let got_str = String::from_utf8_lossy(&our_bam);
        panic!(
            "BAM-mode depad output differs from samtools.\n\
             === samtools ===\n{exp_str}\n=== ours ===\n{got_str}"
        );
    }
}
