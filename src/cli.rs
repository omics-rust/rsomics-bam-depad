use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_bam_depad::{DepadOpts, depad};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(
    name = "rsomics-bam-depad",
    version,
    about,
    long_about = None,
    disable_help_flag = true
)]
pub struct Cli {
    /// Input BAM file.
    pub input: PathBuf,

    /// Output SAM instead of BAM.
    #[arg(short = 's', long = "sam")]
    pub sam_output: bool,

    /// Output uncompressed BAM (cannot be combined with -s).
    #[arg(short = 'u', long = "uncompressed", conflicts_with = "sam_output")]
    pub uncompressed: bool,

    /// Fast (level 1) BAM compression (cannot be combined with -s).
    #[arg(
        short = '1',
        long = "fast-compression",
        conflicts_with = "sam_output",
        conflicts_with = "uncompressed"
    )]
    pub fast_compression: bool,

    /// Padded reference FASTA (used to correct @SQ lengths and as fallback when
    /// no embedded reference is present). The reference sequences must include
    /// the padding gap characters (`*` or `-`).
    #[arg(short = 'T', long = "reference")]
    pub reference: Option<PathBuf>,

    /// Output file path (default stdout).
    #[arg(short = 'o', long = "output", default_value = "-")]
    pub output: String,

    /// Do not add a @PG line to the output header.
    #[arg(long = "no-PG")]
    pub no_pg: bool,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let opts = DepadOpts {
            sam_output: self.sam_output,
            uncompressed: self.uncompressed,
            fast_compression: self.fast_compression,
            reference: self.reference,
            no_pg: self.no_pg,
        };

        let output_path = (self.output != "-").then(|| PathBuf::from(&self.output));
        let stats = depad(&self.input, output_path.as_deref(), &opts)?;

        if self.common.json {
            eprintln!(
                "{}",
                serde_json::to_string(&stats)
                    .map_err(|e| RsomicsError::InvalidInput(format!("JSON: {e}")))?
            );
        }

        Ok(())
    }
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }

    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        self.execute()
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Convert padded BAM to unpadded BAM (removes padding from MSA-derived alignments).",
    origin: Some(Origin {
        upstream: "samtools depad",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["<in.bam> [-s] [-T ref.fa] [-o out.bam]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('s'),
                long: "sam",
                aliases: &[],
                value: None,
                type_hint: None,
                required: false,
                default: Some("BAM"),
                description: "Write output as SAM instead of BAM.",
                why_default: None,
            },
            FlagSpec {
                short: Some('u'),
                long: "uncompressed",
                aliases: &[],
                value: None,
                type_hint: None,
                required: false,
                default: None,
                description: "Write uncompressed BAM output (incompatible with -s).",
                why_default: None,
            },
            FlagSpec {
                short: Some('1'),
                long: "fast-compression",
                aliases: &[],
                value: None,
                type_hint: None,
                required: false,
                default: None,
                description: "Write level-1 compressed BAM (incompatible with -s).",
                why_default: None,
            },
            FlagSpec {
                short: Some('T'),
                long: "reference",
                aliases: &[],
                value: Some("FILE"),
                type_hint: None,
                required: false,
                default: None,
                description: "Padded reference FASTA. Corrects @SQ LN fields and serves as fallback when no embedded reference is present in the BAM.",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("FILE"),
                type_hint: None,
                required: false,
                default: Some("stdout"),
                description: "Output file path.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "no-PG",
                aliases: &[],
                value: None,
                type_hint: None,
                required: false,
                default: None,
                description: "Do not add a @PG line to the output header.",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Convert padded BAM to unpadded BAM (requires embedded reference)",
            command: "rsomics-bam-depad in.padded.bam -o out.bam",
        },
        Example {
            description: "Use a padded FASTA reference (corrects @SQ lengths too)",
            command: "rsomics-bam-depad in.padded.bam -T padded_ref.fa -o out.bam",
        },
        Example {
            description: "Write SAM output",
            command: "rsomics-bam-depad in.padded.bam -s",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
