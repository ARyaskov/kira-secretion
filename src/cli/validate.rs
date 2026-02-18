use std::path::PathBuf;
use std::time::Instant;

use clap::Args;
use tracing::info;

use crate::pipeline::stage1_load::{DatasetCtx, RunMode, run_stage1};

#[derive(Args, Debug)]
pub struct ValidateArgs {
    /// Input 10x directory
    #[arg(long)]
    input: PathBuf,

    /// Output directory
    #[arg(long)]
    out: PathBuf,

    /// Optional metadata TSV
    #[arg(long)]
    meta: Option<PathBuf>,

    /// Skip full nnz line counting
    #[arg(long, default_value_t = true)]
    fast: bool,
}

pub fn handle(args: ValidateArgs) -> anyhow::Result<()> {
    std::fs::create_dir_all(&args.out)?;

    let start = Instant::now();
    info!(stage = "stage1_load", "starting stage");
    let ctx = run_stage1(
        &args.input,
        args.meta.as_deref(),
        &args.out,
        args.fast,
        RunMode::Standalone,
        None,
    )?;
    info!(
        stage = "stage1_load",
        elapsed_ms = start.elapsed().as_millis(),
        "finished stage"
    );

    write_validate(&args.out, &ctx)?;
    write_gene_warnings(&args.out, &ctx)?;
    Ok(())
}

fn write_validate(out_dir: &PathBuf, ctx: &DatasetCtx) -> anyhow::Result<()> {
    let mut lines = Vec::new();
    lines.push(("format", ctx.format.to_string()));
    lines.push(("n_genes", ctx.n_genes.to_string()));
    lines.push(("n_cells", ctx.n_cells.to_string()));
    lines.push(("nnz", ctx.nnz.to_string()));
    lines.push((
        "features_file",
        ctx.features_path.to_string_lossy().to_string(),
    ));
    lines.push((
        "barcodes_file",
        ctx.barcodes_path.to_string_lossy().to_string(),
    ));
    lines.push(("matrix_file", ctx.matrix_path.to_string_lossy().to_string()));
    lines.push(("meta_present", ctx.meta_present.to_string()));
    lines.push(("meta_cells_matched", ctx.meta_cells_matched.to_string()));
    lines.push(("meta_cells_missing", ctx.meta_cells_missing.to_string()));

    let path = out_dir.join("validate.tsv");
    let mut buf = String::new();
    for (k, v) in lines {
        buf.push_str(k);
        buf.push('\t');
        buf.push_str(&v);
        buf.push('\n');
    }
    std::fs::write(path, buf)?;
    Ok(())
}

fn write_gene_warnings(out_dir: &PathBuf, ctx: &DatasetCtx) -> anyhow::Result<()> {
    let path = out_dir.join("gene_mapping_warnings.tsv");
    let mut buf = String::new();
    buf.push_str("symbol\tfirst_row\tdup_row\n");
    for dup in &ctx.duplicate_gene_symbols {
        buf.push_str(&dup.symbol);
        buf.push('\t');
        buf.push_str(&dup.first_row.to_string());
        buf.push('\t');
        buf.push_str(&dup.dup_row.to_string());
        buf.push('\n');
    }
    std::fs::write(path, buf)?;
    Ok(())
}
