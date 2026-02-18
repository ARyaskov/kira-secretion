use std::path::PathBuf;
use std::time::Instant;

use clap::Args;
use tracing::info;

use crate::expr::normalize::Normalization;
use crate::panels::loader::{default_panels_dir, load_panels_from_dir};
use crate::pipeline::stage1_load::{DatasetCtx, RunMode, run_stage1};
use crate::pipeline::stage2_normalize::run_stage2;
use crate::pipeline::stage3_panels::run_stage3_panels;
use crate::pipeline::stage4_axes::run_stage4_axes;
use crate::pipeline::stage5_scores::run_stage5_scores;
use crate::pipeline::stage6_classify::run_stage6_classify;
use crate::pipeline::stage7_report::run_stage7_report;

#[derive(Args, Debug)]
pub struct RunArgs {
    /// Input 10x directory
    #[arg(long)]
    input: PathBuf,

    /// Output directory
    #[arg(long)]
    out: PathBuf,

    /// Mode for downstream processing (placeholder)
    #[arg(long, default_value = "cell")]
    mode: Mode,

    /// Optional metadata TSV
    #[arg(long)]
    meta: Option<PathBuf>,

    /// Input source mode
    #[arg(long, value_enum, default_value = "standalone")]
    pub(crate) run_mode: RunModeArg,

    /// Optional explicit shared cache path (kira-organelle.bin)
    #[arg(long)]
    cache: Option<PathBuf>,
}

#[derive(clap::ValueEnum, Clone, Debug)]
pub enum Mode {
    Cell,
    Sample,
}

#[derive(clap::ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub enum RunModeArg {
    Standalone,
    Pipeline,
}

impl From<RunModeArg> for RunMode {
    fn from(value: RunModeArg) -> Self {
        match value {
            RunModeArg::Standalone => RunMode::Standalone,
            RunModeArg::Pipeline => RunMode::Pipeline,
        }
    }
}

pub fn handle(args: RunArgs) -> anyhow::Result<()> {
    let stage_out = match args.run_mode {
        RunModeArg::Pipeline => args.out.join("kira-secretion"),
        RunModeArg::Standalone => args.out.clone(),
    };
    std::fs::create_dir_all(&stage_out)?;

    let start = Instant::now();
    info!(stage = "stage1_load", "starting stage");
    let ctx = run_stage1(
        &args.input,
        args.meta.as_deref(),
        &stage_out,
        true,
        args.run_mode.into(),
        args.cache.as_deref(),
    )?;
    info!(
        stage = "stage1_load",
        elapsed_ms = start.elapsed().as_millis(),
        "finished stage"
    );

    let start = Instant::now();
    info!(stage = "stage2_normalize", "starting stage");
    let expr_ctx = run_stage2(&ctx, &stage_out, Normalization::default(), true)?;
    info!(
        stage = "stage2_normalize",
        elapsed_ms = start.elapsed().as_millis(),
        nnz = expr_ctx.expr.nnz(),
        "finished stage"
    );

    write_expr_stats(&stage_out, &ctx, &expr_ctx.cell_stats)?;

    let start = Instant::now();
    info!(stage = "stage3_panels", "starting stage");
    let panels_dir = default_panels_dir();
    let panels = load_panels_from_dir(&panels_dir)?;
    if panels.panels.is_empty() {
        anyhow::bail!("no panels loaded");
    }
    let panels_ctx = run_stage3_panels(
        &expr_ctx,
        &panels,
        &ctx.gene_index,
        &ctx.barcodes,
        &stage_out,
    )?;
    let mapped_genes: usize = panels_ctx
        .mappings
        .iter()
        .map(|m| m.mapped.iter().filter(|v| v.is_some()).count())
        .sum();
    info!(
        stage = "stage3_panels",
        elapsed_ms = start.elapsed().as_millis(),
        panels = panels.panels.len(),
        genes = mapped_genes,
        "finished stage"
    );

    let start = Instant::now();
    info!(stage = "stage4_axes", "starting stage");
    let axes_ctx = run_stage4_axes(&ctx, &panels_ctx, &stage_out)?;
    let axis_counts = count_axis_panels(&panels_ctx);
    info!(
        stage = "stage4_axes",
        elapsed_ms = start.elapsed().as_millis(),
        sia = axis_counts.sia,
        eeb_export = axis_counts.eeb_export,
        eeb_degrade = axis_counts.eeb_degrade,
        sli = axis_counts.sli,
        mei = axis_counts.mei,
        ecmi = axis_counts.ecmi,
        apci = axis_counts.apci,
        gdi = axis_counts.gdi,
        "finished stage"
    );

    let start = Instant::now();
    info!(stage = "stage5_scores", "starting stage");
    let scores_ctx = run_stage5_scores(&axes_ctx, &stage_out)?;
    info!(
        stage = "stage5_scores",
        elapsed_ms = start.elapsed().as_millis(),
        "finished stage"
    );

    let start = Instant::now();
    info!(stage = "stage6_classify", "starting stage");
    let classify_ctx = run_stage6_classify(&ctx, &expr_ctx, &axes_ctx, &scores_ctx, &stage_out)?;
    log_regime_counts(&classify_ctx);
    info!(
        stage = "stage6_classify",
        elapsed_ms = start.elapsed().as_millis(),
        "finished stage"
    );

    let start = Instant::now();
    info!(stage = "stage7_report", "starting stage");
    let mode_str = match args.mode {
        Mode::Cell => "cell",
        Mode::Sample => "sample",
    };
    let _summary = run_stage7_report(
        &ctx,
        &expr_ctx,
        &axes_ctx,
        &scores_ctx,
        &classify_ctx,
        &panels_ctx,
        &stage_out,
        mode_str,
        args.run_mode.into(),
        args.meta.as_deref(),
    )?;
    info!(
        stage = "stage7_report",
        elapsed_ms = start.elapsed().as_millis(),
        "finished stage"
    );
    Ok(())
}

struct AxisCounts {
    sia: usize,
    eeb_export: usize,
    eeb_degrade: usize,
    sli: usize,
    mei: usize,
    ecmi: usize,
    apci: usize,
    gdi: usize,
}

fn count_axis_panels(panels_ctx: &crate::pipeline::stage3_panels::PanelsContext) -> AxisCounts {
    let mut counts = AxisCounts {
        sia: 0,
        eeb_export: 0,
        eeb_degrade: 0,
        sli: 0,
        mei: 0,
        ecmi: 0,
        apci: 0,
        gdi: 0,
    };
    for panel in &panels_ctx.panels.panels {
        match panel.axis.as_str() {
            "SIA" => counts.sia += 1,
            "EEB_EXPORT" => counts.eeb_export += 1,
            "EEB_DEGRADE" => counts.eeb_degrade += 1,
            "SLI" => counts.sli += 1,
            "MEI" => counts.mei += 1,
            "ECMI" => counts.ecmi += 1,
            "APCI" => counts.apci += 1,
            "GDI" => counts.gdi += 1,
            _ => {}
        }
    }
    counts
}

fn write_expr_stats(
    out_dir: &PathBuf,
    ctx: &DatasetCtx,
    cell_stats: &[crate::expr::csc::CellStats],
) -> anyhow::Result<()> {
    let path = out_dir.join("expr_stats.tsv");
    let mut buf = String::new();
    buf.push_str("cell_id\tlibsize\tdetected\n");
    for (barcode, stats) in ctx.barcodes.iter().zip(cell_stats.iter()) {
        buf.push_str(barcode);
        buf.push('\t');
        buf.push_str(&stats.libsize.to_string());
        buf.push('\t');
        buf.push_str(&stats.detected.to_string());
        buf.push('\n');
    }
    std::fs::write(path, buf)?;
    Ok(())
}

fn log_regime_counts(ctx: &crate::pipeline::stage6_classify::ClassifyContext) {
    for (regime, count) in &ctx.summary.counts {
        tracing::info!(regime = regime.as_str(), count = *count);
    }
}
