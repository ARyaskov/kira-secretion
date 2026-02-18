use std::io::Write;
use std::path::Path;

use serde::Serialize;
use thiserror::Error;

use crate::model::axes::{AxisConfig, AxisCoverage, AxisValues, saturating_map};
use crate::model::drivers::{format_drivers, format_eeb_drivers, top_k_eeb_drivers, top_k_panels};
use crate::pipeline::stage1_load::DatasetCtx;
use crate::pipeline::stage3_panels::{PanelCellPacked, PanelsContext};

#[derive(Debug, Error)]
pub enum Stage4Error {
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
}

#[derive(Debug, Clone, Serialize)]
pub struct AxisDrivers {
    pub sia: String,
    pub eeb: String,
    pub sli: String,
    pub mei: String,
    pub ecmi: String,
    pub apci: String,
    pub gdi: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct AxesContext {
    pub cell_ids: Vec<String>,
    pub values: Vec<AxisValues>,
    pub coverage: Vec<AxisCoverage>,
    pub drivers: Vec<AxisDrivers>,
    pub stats: AxesSummary,
}

#[derive(Debug, Clone, Serialize)]
pub struct AxisStats {
    pub median: f32,
    pub p90: f32,
    pub p99: f32,
    pub frac_ge_0_65: f32,
    pub frac_ge_0_80: f32,
}

#[derive(Debug, Clone, Serialize)]
pub struct AxisSummaryEntry {
    pub present: bool,
    pub value: AxisStats,
    pub coverage: AxisStats,
}

#[derive(Debug, Clone, Serialize)]
pub struct AxesSummary {
    pub sia: AxisSummaryEntry,
    pub eeb: AxisSummaryEntry,
    pub sli: AxisSummaryEntry,
    pub mei: AxisSummaryEntry,
    pub ecmi: AxisSummaryEntry,
    pub apci: AxisSummaryEntry,
    pub gdi: AxisSummaryEntry,
}

pub fn run_stage4_axes(
    _ctx: &DatasetCtx,
    panels_ctx: &PanelsContext,
    out_dir: &Path,
) -> Result<AxesContext, Stage4Error> {
    let cfg = AxisConfig::default();
    let indices = build_axis_indices(&panels_ctx.panels);

    let mut values = Vec::with_capacity(panels_ctx.cell_ids.len());
    let mut coverage = Vec::with_capacity(panels_ctx.cell_ids.len());
    let mut drivers = Vec::with_capacity(panels_ctx.cell_ids.len());

    let report_path = out_dir.join("axes.tsv");
    let mut writer = std::io::BufWriter::new(std::fs::File::create(&report_path)?);
    writer.write_all(b"cell_id\tSIA\tEEB\tSLI\tMEI\tECMI\tAPCI\tGDI\tcov_SIA\tcov_EEB\tcov_SLI\tcov_MEI\tcov_ECMI\tcov_APCI\tcov_GDI\tdrivers_SIA\tdrivers_EEB\tdrivers_SLI\tdrivers_MEI\tdrivers_ECMI\tdrivers_APCI\tdrivers_GDI\n")?;

    for (cell_idx, cell_id) in panels_ctx.cell_ids.iter().enumerate() {
        let packed = &panels_ctx.per_cell[cell_idx];
        let (vals, cov, drv) = compute_cell_axes(&indices, panels_ctx, packed, &cfg);

        let line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            cell_id,
            format_f32(vals.sia),
            format_f32(vals.eeb),
            format_f32(vals.sli),
            format_f32(vals.mei),
            format_f32(vals.ecmi),
            format_f32(vals.apci),
            format_f32(vals.gdi),
            format_f32(cov.sia),
            format_f32(cov.eeb),
            format_f32(cov.sli),
            format_f32(cov.mei),
            format_f32(cov.ecmi),
            format_f32(cov.apci),
            format_f32(cov.gdi),
            drv.sia,
            drv.eeb,
            drv.sli,
            drv.mei,
            drv.ecmi,
            drv.apci,
            drv.gdi
        );
        writer.write_all(line.as_bytes())?;

        values.push(vals);
        coverage.push(cov);
        drivers.push(drv);
    }

    writer.flush()?;

    let stats = compute_summary(&values, &coverage, &indices);

    Ok(AxesContext {
        cell_ids: panels_ctx.cell_ids.clone(),
        values,
        coverage,
        drivers,
        stats,
    })
}

fn compute_cell_axes(
    indices: &AxisIndices,
    panels_ctx: &PanelsContext,
    packed: &PanelCellPacked,
    cfg: &AxisConfig,
) -> (AxisValues, AxisCoverage, AxisDrivers) {
    let sia_raw = sum_panels(&indices.sia, packed);
    let sli_raw = sum_panels(&indices.sli, packed);
    let mei_raw = sum_panels(&indices.mei, packed);
    let ecmi_raw = sum_panels(&indices.ecmi, packed);
    let gdi_raw = sum_panels(&indices.gdi, packed);

    let export_raw = sum_panels(&indices.eeb_export, packed);
    let degrade_raw = sum_panels(&indices.eeb_degrade, packed);
    let denom = cfg.epsilon + export_raw + degrade_raw;
    let mut eeb = if denom > 0.0 {
        (export_raw - degrade_raw) / denom
    } else {
        0.0
    };
    if eeb > 1.0 {
        eeb = 1.0;
    } else if eeb < -1.0 {
        eeb = -1.0;
    }

    let apci_present = !indices.apci.is_empty();
    let apci_raw = if apci_present {
        sum_panels(&indices.apci, packed)
    } else {
        0.0
    };

    let sia = saturating_map(sia_raw, cfg.k);
    let sli = saturating_map(sli_raw, cfg.k);
    let mei = saturating_map(mei_raw, cfg.k);
    let ecmi = saturating_map(ecmi_raw, cfg.k);
    let gdi = saturating_map(gdi_raw, cfg.k);
    let apci = if apci_present {
        saturating_map(apci_raw, cfg.k)
    } else {
        f32::NAN
    };

    let cov_sia = coverage_axis(&indices.sia, panels_ctx, packed);
    let cov_sli = coverage_axis(&indices.sli, panels_ctx, packed);
    let cov_mei = coverage_axis(&indices.mei, panels_ctx, packed);
    let cov_ecmi = coverage_axis(&indices.ecmi, panels_ctx, packed);
    let cov_gdi = coverage_axis(&indices.gdi, panels_ctx, packed);
    let cov_eeb = coverage_axis_union(
        &indices.eeb_export,
        &indices.eeb_degrade,
        panels_ctx,
        packed,
    );
    let cov_apci = if apci_present {
        coverage_axis(&indices.apci, panels_ctx, packed)
    } else {
        0.0
    };

    let drivers_sia = drivers_for_axis(&indices.sia, panels_ctx, packed, 3);
    let drivers_sli = drivers_for_axis(&indices.sli, panels_ctx, packed, 3);
    let drivers_mei = drivers_for_axis(&indices.mei, panels_ctx, packed, 3);
    let drivers_ecmi = drivers_for_axis(&indices.ecmi, panels_ctx, packed, 3);
    let drivers_gdi = drivers_for_axis(&indices.gdi, panels_ctx, packed, 3);
    let drivers_apci = if apci_present {
        drivers_for_axis(&indices.apci, panels_ctx, packed, 3)
    } else {
        ".".to_string()
    };

    let drivers_eeb = drivers_for_eeb(
        &indices.eeb_export,
        &indices.eeb_degrade,
        panels_ctx,
        packed,
    );

    (
        AxisValues {
            sia,
            eeb,
            sli,
            mei,
            ecmi,
            apci,
            gdi,
        },
        AxisCoverage {
            sia: cov_sia,
            eeb: cov_eeb,
            sli: cov_sli,
            mei: cov_mei,
            ecmi: cov_ecmi,
            apci: cov_apci,
            gdi: cov_gdi,
        },
        AxisDrivers {
            sia: drivers_sia,
            eeb: drivers_eeb,
            sli: drivers_sli,
            mei: drivers_mei,
            ecmi: drivers_ecmi,
            apci: drivers_apci,
            gdi: drivers_gdi,
        },
    )
}

fn sum_panels(indices: &[usize], packed: &PanelCellPacked) -> f32 {
    let mut sum = 0.0;
    for idx in indices {
        sum += packed.sums[*idx];
    }
    sum
}

fn coverage_axis(indices: &[usize], panels_ctx: &PanelsContext, packed: &PanelCellPacked) -> f32 {
    let (required_total, required_missing) = coverage_counts(indices, panels_ctx, packed);
    if required_total == 0 {
        1.0
    } else {
        let cov = 1.0 - (required_missing as f32 / required_total as f32);
        cov.min(1.0).max(0.0)
    }
}

fn coverage_axis_union(
    export_idx: &[usize],
    degrade_idx: &[usize],
    panels_ctx: &PanelsContext,
    packed: &PanelCellPacked,
) -> f32 {
    let (total_a, missing_a) = coverage_counts(export_idx, panels_ctx, packed);
    let (total_b, missing_b) = coverage_counts(degrade_idx, panels_ctx, packed);
    let total = total_a + total_b;
    let missing = missing_a + missing_b;
    if total == 0 {
        1.0
    } else {
        let cov = 1.0 - (missing as f32 / total as f32);
        cov.min(1.0).max(0.0)
    }
}

fn coverage_counts(
    indices: &[usize],
    panels_ctx: &PanelsContext,
    packed: &PanelCellPacked,
) -> (u32, u32) {
    let mut total = 0u32;
    let mut missing = 0u32;
    for idx in indices {
        total += panels_ctx.mappings[*idx].required_total as u32;
        missing += packed.required_missing[*idx];
    }
    (total, missing)
}

fn drivers_for_axis(
    indices: &[usize],
    panels_ctx: &PanelsContext,
    packed: &PanelCellPacked,
    k: usize,
) -> String {
    if indices.is_empty() {
        return ".".to_string();
    }
    let mut ids = Vec::with_capacity(indices.len());
    let mut vals = Vec::with_capacity(indices.len());
    for idx in indices {
        ids.push(panels_ctx.panels.panels[*idx].id.clone());
        vals.push(packed.sums[*idx]);
    }
    let drivers = top_k_panels(&ids, &vals, k);
    format_drivers(&drivers)
}

fn drivers_for_eeb(
    export_idx: &[usize],
    degrade_idx: &[usize],
    panels_ctx: &PanelsContext,
    packed: &PanelCellPacked,
) -> String {
    let mut export_ids = Vec::with_capacity(export_idx.len());
    let mut export_vals = Vec::with_capacity(export_idx.len());
    for idx in export_idx {
        export_ids.push(panels_ctx.panels.panels[*idx].id.clone());
        export_vals.push(packed.sums[*idx]);
    }
    let mut degrade_ids = Vec::with_capacity(degrade_idx.len());
    let mut degrade_vals = Vec::with_capacity(degrade_idx.len());
    for idx in degrade_idx {
        degrade_ids.push(panels_ctx.panels.panels[*idx].id.clone());
        degrade_vals.push(packed.sums[*idx]);
    }

    let (export, degrade) =
        top_k_eeb_drivers(&export_ids, &export_vals, &degrade_ids, &degrade_vals, 2);
    format_eeb_drivers(&export, &degrade)
}

#[derive(Debug, Clone)]
struct AxisIndices {
    sia: Vec<usize>,
    eeb_export: Vec<usize>,
    eeb_degrade: Vec<usize>,
    sli: Vec<usize>,
    mei: Vec<usize>,
    ecmi: Vec<usize>,
    apci: Vec<usize>,
    gdi: Vec<usize>,
}

fn build_axis_indices(panels: &crate::panels::defs::PanelSet) -> AxisIndices {
    let mut indices = AxisIndices {
        sia: Vec::new(),
        eeb_export: Vec::new(),
        eeb_degrade: Vec::new(),
        sli: Vec::new(),
        mei: Vec::new(),
        ecmi: Vec::new(),
        apci: Vec::new(),
        gdi: Vec::new(),
    };

    for (idx, panel) in panels.panels.iter().enumerate() {
        match panel.axis.as_str() {
            "SIA" => indices.sia.push(idx),
            "EEB_EXPORT" => indices.eeb_export.push(idx),
            "EEB_DEGRADE" => indices.eeb_degrade.push(idx),
            "SLI" => indices.sli.push(idx),
            "MEI" => indices.mei.push(idx),
            "ECMI" => indices.ecmi.push(idx),
            "APCI" => indices.apci.push(idx),
            "GDI" => indices.gdi.push(idx),
            _ => {}
        }
    }

    indices
}

fn format_f32(value: f32) -> String {
    if value.is_nan() {
        "nan".to_string()
    } else {
        format!("{:.6}", value)
    }
}

fn compute_summary(
    values: &[AxisValues],
    coverage: &[AxisCoverage],
    indices: &AxisIndices,
) -> AxesSummary {
    AxesSummary {
        sia: summary_entry(
            values.iter().map(|v| v.sia),
            coverage.iter().map(|c| c.sia),
            true,
        ),
        eeb: summary_entry(
            values.iter().map(|v| v.eeb),
            coverage.iter().map(|c| c.eeb),
            true,
        ),
        sli: summary_entry(
            values.iter().map(|v| v.sli),
            coverage.iter().map(|c| c.sli),
            true,
        ),
        mei: summary_entry(
            values.iter().map(|v| v.mei),
            coverage.iter().map(|c| c.mei),
            true,
        ),
        ecmi: summary_entry(
            values.iter().map(|v| v.ecmi),
            coverage.iter().map(|c| c.ecmi),
            true,
        ),
        apci: summary_entry(
            values.iter().map(|v| v.apci),
            coverage.iter().map(|c| c.apci),
            !indices.apci.is_empty(),
        ),
        gdi: summary_entry(
            values.iter().map(|v| v.gdi),
            coverage.iter().map(|c| c.gdi),
            true,
        ),
    }
}

fn summary_entry<I1, I2>(values: I1, coverage: I2, present: bool) -> AxisSummaryEntry
where
    I1: Iterator<Item = f32>,
    I2: Iterator<Item = f32>,
{
    let mut vals: Vec<f32> = values.filter(|v| !v.is_nan()).collect();
    let mut covs: Vec<f32> = coverage.filter(|v| !v.is_nan()).collect();
    let value_stats = stats_from_vec(&mut vals);
    let coverage_stats = stats_from_vec(&mut covs);
    AxisSummaryEntry {
        present,
        value: value_stats,
        coverage: coverage_stats,
    }
}

fn stats_from_vec(values: &mut Vec<f32>) -> AxisStats {
    if values.is_empty() {
        return AxisStats {
            median: f32::NAN,
            p90: f32::NAN,
            p99: f32::NAN,
            frac_ge_0_65: 0.0,
            frac_ge_0_80: 0.0,
        };
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median = percentile(values, 0.5);
    let p90 = percentile(values, 0.9);
    let p99 = percentile(values, 0.99);
    let frac_ge_0_65 = fraction_ge(values, 0.65);
    let frac_ge_0_80 = fraction_ge(values, 0.80);
    AxisStats {
        median,
        p90,
        p99,
        frac_ge_0_65,
        frac_ge_0_80,
    }
}

fn percentile(values: &[f32], p: f32) -> f32 {
    if values.is_empty() {
        return f32::NAN;
    }
    let n = values.len();
    let idx = ((p * (n as f32 - 1.0)).floor() as usize).min(n - 1);
    values[idx]
}

fn fraction_ge(values: &[f32], threshold: f32) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    let mut count = 0usize;
    for v in values {
        if *v >= threshold {
            count += 1;
        }
    }
    count as f32 / values.len() as f32
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage4_axes.rs"]
mod tests;
