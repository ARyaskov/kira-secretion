use std::io::Write;
use std::path::Path;

use thiserror::Error;

use crate::model::drivers::top_k_components;
use crate::model::scores::{WeightsDefault, clamp01, pos_eeb};
use crate::pipeline::stage4_axes::AxesContext;

#[derive(Debug, Error)]
pub enum Stage5Error {
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct CompositeStats {
    pub median: f32,
    pub p90: f32,
    pub p99: f32,
    pub frac_ge_0_65: f32,
    pub frac_ge_0_80: f32,
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct CompositesSummary {
    pub oii: CompositeStats,
    pub iai: CompositeStats,
    pub esi: CompositeStats,
}

#[derive(Debug, Clone)]
pub struct ScoresContext {
    pub oii: Vec<f32>,
    pub iai: Vec<f32>,
    pub esi: Vec<f32>,
    pub cov_oii: Vec<f32>,
    pub cov_iai: Vec<f32>,
    pub cov_esi: Vec<f32>,
    pub drivers_oii: Vec<String>,
    pub drivers_iai: Vec<String>,
    pub drivers_esi: Vec<String>,
    pub summary: CompositesSummary,
}

pub fn run_stage5_scores(
    axes_ctx: &AxesContext,
    out_dir: &Path,
) -> Result<ScoresContext, Stage5Error> {
    let weights = WeightsDefault::default();

    let mut oii = Vec::with_capacity(axes_ctx.values.len());
    let mut iai = Vec::with_capacity(axes_ctx.values.len());
    let mut esi = Vec::with_capacity(axes_ctx.values.len());
    let mut cov_oii = Vec::with_capacity(axes_ctx.values.len());
    let mut cov_iai = Vec::with_capacity(axes_ctx.values.len());
    let mut cov_esi = Vec::with_capacity(axes_ctx.values.len());
    let mut drivers_oii = Vec::with_capacity(axes_ctx.values.len());
    let mut drivers_iai = Vec::with_capacity(axes_ctx.values.len());
    let mut drivers_esi = Vec::with_capacity(axes_ctx.values.len());

    let out_path = out_dir.join("composites.tsv");
    let mut writer = std::io::BufWriter::new(std::fs::File::create(&out_path)?);
    writer.write_all(b"cell_id\tOII\tIAI\tESI\tcov_OII\tcov_IAI\tcov_ESI\tdrivers_OII\tdrivers_IAI\tdrivers_ESI\n")?;

    for (idx, cell_id) in axes_ctx.cell_ids.iter().enumerate() {
        let v = &axes_ctx.values[idx];
        let cov = &axes_ctx.coverage[idx];

        let eeb_pos = pos_eeb(v.eeb);

        let oii_val = clamp01(
            weights.oii.sia * v.sia
                + weights.oii.pos_eeb * eeb_pos
                + weights.oii.sli * v.sli
                + weights.oii.mei * v.mei
                + weights.oii.ecmi * v.ecmi
                + weights.oii.gdi * v.gdi,
        );

        let (iai_val, iai_driver) = if v.apci.is_nan() {
            let val = clamp01(
                weights.iai_no_apci.mei * v.mei
                    + weights.iai_no_apci.gdi * v.gdi
                    + weights.iai_no_apci.sia * v.sia
                    + weights.iai_no_apci.pos_eeb * eeb_pos,
            );
            let names = ["MEI", "GDI", "SIA", "EEB_POS"];
            let contribs = [
                weights.iai_no_apci.mei * v.mei,
                weights.iai_no_apci.gdi * v.gdi,
                weights.iai_no_apci.sia * v.sia,
                weights.iai_no_apci.pos_eeb * eeb_pos,
            ];
            (val, top_k_components(&names, &contribs, 3))
        } else {
            let val = clamp01(
                weights.iai_with_apci.mei * v.mei
                    + weights.iai_with_apci.gdi * v.gdi
                    + weights.iai_with_apci.apci * v.apci
                    + weights.iai_with_apci.sia * v.sia
                    + weights.iai_with_apci.pos_eeb * eeb_pos,
            );
            let names = ["MEI", "GDI", "APCI", "SIA", "EEB_POS"];
            let contribs = [
                weights.iai_with_apci.mei * v.mei,
                weights.iai_with_apci.gdi * v.gdi,
                weights.iai_with_apci.apci * v.apci,
                weights.iai_with_apci.sia * v.sia,
                weights.iai_with_apci.pos_eeb * eeb_pos,
            ];
            (val, top_k_components(&names, &contribs, 3))
        };

        let esi_val = clamp01(
            weights.esi.ecmi * v.ecmi
                + weights.esi.mei * v.mei
                + weights.esi.pos_eeb * eeb_pos
                + weights.esi.sli * v.sli,
        );

        let oii_driver = {
            let names = ["SIA", "EEB_POS", "SLI", "MEI", "ECMI", "GDI"];
            let contribs = [
                weights.oii.sia * v.sia,
                weights.oii.pos_eeb * eeb_pos,
                weights.oii.sli * v.sli,
                weights.oii.mei * v.mei,
                weights.oii.ecmi * v.ecmi,
                weights.oii.gdi * v.gdi,
            ];
            top_k_components(&names, &contribs, 3)
        };
        let esi_driver = {
            let names = ["ECMI", "MEI", "EEB_POS", "SLI"];
            let contribs = [
                weights.esi.ecmi * v.ecmi,
                weights.esi.mei * v.mei,
                weights.esi.pos_eeb * eeb_pos,
                weights.esi.sli * v.sli,
            ];
            top_k_components(&names, &contribs, 3)
        };

        let cov_oii_val = weighted_cov_oii(cov, &weights);
        let cov_esi_val = weighted_cov_esi(cov, &weights);
        let cov_iai_val = if v.apci.is_nan() {
            weighted_cov_iai_no_apci(cov, &weights)
        } else {
            weighted_cov_iai(cov, &weights)
        };

        oii.push(oii_val);
        iai.push(iai_val);
        esi.push(esi_val);
        cov_oii.push(cov_oii_val);
        cov_iai.push(cov_iai_val);
        cov_esi.push(cov_esi_val);
        drivers_oii.push(oii_driver.clone());
        drivers_iai.push(iai_driver.clone());
        drivers_esi.push(esi_driver.clone());

        let line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            cell_id,
            format_f32(oii_val),
            format_f32(iai_val),
            format_f32(esi_val),
            format_f32(cov_oii_val),
            format_f32(cov_iai_val),
            format_f32(cov_esi_val),
            oii_driver,
            iai_driver,
            esi_driver
        );
        writer.write_all(line.as_bytes())?;
    }

    writer.flush()?;

    let summary = CompositesSummary {
        oii: summary_stats(&oii),
        iai: summary_stats(&iai),
        esi: summary_stats(&esi),
    };

    Ok(ScoresContext {
        oii,
        iai,
        esi,
        cov_oii,
        cov_iai,
        cov_esi,
        drivers_oii,
        drivers_iai,
        drivers_esi,
        summary,
    })
}

fn weighted_cov_oii(cov: &crate::model::axes::AxisCoverage, w: &WeightsDefault) -> f32 {
    let weights = [
        w.oii.sia,
        w.oii.pos_eeb,
        w.oii.sli,
        w.oii.mei,
        w.oii.ecmi,
        w.oii.gdi,
    ];
    let values = [cov.sia, cov.eeb, cov.sli, cov.mei, cov.ecmi, cov.gdi];
    weighted_cov(&weights, &values)
}

fn weighted_cov_esi(cov: &crate::model::axes::AxisCoverage, w: &WeightsDefault) -> f32 {
    let weights = [w.esi.ecmi, w.esi.mei, w.esi.pos_eeb, w.esi.sli];
    let values = [cov.ecmi, cov.mei, cov.eeb, cov.sli];
    weighted_cov(&weights, &values)
}

fn weighted_cov_iai(cov: &crate::model::axes::AxisCoverage, w: &WeightsDefault) -> f32 {
    let weights = [
        w.iai_with_apci.mei,
        w.iai_with_apci.gdi,
        w.iai_with_apci.apci,
        w.iai_with_apci.sia,
        w.iai_with_apci.pos_eeb,
    ];
    let values = [cov.mei, cov.gdi, cov.apci, cov.sia, cov.eeb];
    weighted_cov(&weights, &values)
}

fn weighted_cov_iai_no_apci(cov: &crate::model::axes::AxisCoverage, w: &WeightsDefault) -> f32 {
    let weights = [
        w.iai_no_apci.mei,
        w.iai_no_apci.gdi,
        w.iai_no_apci.sia,
        w.iai_no_apci.pos_eeb,
    ];
    let values = [cov.mei, cov.gdi, cov.sia, cov.eeb];
    weighted_cov(&weights, &values)
}

fn weighted_cov(weights: &[f32], values: &[f32]) -> f32 {
    let mut sum_w = 0.0;
    let mut sum = 0.0;
    for (w, v) in weights.iter().zip(values.iter()) {
        sum_w += *w;
        sum += *w * *v;
    }
    if sum_w == 0.0 {
        0.0
    } else {
        (sum / sum_w).min(1.0).max(0.0)
    }
}

fn summary_stats(values: &[f32]) -> CompositeStats {
    let mut vals: Vec<f32> = values.iter().copied().filter(|v| !v.is_nan()).collect();
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median = percentile(&vals, 0.5);
    let p90 = percentile(&vals, 0.9);
    let p99 = percentile(&vals, 0.99);
    let frac_ge_0_65 = fraction_ge(&vals, 0.65);
    let frac_ge_0_80 = fraction_ge(&vals, 0.80);
    CompositeStats {
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

fn format_f32(value: f32) -> String {
    if value.is_nan() {
        "nan".to_string()
    } else {
        format!("{:.6}", value)
    }
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage5_scores.rs"]
mod tests;
