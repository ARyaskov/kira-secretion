use std::io::Write;
use std::path::Path;

use thiserror::Error;

use crate::model::flags::Flags;
use crate::model::regimes::{Regime, RuleId};
use crate::model::scores::pos_eeb;
use crate::model::thresholds::Thresholds;
use crate::pipeline::stage1_load::DatasetCtx;
use crate::pipeline::stage2_normalize::ExprContext;
use crate::pipeline::stage4_axes::AxesContext;
use crate::pipeline::stage5_scores::ScoresContext;

#[derive(Debug, Error)]
pub enum Stage6Error {
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
}

#[derive(Debug, Clone)]
pub struct ClassifyContext {
    pub regimes: Vec<Regime>,
    pub rule_ids: Vec<RuleId>,
    pub flags: Vec<Flags>,
    pub summary: RegimeSummary,
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct RegimeSummary {
    pub counts: Vec<(Regime, usize)>,
    pub fractions: Vec<(Regime, f32)>,
    pub flagged_fractions: Vec<(String, f32)>,
}

pub fn run_stage6_classify(
    dataset: &DatasetCtx,
    expr: &ExprContext,
    axes: &AxesContext,
    scores: &ScoresContext,
    out_dir: &Path,
) -> Result<ClassifyContext, Stage6Error> {
    let thresholds = Thresholds::default();
    let n = dataset.n_cells;

    let mut regimes = Vec::with_capacity(n);
    let mut rule_ids = Vec::with_capacity(n);
    let mut flags = Vec::with_capacity(n);

    let cell_ids = &dataset.barcodes;

    let out_path = out_dir.join("classify.tsv");
    let mut writer = std::io::BufWriter::new(std::fs::File::create(&out_path)?);
    writer.write_all(b"cell_id\tregime\trule_id\tflags\n")?;

    for idx in 0..n {
        let axis = &axes.values[idx];
        let cov = &axes.coverage[idx];
        let comp_oii = scores.oii[idx];
        let comp_esi = scores.esi[idx];

        let mut f = Flags::empty();
        let cell_stats = &expr.cell_stats[idx];
        if cell_stats.libsize < thresholds.low_counts as u64 {
            f.set(Flags::LOW_COUNTS);
        }
        if cell_stats.detected < thresholds.few_detected {
            f.set(Flags::FEW_DETECTED_GENES);
        }
        if cov.sia < thresholds.cov_min
            || cov.eeb < thresholds.cov_min
            || cov.sli < thresholds.cov_min
            || cov.mei < thresholds.cov_min
            || cov.ecmi < thresholds.cov_min
            || cov.gdi < thresholds.cov_min
            || (!axis.apci.is_nan() && cov.apci < thresholds.cov_min)
        {
            f.set(Flags::LOW_CONFIDENCE);
        }
        let eeb_pos = pos_eeb(axis.eeb);
        if f.contains(Flags::FEW_DETECTED_GENES)
            && axis.gdi >= thresholds.ambient_gdi
            && axis.sia < thresholds.ambient_sia
        {
            f.set(Flags::HIGH_AMBIENT_RISK);
        }

        let (regime, rule) = classify_cell(axis, eeb_pos, comp_oii, comp_esi, &thresholds);

        regimes.push(regime);
        rule_ids.push(rule);
        flags.push(f);

        let line = format!(
            "{}\t{}\t{}\t{}\n",
            cell_ids[idx],
            regime.as_str(),
            rule.as_str(),
            f.to_csv()
        );
        writer.write_all(line.as_bytes())?;
    }

    writer.flush()?;

    let summary = summarize(&regimes, &flags);

    Ok(ClassifyContext {
        regimes,
        rule_ids,
        flags,
        summary,
    })
}

fn classify_cell(
    axis: &crate::model::axes::AxisValues,
    pos_eeb: f32,
    oii: f32,
    esi: f32,
    t: &Thresholds,
) -> (Regime, RuleId) {
    if axis.sia < t.sia_low
        && pos_eeb < t.pos_eeb_low
        && axis.mei < 0.45
        && axis.ecmi < 0.45
        && axis.gdi < 0.50
    {
        return (Regime::SelfPreserving, RuleId::R1SelfPreserving);
    }
    if axis.sli >= t.sli_hi && axis.sia >= 0.45 {
        return (
            Regime::SecretoryLysosomeActive,
            RuleId::R2SecretoryLysosomeActive,
        );
    }
    if pos_eeb >= t.pos_eeb_hi && axis.sia >= t.sia_hi && oii >= 0.60 {
        return (Regime::ExportDominant, RuleId::R3ExportDominant);
    }
    if axis.mei >= t.mei_hi
        && (pos_eeb >= t.pos_eeb_mid || axis.sia >= t.sia_hi)
        && axis.gdi < t.gdi_hi
    {
        return (Regime::MetabolicSuppressive, RuleId::R4MetabolicSuppressive);
    }
    if axis.gdi >= t.gdi_hi && axis.sia >= t.sia_mid {
        return (Regime::InflammatorySignaler, RuleId::R5InflammatorySignaler);
    }
    if !axis.apci.is_nan() && axis.apci >= t.apci_hi && (axis.sia >= 0.45 || axis.gdi >= 0.60) {
        return (Regime::PresentationHigh, RuleId::R6PresentationHigh);
    }
    if (oii >= t.oii_hi && esi >= t.esi_hi) || esi >= t.esi_very {
        return (Regime::EnvironmentShaping, RuleId::R7EnvironmentShaping);
    }
    (Regime::Unclassified, RuleId::R0Unclassified)
}

fn summarize(regimes: &[Regime], flags: &[Flags]) -> RegimeSummary {
    let mut counts = Vec::new();
    let mut fractions = Vec::new();
    let n = regimes.len() as f32;

    for r in Regime::ordered() {
        let c = regimes.iter().filter(|v| **v == *r).count();
        counts.push((*r, c));
        let frac = if n == 0.0 { 0.0 } else { c as f32 / n };
        fractions.push((*r, frac));
    }

    let mut flagged = Vec::new();
    let flags_list = [
        ("LOW_CONFIDENCE", Flags::LOW_CONFIDENCE),
        ("FEW_DETECTED_GENES", Flags::FEW_DETECTED_GENES),
        ("LOW_COUNTS", Flags::LOW_COUNTS),
        ("HIGH_AMBIENT_RISK", Flags::HIGH_AMBIENT_RISK),
    ];
    for (name, bit) in flags_list {
        let c = flags.iter().filter(|f| f.contains(bit)).count();
        let frac = if n == 0.0 { 0.0 } else { c as f32 / n };
        flagged.push((name.to_string(), frac));
    }

    RegimeSummary {
        counts,
        fractions,
        flagged_fractions: flagged,
    }
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage6_classify.rs"]
mod tests;
