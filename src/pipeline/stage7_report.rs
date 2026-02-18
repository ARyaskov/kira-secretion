use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt::Write as FmtWrite;
use std::io::{BufRead, BufWriter, Write};
use std::path::Path;

use serde::Serialize;
use serde_json::json;
use thiserror::Error;

use crate::input::open_reader;
use crate::model::flags::Flags;
use crate::model::regimes::Regime;
use crate::model::scores::pos_eeb;
use crate::pipeline::stage1_load::DatasetCtx;
use crate::pipeline::stage1_load::RunMode;
use crate::pipeline::stage2_normalize::ExprContext;
use crate::pipeline::stage3_panels::PanelsContext;
use crate::pipeline::stage4_axes::AxesContext;
use crate::pipeline::stage5_scores::ScoresContext;
use crate::pipeline::stage6_classify::ClassifyContext;
use crate::report::text::render_report;
use crate::simd;

#[derive(Debug, Error)]
pub enum Stage7Error {
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
    #[error("json error: {0}")]
    Json(#[from] serde_json::Error),
}

#[derive(Debug, Clone, Serialize)]
pub struct FinalSummary {
    pub tool: ToolSummary,
    pub input: InputSummary,
    pub distributions: DistributionSummary,
    pub regimes: RegimeSummary,
    pub qc: QcSummary,
}

#[derive(Debug, Clone, Serialize)]
pub struct ToolSummary {
    pub name: String,
    pub version: String,
    pub simd: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct InputSummary {
    pub n_cells: usize,
    pub species: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct DistributionSummary {
    pub secretory_load: Quantiles,
    pub er_golgi_pressure: Quantiles,
    pub stress_secretion_index: Quantiles,
}

#[derive(Debug, Clone, Serialize)]
pub struct Quantiles {
    pub median: f32,
    pub p90: f32,
    pub p99: f32,
}

#[derive(Debug, Clone, Serialize)]
pub struct RegimeSummary {
    pub counts: BTreeMap<String, usize>,
    pub fractions: BTreeMap<String, f32>,
}

#[derive(Debug, Clone, Serialize)]
pub struct QcSummary {
    pub low_confidence_fraction: f32,
    pub low_secretory_signal_fraction: f32,
}

#[derive(Debug, Clone)]
struct CellOutput {
    barcode: String,
    sample: String,
    condition: String,
    species: String,
    libsize: u64,
    nnz: u32,
    expressed_genes: u32,
    secretory_load: f32,
    exocytosis_bias: f32,
    vesicle_traffic_intensity: f32,
    er_golgi_pressure: f32,
    paracrine_signal_potential: f32,
    stress_secretion_index: f32,
    regime: String,
    flags: String,
    confidence: f32,
    low_confidence: bool,
    low_secretory_signal: bool,
}

#[derive(Debug, Clone, Default)]
struct MetaColumns {
    sample: Vec<String>,
    condition: Vec<String>,
    species: Vec<String>,
}

const PIPELINE_REGIMES: [&str; 6] = [
    "HomeostaticSecretion",
    "AdaptiveSecretion",
    "InflammatorySecretion",
    "HypersecretoryState",
    "SecretoryCollapse",
    "Unclassified",
];

pub fn run_stage7_report(
    dataset: &DatasetCtx,
    expr: &ExprContext,
    axes: &AxesContext,
    scores: &ScoresContext,
    classify: &ClassifyContext,
    panels: &PanelsContext,
    out_dir: &Path,
    _mode: &str,
    run_mode: RunMode,
    meta_path: Option<&Path>,
) -> Result<FinalSummary, Stage7Error> {
    std::fs::create_dir_all(out_dir)?;

    let meta = if let Some(path) = meta_path {
        read_meta_columns(path, &dataset.barcodes)?
    } else {
        MetaColumns {
            sample: vec![".".to_string(); dataset.n_cells],
            condition: vec![".".to_string(); dataset.n_cells],
            species: vec!["unknown".to_string(); dataset.n_cells],
        }
    };

    let mut rows = Vec::with_capacity(dataset.n_cells);
    for i in 0..dataset.n_cells {
        let axis = &axes.values[i];
        let cov = &axes.coverage[i];
        let exo_bias = clamp01(pos_eeb(axis.eeb));
        let secretory_load = clamp01(scores.oii[i]);
        let vesicle = clamp01(axis.sli);
        let er_golgi = clamp01(axis.sia);
        let paracrine = clamp01(scores.esi[i]);
        let stress = clamp01(axis.gdi);

        let confidence = clamp01(
            cov.sia
                .min(cov.eeb)
                .min(cov.sli)
                .min(cov.mei)
                .min(cov.ecmi)
                .min(cov.gdi)
                .min(scores.cov_oii[i])
                .min(scores.cov_esi[i]),
        );

        let regime = to_pipeline_regime(classify.regimes[i], secretory_load, stress, paracrine);

        let mut flag_set = Vec::new();
        let low_conf = classify.flags[i].contains(Flags::LOW_CONFIDENCE) || confidence < 0.60;
        let low_sig = secretory_load < 0.20 || vesicle < 0.20;
        if low_conf {
            flag_set.push("LOW_CONFIDENCE");
        }
        if low_sig {
            flag_set.push("LOW_SECRETORY_SIGNAL");
        }
        let flags = if flag_set.is_empty() {
            ".".to_string()
        } else {
            flag_set.join(",")
        };

        rows.push(CellOutput {
            barcode: dataset.barcodes[i].clone(),
            sample: meta.sample[i].clone(),
            condition: meta.condition[i].clone(),
            species: meta.species[i].clone(),
            libsize: expr.cell_stats[i].libsize,
            nnz: expr.cell_stats[i].detected,
            expressed_genes: expr.cell_stats[i].detected,
            secretory_load,
            exocytosis_bias: exo_bias,
            vesicle_traffic_intensity: vesicle,
            er_golgi_pressure: er_golgi,
            paracrine_signal_potential: paracrine,
            stress_secretion_index: stress,
            regime: regime.to_string(),
            flags,
            confidence,
            low_confidence: low_conf,
            low_secretory_signal: low_sig,
        });
    }

    let mut sorted_rows = rows.clone();
    sorted_rows.sort_by(|a, b| a.barcode.cmp(&b.barcode));
    write_secretion_tsv(out_dir, &sorted_rows)?;
    write_panels_report(out_dir, panels)?;

    let summary = build_summary(&rows);
    write_summary_json(out_dir, &summary)?;
    if run_mode == RunMode::Pipeline {
        write_pipeline_step_json(out_dir)?;
    }

    std::fs::write(out_dir.join("report.txt"), render_report(&summary))?;

    Ok(summary)
}

fn write_secretion_tsv(out_dir: &Path, rows: &[CellOutput]) -> Result<(), Stage7Error> {
    let mut writer = BufWriter::new(std::fs::File::create(out_dir.join("secretion.tsv"))?);
    writer.write_all(b"barcode\tsample\tcondition\tspecies\tlibsize\tnnz\texpressed_genes\tsecretory_load\texocytosis_bias\tvesicle_traffic_intensity\ter_golgi_pressure\tparacrine_signal_potential\tstress_secretion_index\tregime\tflags\tconfidence\n")?;

    for row in rows {
        let line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            row.barcode,
            row.sample,
            row.condition,
            row.species,
            row.libsize,
            row.nnz,
            row.expressed_genes,
            fmt6(row.secretory_load),
            fmt6(row.exocytosis_bias),
            fmt6(row.vesicle_traffic_intensity),
            fmt6(row.er_golgi_pressure),
            fmt6(row.paracrine_signal_potential),
            fmt6(row.stress_secretion_index),
            row.regime,
            row.flags,
            fmt6(row.confidence),
        );
        writer.write_all(line.as_bytes())?;
    }
    writer.flush()?;
    Ok(())
}

fn write_summary_json(out_dir: &Path, summary: &FinalSummary) -> Result<(), Stage7Error> {
    fn push_quoted(buf: &mut String, s: &str) -> Result<(), Stage7Error> {
        buf.push_str(&serde_json::to_string(s)?);
        Ok(())
    }

    let mut out = String::with_capacity(2048);
    out.push_str("{\n");
    out.push_str("  \"tool\": {\n");
    out.push_str("    \"name\": ");
    push_quoted(&mut out, &summary.tool.name)?;
    out.push_str(",\n");
    out.push_str("    \"version\": ");
    push_quoted(&mut out, &summary.tool.version)?;
    out.push_str(",\n");
    out.push_str("    \"simd\": ");
    push_quoted(&mut out, &summary.tool.simd)?;
    out.push_str("\n");
    out.push_str("  },\n");
    out.push_str("  \"input\": {\n");
    let _ = write!(out, "    \"n_cells\": {},\n", summary.input.n_cells);
    out.push_str("    \"species\": ");
    push_quoted(&mut out, &summary.input.species)?;
    out.push_str("\n");
    out.push_str("  },\n");
    out.push_str("  \"distributions\": {\n");
    out.push_str("    \"secretory_load\": {");
    push_quantiles_json(&mut out, &summary.distributions.secretory_load);
    out.push_str("},\n");
    out.push_str("    \"er_golgi_pressure\": {");
    push_quantiles_json(&mut out, &summary.distributions.er_golgi_pressure);
    out.push_str("},\n");
    out.push_str("    \"stress_secretion_index\": {");
    push_quantiles_json(&mut out, &summary.distributions.stress_secretion_index);
    out.push_str("}\n");
    out.push_str("  },\n");
    out.push_str("  \"regimes\": {\n");
    out.push_str("    \"counts\": {\n");
    let mut counts_iter = summary.regimes.counts.iter().peekable();
    while let Some((name, count)) = counts_iter.next() {
        out.push_str("      ");
        push_quoted(&mut out, name)?;
        let _ = write!(out, ": {}", count);
        if counts_iter.peek().is_some() {
            out.push(',');
        }
        out.push('\n');
    }
    out.push_str("    },\n");
    out.push_str("    \"fractions\": {\n");
    let mut fracs_iter = summary.regimes.fractions.iter().peekable();
    while let Some((name, frac)) = fracs_iter.next() {
        out.push_str("      ");
        push_quoted(&mut out, name)?;
        let _ = write!(out, ": {}", fmt6(*frac));
        if fracs_iter.peek().is_some() {
            out.push(',');
        }
        out.push('\n');
    }
    out.push_str("    }\n");
    out.push_str("  },\n");
    out.push_str("  \"qc\": {\n");
    let _ = write!(
        out,
        "    \"low_confidence_fraction\": {},\n",
        fmt6(summary.qc.low_confidence_fraction)
    );
    let _ = write!(
        out,
        "    \"low_secretory_signal_fraction\": {}\n",
        fmt6(summary.qc.low_secretory_signal_fraction)
    );
    out.push_str("  }\n");
    out.push_str("}\n");
    std::fs::write(out_dir.join("summary.json"), out)?;
    Ok(())
}

fn push_quantiles_json(buf: &mut String, q: &Quantiles) {
    let _ = write!(
        buf,
        "\"median\": {}, \"p90\": {}, \"p99\": {}",
        fmt6(q.median),
        fmt6(q.p90),
        fmt6(q.p99),
    );
}

fn write_pipeline_step_json(out_dir: &Path) -> Result<(), Stage7Error> {
    let pipeline_step = json!({
        "tool": {
            "name": "kira-secretion",
            "stage": "secretion",
            "version": env!("CARGO_PKG_VERSION")
        },
        "artifacts": {
            "summary": "summary.json",
            "primary_metrics": "secretion.tsv",
            "panels": "panels_report.tsv"
        },
        "cell_metrics": {
            "file": "secretion.tsv",
            "id_column": "barcode",
            "regime_column": "regime",
            "confidence_column": "confidence",
            "flag_column": "flags"
        },
        "regimes": PIPELINE_REGIMES
    });
    std::fs::write(
        out_dir.join("pipeline_step.json"),
        serde_json::to_string_pretty(&pipeline_step)?,
    )?;
    Ok(())
}

fn write_panels_report(out_dir: &Path, panels: &PanelsContext) -> Result<(), Stage7Error> {
    let mut writer = BufWriter::new(std::fs::File::create(out_dir.join("panels_report.tsv"))?);
    writer.write_all(b"panel_id\tpanel_name\tpanel_group\tpanel_size_defined\tpanel_size_mappable\tmissing_genes\tcoverage_median\tcoverage_p10\tsum_median\tsum_p90\tsum_p99\n")?;

    for (panel_idx, panel) in panels.panels.panels.iter().enumerate() {
        let mapping = &panels.mappings[panel_idx];
        let mut coverages = Vec::with_capacity(panels.per_cell.len());
        let mut sums = Vec::with_capacity(panels.per_cell.len());

        for cell in &panels.per_cell {
            sums.push(cell.sums[panel_idx]);
            let req_total = mapping.required_total as u32;
            let missing = cell.required_missing[panel_idx];
            let cov = if req_total == 0 {
                1.0
            } else {
                1.0 - (missing as f32 / req_total as f32)
            };
            coverages.push(clamp01(cov));
        }

        let mut missing = Vec::new();
        for (gene_pos, mapped) in mapping.mapped.iter().enumerate() {
            if mapped.is_none() {
                missing.push(panel.genes[gene_pos].symbol.clone());
            }
        }

        coverages.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        sums.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            panel.id,
            panel.description,
            panel.axis,
            panel.genes.len(),
            mapping.mapped.iter().filter(|m| m.is_some()).count(),
            if missing.is_empty() {
                ".".to_string()
            } else {
                missing.join(",")
            },
            fmt6(percentile(&coverages, 0.5)),
            fmt6(percentile(&coverages, 0.10)),
            fmt6(percentile(&sums, 0.5)),
            fmt6(percentile(&sums, 0.90)),
            fmt6(percentile(&sums, 0.99)),
        );
        writer.write_all(line.as_bytes())?;
    }

    writer.flush()?;
    Ok(())
}

fn read_meta_columns(path: &Path, barcodes: &[String]) -> Result<MetaColumns, Stage7Error> {
    let mut sample = vec![".".to_string(); barcodes.len()];
    let mut condition = vec![".".to_string(); barcodes.len()];
    let mut species = vec!["unknown".to_string(); barcodes.len()];

    let mut index: HashMap<&str, usize> = HashMap::new();
    for (i, bc) in barcodes.iter().enumerate() {
        index.insert(bc.as_str(), i);
    }

    let mut seen: HashSet<String> = HashSet::new();
    let mut reader = open_reader(path).map_err(|e| std::io::Error::other(e.to_string()))?;

    let mut header = String::new();
    if reader.read_line(&mut header)? == 0 {
        return Ok(MetaColumns {
            sample,
            condition,
            species,
        });
    }

    let cols: Vec<&str> = header.trim_end_matches(['\n', '\r']).split('\t').collect();
    let cell_idx = cols.iter().position(|c| *c == "cell_id");
    let sample_idx = cols.iter().position(|c| *c == "sample_id");
    let cond_idx = cols.iter().position(|c| *c == "condition");
    let species_idx = cols.iter().position(|c| *c == "species");

    let Some(cell_col) = cell_idx else {
        return Ok(MetaColumns {
            sample,
            condition,
            species,
        });
    };

    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        let raw = line.trim_end_matches(['\n', '\r']);
        if raw.is_empty() {
            continue;
        }
        let parts: Vec<&str> = raw.split('\t').collect();
        if cell_col >= parts.len() {
            continue;
        }
        let cell = parts[cell_col];
        if cell.is_empty() || !seen.insert(cell.to_string()) {
            continue;
        }
        let Some(&i) = index.get(cell) else {
            continue;
        };

        if let Some(idx) = sample_idx {
            if idx < parts.len() && !parts[idx].is_empty() {
                sample[i] = parts[idx].to_string();
            }
        }
        if let Some(idx) = cond_idx {
            if idx < parts.len() && !parts[idx].is_empty() {
                condition[i] = parts[idx].to_string();
            }
        }
        if let Some(idx) = species_idx {
            if idx < parts.len() && !parts[idx].is_empty() {
                species[i] = normalize_species(parts[idx]);
            }
        }
    }

    Ok(MetaColumns {
        sample,
        condition,
        species,
    })
}

fn normalize_species(s: &str) -> String {
    let x = s.trim().to_ascii_lowercase();
    if x.contains("human") || x == "hs" || x == "homo_sapiens" {
        "human".to_string()
    } else if x.contains("mouse") || x == "mm" || x == "mus_musculus" {
        "mouse".to_string()
    } else {
        "unknown".to_string()
    }
}

fn build_summary(rows: &[CellOutput]) -> FinalSummary {
    let species = rows
        .iter()
        .find(|r| r.species == "human" || r.species == "mouse")
        .map(|r| r.species.clone())
        .unwrap_or_else(|| "unknown".to_string());

    let secretory: Vec<f32> = rows.iter().map(|r| r.secretory_load).collect();
    let er_golgi: Vec<f32> = rows.iter().map(|r| r.er_golgi_pressure).collect();
    let stress: Vec<f32> = rows.iter().map(|r| r.stress_secretion_index).collect();

    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    for name in PIPELINE_REGIMES {
        counts.insert(name.to_string(), 0);
    }
    for row in rows {
        if let Some(c) = counts.get_mut(&row.regime) {
            *c += 1;
        }
    }

    let n = rows.len() as f32;
    let mut fracs = BTreeMap::new();
    for (name, count) in &counts {
        fracs.insert(name.clone(), if n == 0.0 { 0.0 } else { *count as f32 / n });
    }

    let low_conf_count = rows.iter().filter(|r| r.low_confidence).count() as f32;
    let low_sig_count = rows.iter().filter(|r| r.low_secretory_signal).count() as f32;

    FinalSummary {
        tool: ToolSummary {
            name: "kira-secretion".to_string(),
            version: env!("CARGO_PKG_VERSION").to_string(),
            simd: simd_name(),
        },
        input: InputSummary {
            n_cells: rows.len(),
            species,
        },
        distributions: DistributionSummary {
            secretory_load: stats(&secretory),
            er_golgi_pressure: stats(&er_golgi),
            stress_secretion_index: stats(&stress),
        },
        regimes: RegimeSummary {
            counts,
            fractions: fracs,
        },
        qc: QcSummary {
            low_confidence_fraction: if n == 0.0 { 0.0 } else { low_conf_count / n },
            low_secretory_signal_fraction: if n == 0.0 { 0.0 } else { low_sig_count / n },
        },
    }
}

fn simd_name() -> String {
    simd::backend_name().to_string()
}

fn stats(values: &[f32]) -> Quantiles {
    let mut vals: Vec<f32> = values.iter().copied().filter(|v| v.is_finite()).collect();
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    Quantiles {
        median: percentile(&vals, 0.5),
        p90: percentile(&vals, 0.9),
        p99: percentile(&vals, 0.99),
    }
}

fn percentile(values: &[f32], p: f32) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    let idx = ((p * (values.len() as f32 - 1.0)).floor() as usize).min(values.len() - 1);
    values[idx]
}

fn to_pipeline_regime(
    old: Regime,
    secretory_load: f32,
    stress: f32,
    paracrine: f32,
) -> &'static str {
    if secretory_load < 0.20 {
        return "SecretoryCollapse";
    }
    if secretory_load >= 0.80 && stress >= 0.75 {
        return "HypersecretoryState";
    }
    if stress >= 0.75 {
        return "InflammatorySecretion";
    }

    match old {
        Regime::SelfPreserving => "HomeostaticSecretion",
        Regime::InflammatorySignaler => "InflammatorySecretion",
        Regime::MetabolicSuppressive => "SecretoryCollapse",
        Regime::Unclassified => {
            if paracrine >= 0.65 {
                "AdaptiveSecretion"
            } else {
                "Unclassified"
            }
        }
        _ => "AdaptiveSecretion",
    }
}

fn fmt6(v: f32) -> String {
    if v.is_finite() {
        format!("{:.6}", clamp01(v))
    } else {
        "0.000000".to_string()
    }
}

fn clamp01(v: f32) -> f32 {
    v.clamp(0.0, 1.0)
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage7_report.rs"]
mod tests;
