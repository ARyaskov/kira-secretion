use std::io::Write;
use std::path::Path;

use thiserror::Error;

use crate::expr::csc::CellStats;
use crate::input::InputError;
use crate::input::features::GeneIndex;
use crate::panels::defs::PanelSet;
use crate::panels::mapping::{GeneMapping, MappingWarning, map_panel};
use crate::pipeline::stage2_normalize::ExprContext;

#[derive(Debug, Error)]
pub enum Stage3Error {
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
    #[error("input error: {0}")]
    Input(#[from] InputError),
}

#[derive(Debug, Clone)]
pub struct PanelAccum {
    pub sum: f32,
    pub hits: u32,
}

#[derive(Debug, Clone)]
pub struct PanelCellResult {
    pub panel_id: String,
    pub sum: f32,
    pub hits: u32,
    pub coverage: f32,
    pub required_missing: u32,
}

#[derive(Debug, Clone)]
pub struct PanelCellPacked {
    pub sums: Vec<f32>,
    pub hits: Vec<u32>,
    pub required_missing: Vec<u32>,
}

#[derive(Debug, Clone)]
pub struct PanelsContext {
    pub panels: PanelSet,
    pub mappings: Vec<GeneMapping>,
    pub warnings: Vec<MappingWarning>,
    pub cell_ids: Vec<String>,
    pub per_cell: Vec<PanelCellPacked>,
}

pub fn run_stage3_panels(
    expr: &ExprContext,
    panels: &PanelSet,
    gene_index: &GeneIndex,
    cell_ids: &[String],
    out_dir: &Path,
) -> Result<PanelsContext, Stage3Error> {
    let (mappings, warnings, reverse_index) =
        build_mappings(panels, gene_index, expr.expr.n_genes());
    let mut per_cell = Vec::with_capacity(cell_ids.len());

    let report_path = out_dir.join("panels_report.tsv");
    let mut writer = std::io::BufWriter::new(std::fs::File::create(&report_path)?);

    write_warnings(&mut writer, &warnings)?;
    writer.write_all(b"cell_id\tpanel_id\taxis\tsum\thits\tcoverage\trequired_missing\n")?;

    for (cell_idx, barcode) in cell_ids.iter().enumerate() {
        let mut accums = vec![PanelAccum { sum: 0.0, hits: 0 }; panels.panels.len()];
        let mut last_row_hit = vec![u32::MAX; panels.panels.len()];
        let cell_stats: &CellStats = &expr.cell_stats[cell_idx];
        let inv_denom = if expr.normalization.enabled {
            expr.normalization.scale / (cell_stats.libsize as f32 + expr.normalization.epsilon)
        } else {
            1.0
        };

        expr.expr.for_each_cell_raw(cell_idx, |row, raw_value| {
            let row_usize = row as usize;
            if row_usize >= reverse_index.len() || reverse_index[row_usize].is_empty() {
                return;
            }
            let value = if expr.normalization.enabled {
                (raw_value as f32 * inv_denom).ln_1p()
            } else {
                raw_value as f32
            };
            for (panel_idx, weight) in &reverse_index[row_usize] {
                let acc = &mut accums[*panel_idx];
                acc.sum += value * *weight;
                if last_row_hit[*panel_idx] != row {
                    acc.hits += 1;
                    last_row_hit[*panel_idx] = row;
                }
            }
        });

        let mut required_missing = vec![0u32; panels.panels.len()];
        for (panel_idx, panel) in panels.panels.iter().enumerate() {
            let required_total = mappings[panel_idx].required_total as u32;
            let hits = accums[panel_idx].hits;
            let coverage = if required_total == 0 {
                0.0
            } else {
                (hits as f32 / required_total as f32).min(1.0).max(0.0)
            };
            let missing = required_total.saturating_sub(hits.min(required_total));
            required_missing[panel_idx] = missing;

            let sum = accums[panel_idx].sum;
            let line = format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                barcode,
                panel.id,
                panel.axis,
                format_f32(sum),
                hits,
                format_f32(coverage),
                missing
            );
            writer.write_all(line.as_bytes())?;
        }

        per_cell.push(PanelCellPacked {
            sums: accums.iter().map(|a| a.sum).collect(),
            hits: accums.iter().map(|a| a.hits).collect(),
            required_missing,
        });
    }

    writer.flush()?;

    Ok(PanelsContext {
        panels: panels.clone(),
        mappings,
        warnings,
        cell_ids: cell_ids.to_vec(),
        per_cell,
    })
}

fn build_mappings(
    panels: &PanelSet,
    gene_index: &GeneIndex,
    n_genes: usize,
) -> (
    Vec<GeneMapping>,
    Vec<MappingWarning>,
    Vec<Vec<(usize, f32)>>,
) {
    let mut mappings = Vec::with_capacity(panels.panels.len());
    let mut warnings = Vec::new();
    let mut reverse_index: Vec<Vec<(usize, f32)>> = vec![Vec::new(); n_genes];

    for (panel_idx, panel) in panels.panels.iter().enumerate() {
        let (mapping, warning) = map_panel(panel, gene_index);
        if let Some(w) = warning {
            warnings.push(w);
        }

        let weights = panel.weights.as_ref();
        for (gene_pos, mapped) in mapping.mapped.iter().enumerate() {
            if let Some(row) = mapped {
                let weight = weights
                    .and_then(|w| w.get(gene_pos).copied())
                    .unwrap_or(1.0);
                let row_usize = *row as usize;
                if row_usize < reverse_index.len() {
                    reverse_index[row_usize].push((panel_idx, weight));
                }
            }
        }

        mappings.push(mapping);
    }

    (mappings, warnings, reverse_index)
}

fn format_f32(value: f32) -> String {
    format!("{:.6}", value)
}

fn write_warnings(
    writer: &mut dyn std::io::Write,
    warnings: &[MappingWarning],
) -> Result<(), std::io::Error> {
    if warnings.is_empty() {
        return Ok(());
    }
    writer.write_all(b"# warnings: missing required genes\n")?;
    for warn in warnings {
        let mut line = String::new();
        line.push_str("# ");
        line.push_str(&warn.panel_id);
        line.push(':');
        for (i, gene) in warn.missing_required.iter().enumerate() {
            if i > 0 {
                line.push(',');
            }
            line.push_str(gene);
        }
        line.push('\n');
        writer.write_all(line.as_bytes())?;
    }
    Ok(())
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage3_panels.rs"]
mod tests;
