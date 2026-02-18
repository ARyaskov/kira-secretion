use crate::input::features::GeneIndex;
use crate::panels::defs::PanelDef;

#[derive(Debug, Clone)]
pub struct GeneMapping {
    pub panel_id: String,
    pub mapped: Vec<Option<u32>>,
    pub required_hits: usize,
    pub required_total: usize,
}

#[derive(Debug, Clone)]
pub struct MappingWarning {
    pub panel_id: String,
    pub missing_required: Vec<String>,
}

pub fn map_panel(
    panel: &PanelDef,
    gene_index: &GeneIndex,
) -> (GeneMapping, Option<MappingWarning>) {
    let mut mapped = Vec::with_capacity(panel.genes.len());
    for gene in &panel.genes {
        if let Some(row) = gene_index.first_index_by_symbol.get(&gene.symbol) {
            mapped.push(Some((*row as u32) - 1));
        } else {
            mapped.push(None);
        }
    }

    let mut required_hits = 0usize;
    let mut missing_required = Vec::new();
    for req in &panel.required {
        if gene_index.first_index_by_symbol.contains_key(req) {
            required_hits += 1;
        } else {
            missing_required.push(req.clone());
        }
    }

    let warning = if missing_required.is_empty() {
        None
    } else {
        Some(MappingWarning {
            panel_id: panel.id.clone(),
            missing_required,
        })
    };

    (
        GeneMapping {
            panel_id: panel.id.clone(),
            mapped,
            required_hits,
            required_total: panel.required.len(),
        },
        warning,
    )
}

#[cfg(test)]
#[path = "../../tests/src_inline/panels/mapping.rs"]
mod tests;
