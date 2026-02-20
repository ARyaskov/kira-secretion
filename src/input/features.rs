use std::collections::HashMap;
use std::path::Path;

use kira_scio::api::{Reader, ReaderOptions};
use kira_scio::detect::DetectedFormat;

use crate::input::InputError;

#[derive(Debug, Clone)]
pub struct FeatureRow {
    pub id: String,
    pub symbol: String,
}

#[derive(Debug, Clone)]
pub struct DuplicateGene {
    pub symbol: String,
    pub first_row: usize,
    pub dup_row: usize,
}

#[derive(Debug, Clone)]
pub struct GeneIndex {
    pub rows: Vec<FeatureRow>,
    pub duplicates: Vec<DuplicateGene>,
    pub first_index_by_symbol: HashMap<String, usize>,
}

pub fn read_features(path: &Path) -> Result<GeneIndex, InputError> {
    let md = Reader::with_options(
        path,
        ReaderOptions {
            force_format: Some(DetectedFormat::Mtx10x),
            strict: true,
        },
    )
    .read_metadata()
    .map_err(|e| InputError::InvalidTsvRow {
        line: 0,
        reason: e.message,
    })?;

    let rows = md
        .gene_symbols
        .into_iter()
        .enumerate()
        .map(|(idx, symbol)| FeatureRow {
            id: md
                .gene_ids
                .get(idx)
                .cloned()
                .unwrap_or_else(|| format!("GENE_{}", idx + 1)),
            symbol,
        })
        .collect::<Vec<_>>();

    if rows.is_empty() {
        return Err(InputError::InvalidTsvRow {
            line: 0,
            reason: "no feature rows found".to_string(),
        });
    }

    Ok(build_gene_index(rows))
}

pub fn build_gene_index(rows: Vec<FeatureRow>) -> GeneIndex {
    let mut duplicates = Vec::new();
    let mut first_index_by_symbol: HashMap<String, usize> = HashMap::new();
    for (idx, row) in rows.iter().enumerate() {
        let row_no = idx + 1;
        if let Some(first_row) = first_index_by_symbol.get(&row.symbol).copied() {
            duplicates.push(DuplicateGene {
                symbol: row.symbol.clone(),
                first_row,
                dup_row: row_no,
            });
        } else {
            first_index_by_symbol.insert(row.symbol.clone(), row_no);
        }
    }
    GeneIndex {
        rows,
        duplicates,
        first_index_by_symbol,
    }
}

#[cfg(test)]
#[path = "../../tests/src_inline/input/features.rs"]
mod tests;
