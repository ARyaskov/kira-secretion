use std::collections::HashMap;
use std::path::Path;

use crate::input::{InputError, open_reader};

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
    let mut reader = open_reader(path)?;
    let mut line = String::new();
    let mut rows = Vec::new();
    let mut duplicates = Vec::new();
    let mut first_index_by_symbol: HashMap<String, usize> = HashMap::new();

    let mut line_no = 0usize;
    loop {
        line.clear();
        let read = reader.read_line(&mut line)?;
        if read == 0 {
            break;
        }
        line_no += 1;
        let value = line.trim_end_matches(['\n', '\r']);
        if value.is_empty() {
            return Err(InputError::InvalidTsvRow {
                line: line_no,
                reason: "empty row".to_string(),
            });
        }
        let mut parts = value.split('\t');
        let id = parts.next().unwrap_or("").to_string();
        let symbol = parts.next().unwrap_or("").to_string();
        if id.is_empty() || symbol.is_empty() {
            return Err(InputError::InvalidTsvRow {
                line: line_no,
                reason: "expected at least 2 columns".to_string(),
            });
        }

        let row_index = line_no;
        if let Some(first_row) = first_index_by_symbol.get(&symbol).copied() {
            duplicates.push(DuplicateGene {
                symbol: symbol.clone(),
                first_row,
                dup_row: row_index,
            });
        } else {
            first_index_by_symbol.insert(symbol.clone(), row_index);
        }

        rows.push(FeatureRow { id, symbol });
    }

    if rows.is_empty() {
        return Err(InputError::InvalidTsvRow {
            line: 0,
            reason: "no feature rows found".to_string(),
        });
    }

    Ok(GeneIndex {
        rows,
        duplicates,
        first_index_by_symbol,
    })
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
