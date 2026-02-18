use std::collections::{HashMap, HashSet};
use std::path::Path;

use crate::input::{InputError, open_reader};

#[derive(Debug, Default, Clone)]
pub struct MetaStats {
    pub matched: usize,
    pub missing: usize,
    pub duplicate_rows: usize,
    pub sample_counts: Option<HashMap<String, usize>>,
}

pub fn read_meta(path: &Path, barcodes: &[String]) -> Result<MetaStats, InputError> {
    let mut reader = open_reader(path)?;
    let mut line = String::new();

    let mut header_line = String::new();
    let read = reader.read_line(&mut header_line)?;
    if read == 0 {
        return Err(InputError::InvalidTsvRow {
            line: 0,
            reason: "empty meta file".to_string(),
        });
    }

    let header = header_line.trim_end_matches(['\n', '\r']);
    let columns: Vec<&str> = header.split('\t').collect();
    let cell_idx = columns
        .iter()
        .position(|c| *c == "cell_id")
        .ok_or_else(|| InputError::MissingMetaColumn("cell_id".to_string()))?;
    let sample_idx = columns.iter().position(|c| *c == "sample_id");

    let barcode_set: HashSet<&str> = barcodes.iter().map(|s| s.as_str()).collect();
    let mut seen_cells: HashSet<String> = HashSet::new();

    let mut stats = MetaStats::default();
    if sample_idx.is_some() {
        stats.sample_counts = Some(HashMap::new());
    }

    let mut line_no = 1usize;
    loop {
        line.clear();
        let read = reader.read_line(&mut line)?;
        if read == 0 {
            break;
        }
        line_no += 1;
        let value = line.trim_end_matches(['\n', '\r']);
        if value.is_empty() {
            continue;
        }
        let parts: Vec<&str> = value.split('\t').collect();
        if cell_idx >= parts.len() {
            return Err(InputError::MissingMetaCellId(line_no));
        }
        let cell_id = parts[cell_idx];
        if cell_id.is_empty() {
            return Err(InputError::MissingMetaCellId(line_no));
        }
        if !seen_cells.insert(cell_id.to_string()) {
            stats.duplicate_rows += 1;
            continue;
        }
        if barcode_set.contains(cell_id) {
            stats.matched += 1;
        } else {
            stats.missing += 1;
        }
        if let (Some(sample_col), Some(counts)) = (sample_idx, stats.sample_counts.as_mut()) {
            if sample_col < parts.len() {
                let sample_id = parts[sample_col];
                if !sample_id.is_empty() {
                    *counts.entry(sample_id.to_string()).or_insert(0) += 1;
                }
            }
        }
    }

    Ok(stats)
}

pub fn read_meta_mapping(
    path: &Path,
    barcodes: &[String],
) -> Result<(Vec<String>, MetaStats), InputError> {
    let mut reader = open_reader(path)?;
    let mut line = String::new();

    let mut header_line = String::new();
    let read = reader.read_line(&mut header_line)?;
    if read == 0 {
        return Err(InputError::InvalidTsvRow {
            line: 0,
            reason: "empty meta file".to_string(),
        });
    }

    let header = header_line.trim_end_matches(['\n', '\r']);
    let columns: Vec<&str> = header.split('\t').collect();
    let cell_idx = columns
        .iter()
        .position(|c| *c == "cell_id")
        .ok_or_else(|| InputError::MissingMetaColumn("cell_id".to_string()))?;
    let sample_idx = columns.iter().position(|c| *c == "sample_id");

    let mut index_by_cell: HashMap<&str, usize> = HashMap::new();
    for (i, c) in barcodes.iter().enumerate() {
        index_by_cell.insert(c.as_str(), i);
    }
    let mut sample_ids = vec![".".to_string(); barcodes.len()];
    let mut seen_cells: HashSet<String> = HashSet::new();

    let mut stats = MetaStats::default();
    if sample_idx.is_some() {
        stats.sample_counts = Some(HashMap::new());
    }

    let mut line_no = 1usize;
    loop {
        line.clear();
        let read = reader.read_line(&mut line)?;
        if read == 0 {
            break;
        }
        line_no += 1;
        let value = line.trim_end_matches(['\n', '\r']);
        if value.is_empty() {
            continue;
        }
        let parts: Vec<&str> = value.split('\t').collect();
        if cell_idx >= parts.len() {
            return Err(InputError::MissingMetaCellId(line_no));
        }
        let cell_id = parts[cell_idx];
        if cell_id.is_empty() {
            return Err(InputError::MissingMetaCellId(line_no));
        }
        if !seen_cells.insert(cell_id.to_string()) {
            stats.duplicate_rows += 1;
            continue;
        }
        if let Some(&idx) = index_by_cell.get(cell_id) {
            stats.matched += 1;
            if let (Some(sample_col), Some(counts)) = (sample_idx, stats.sample_counts.as_mut()) {
                if sample_col < parts.len() {
                    let sample_id = parts[sample_col];
                    if !sample_id.is_empty() {
                        sample_ids[idx] = sample_id.to_string();
                        *counts.entry(sample_id.to_string()).or_insert(0) += 1;
                    }
                }
            }
        } else {
            stats.missing += 1;
        }
    }

    Ok((sample_ids, stats))
}
