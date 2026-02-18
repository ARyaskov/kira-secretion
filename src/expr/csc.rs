use std::path::Path;

use crate::expr::normalize::Normalization;
use crate::input::InputError;
use crate::input::mtx::{MatrixHeader, read_entries};

#[derive(Debug, Clone)]
pub struct ExprCsc {
    pub n_genes: usize,
    pub n_cells: usize,
    pub nnz: usize,
    pub col_ptr: Vec<u64>,
    pub row_idx: Vec<u32>,
    pub values: Vec<u32>,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct CellStats {
    pub libsize: u64,
    pub detected: u32,
}

impl ExprCsc {
    pub fn from_mtx(
        path: &Path,
        n_genes: usize,
        n_cells: usize,
        fast: bool,
    ) -> Result<(Self, Vec<CellStats>), InputError> {
        let (header, mut entries) = read_entries(path)?;
        validate_header(&header, n_genes, n_cells, fast)?;
        if !fast && header.nnz != entries.len() {
            return Err(InputError::InvalidMtxDimensions(
                "nnz count does not match header".to_string(),
            ));
        }

        entries.sort_by(|a, b| match a.0.cmp(&b.0) {
            std::cmp::Ordering::Equal => a.1.cmp(&b.1),
            other => other,
        });

        let mut col_counts = vec![0u64; n_cells];
        for (col, _row, _val) in &entries {
            let col_usize = *col as usize;
            if col_usize >= n_cells {
                return Err(InputError::InvalidMtxDimensions(
                    "column index out of bounds".to_string(),
                ));
            }
            col_counts[col_usize] += 1;
        }

        let mut col_ptr = vec![0u64; n_cells + 1];
        for i in 0..n_cells {
            col_ptr[i + 1] = col_ptr[i] + col_counts[i];
        }

        let nnz = entries.len();
        let mut row_idx = Vec::with_capacity(nnz);
        let mut values = Vec::with_capacity(nnz);
        let mut stats = vec![CellStats::default(); n_cells];

        let mut current_col: Option<u32> = None;
        let mut last_row: u32 = 0;
        for (col, row, val) in entries {
            if !fast {
                if row as usize >= n_genes {
                    return Err(InputError::InvalidMtxDimensions(
                        "row index out of bounds".to_string(),
                    ));
                }
                if col as usize >= n_cells {
                    return Err(InputError::InvalidMtxDimensions(
                        "column index out of bounds".to_string(),
                    ));
                }
            }

            if current_col != Some(col) {
                current_col = Some(col);
                last_row = row;
                stats[col as usize].detected += 1;
            } else if row != last_row {
                stats[col as usize].detected += 1;
                last_row = row;
            }

            stats[col as usize].libsize += val as u64;
            row_idx.push(row);
            values.push(val);
        }

        Ok((
            ExprCsc {
                n_genes,
                n_cells,
                nnz,
                col_ptr,
                row_idx,
                values,
            },
            stats,
        ))
    }

    pub fn iter_cell_norm<'a>(
        &'a self,
        cell_idx: usize,
        norm: &'a Normalization,
        cell_stats: &'a CellStats,
    ) -> impl Iterator<Item = (u32, f32)> + 'a {
        let start = self.col_ptr[cell_idx] as usize;
        let end = self.col_ptr[cell_idx + 1] as usize;
        let rows = &self.row_idx[start..end];
        let vals = &self.values[start..end];

        rows.iter()
            .copied()
            .zip(vals.iter().copied())
            .map(move |(row, v)| {
                let raw = v as f32;
                if norm.enabled {
                    let denom = cell_stats.libsize as f32 + norm.epsilon;
                    let scaled = raw * (norm.scale / denom);
                    (row, scaled.ln_1p())
                } else {
                    (row, raw)
                }
            })
    }

    pub fn iter_cell_raw<'a>(&'a self, cell_idx: usize) -> impl Iterator<Item = (u32, u32)> + 'a {
        let start = self.col_ptr[cell_idx] as usize;
        let end = self.col_ptr[cell_idx + 1] as usize;
        let rows = &self.row_idx[start..end];
        let vals = &self.values[start..end];
        rows.iter().copied().zip(vals.iter().copied())
    }
}

fn validate_header(
    header: &MatrixHeader,
    n_genes: usize,
    n_cells: usize,
    fast: bool,
) -> Result<(), InputError> {
    if !fast {
        if header.n_rows != n_genes || header.n_cols != n_cells {
            return Err(InputError::InvalidMtxDimensions(
                "matrix dims do not match stage1".to_string(),
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "../../tests/src_inline/expr/csc.rs"]
mod tests;
