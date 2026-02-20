use std::path::Path;

use kira_scio::api::{Reader, ReaderOptions};
use kira_scio::detect::DetectedFormat;

use crate::input::InputError;

#[derive(Debug, Clone, Copy)]
pub struct MatrixHeader {
    pub n_rows: usize,
    pub n_cols: usize,
    pub nnz: usize,
}

pub fn read_header(path: &Path) -> Result<MatrixHeader, InputError> {
    let matrix = Reader::with_options(
        path,
        ReaderOptions {
            force_format: Some(DetectedFormat::Mtx10x),
            strict: true,
        },
    )
    .read_matrix()
    .map_err(|e| InputError::InvalidMtxHeader(e.message))?;

    Ok(MatrixHeader {
        n_rows: matrix.n_genes,
        n_cols: matrix.n_cells,
        nnz: matrix.values.len(),
    })
}

pub fn count_nnz_lines(path: &Path) -> Result<usize, InputError> {
    Ok(read_header(path)?.nnz)
}

pub fn read_entries(path: &Path) -> Result<(MatrixHeader, Vec<(u32, u32, u32)>), InputError> {
    let matrix = Reader::with_options(
        path,
        ReaderOptions {
            force_format: Some(DetectedFormat::Mtx10x),
            strict: true,
        },
    )
    .read_matrix()
    .map_err(|e| InputError::InvalidMtxHeader(e.message))?;

    let header = MatrixHeader {
        n_rows: matrix.n_genes,
        n_cols: matrix.n_cells,
        nnz: matrix.values.len(),
    };

    let mut entries = Vec::with_capacity(header.nnz);
    for (col, w) in matrix.col_ptr.windows(2).enumerate() {
        for idx in w[0]..w[1] {
            let row = matrix.row_idx[idx];
            let value = matrix.values[idx];
            if value < 0.0 || value.fract().abs() > 1e-6 {
                return Err(InputError::InvalidMtxDimensions(
                    "non-integer matrix value".to_string(),
                ));
            }
            entries.push((col as u32, row as u32, value as u32));
        }
    }

    Ok((header, entries))
}

#[cfg(test)]
#[path = "../../tests/src_inline/input/mtx.rs"]
mod tests;
