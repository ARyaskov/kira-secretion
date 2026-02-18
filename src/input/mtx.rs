use std::path::Path;

use crate::input::{InputError, open_reader};

#[derive(Debug, Clone, Copy)]
pub struct MatrixHeader {
    pub n_rows: usize,
    pub n_cols: usize,
    pub nnz: usize,
}

pub fn read_header(path: &Path) -> Result<MatrixHeader, InputError> {
    let mut reader = open_reader(path)?;
    let mut line = String::new();
    let read = reader.read_line(&mut line)?;
    if read == 0 {
        return Err(InputError::InvalidMtxHeader("empty file".to_string()));
    }
    let header = line.trim_end_matches(['\n', '\r']);
    let parts: Vec<&str> = header.split_whitespace().collect();
    if parts.len() < 5 {
        return Err(InputError::InvalidMtxHeader(
            "expected MatrixMarket banner".to_string(),
        ));
    }
    if !parts[0].eq_ignore_ascii_case("%%MatrixMarket")
        || !parts[1].eq_ignore_ascii_case("matrix")
        || !parts[2].eq_ignore_ascii_case("coordinate")
        || !parts[3].eq_ignore_ascii_case("integer")
    {
        return Err(InputError::InvalidMtxHeader(
            "unsupported MatrixMarket format".to_string(),
        ));
    }

    let dims = read_dims(&mut reader)?;
    Ok(dims)
}

pub fn count_nnz_lines(path: &Path) -> Result<usize, InputError> {
    let mut reader = open_reader(path)?;
    let mut line = String::new();
    let read = reader.read_line(&mut line)?;
    if read == 0 {
        return Err(InputError::InvalidMtxHeader("empty file".to_string()));
    }

    let _ = read_dims(&mut reader)?;

    let mut count = 0usize;
    loop {
        line.clear();
        let read = reader.read_line(&mut line)?;
        if read == 0 {
            break;
        }
        let value = line.trim_end_matches(['\n', '\r']);
        if value.is_empty() || value.starts_with('%') {
            continue;
        }
        count += 1;
    }
    Ok(count)
}

pub fn read_entries(path: &Path) -> Result<(MatrixHeader, Vec<(u32, u32, u32)>), InputError> {
    let mut reader = open_reader(path)?;
    let mut line = String::new();

    let read = reader.read_line(&mut line)?;
    if read == 0 {
        return Err(InputError::InvalidMtxHeader("empty file".to_string()));
    }
    let header = line.trim_end_matches(['\n', '\r']);
    let parts: Vec<&str> = header.split_whitespace().collect();
    if parts.len() < 5 {
        return Err(InputError::InvalidMtxHeader(
            "expected MatrixMarket banner".to_string(),
        ));
    }
    if !parts[0].eq_ignore_ascii_case("%%MatrixMarket")
        || !parts[1].eq_ignore_ascii_case("matrix")
        || !parts[2].eq_ignore_ascii_case("coordinate")
        || !parts[3].eq_ignore_ascii_case("integer")
    {
        return Err(InputError::InvalidMtxHeader(
            "unsupported MatrixMarket format".to_string(),
        ));
    }

    let header = read_dims(&mut reader)?;
    let mut entries = Vec::with_capacity(header.nnz);

    loop {
        line.clear();
        let read = reader.read_line(&mut line)?;
        if read == 0 {
            break;
        }
        let value = line.trim_end_matches(['\n', '\r']);
        if value.is_empty() || value.starts_with('%') {
            continue;
        }
        let parts: Vec<&str> = value.split_whitespace().collect();
        if parts.len() < 3 {
            return Err(InputError::InvalidTsvRow {
                line: 0,
                reason: "invalid mtx entry".to_string(),
            });
        }
        let row: u32 = parts[0]
            .parse::<u32>()
            .map_err(|_| InputError::InvalidMtxDimensions("invalid row".to_string()))?;
        let col: u32 = parts[1]
            .parse::<u32>()
            .map_err(|_| InputError::InvalidMtxDimensions("invalid col".to_string()))?;
        let val: u32 = parts[2]
            .parse::<u32>()
            .map_err(|_| InputError::InvalidMtxDimensions("invalid value".to_string()))?;
        if row == 0 || col == 0 {
            return Err(InputError::InvalidMtxDimensions(
                "matrix indices must be 1-based".to_string(),
            ));
        }
        entries.push((col - 1, row - 1, val));
    }

    Ok((header, entries))
}

fn read_dims(reader: &mut dyn std::io::BufRead) -> Result<MatrixHeader, InputError> {
    let mut line = String::new();
    loop {
        line.clear();
        let read = reader.read_line(&mut line)?;
        if read == 0 {
            return Err(InputError::InvalidMtxDimensions(
                "missing dimensions".to_string(),
            ));
        }
        let value = line.trim_end_matches(['\n', '\r']);
        if value.is_empty() || value.starts_with('%') {
            continue;
        }
        let parts: Vec<&str> = value.split_whitespace().collect();
        if parts.len() < 3 {
            return Err(InputError::InvalidMtxDimensions(
                "expected three integers".to_string(),
            ));
        }
        let n_rows: usize = parts[0]
            .parse()
            .map_err(|_| InputError::InvalidMtxDimensions("invalid rows".to_string()))?;
        let n_cols: usize = parts[1]
            .parse()
            .map_err(|_| InputError::InvalidMtxDimensions("invalid cols".to_string()))?;
        let nnz: usize = parts[2]
            .parse()
            .map_err(|_| InputError::InvalidMtxDimensions("invalid nnz".to_string()))?;
        if n_rows == 0 || n_cols == 0 {
            return Err(InputError::InvalidMtxDimensions(
                "zero dimensions".to_string(),
            ));
        }
        return Ok(MatrixHeader {
            n_rows,
            n_cols,
            nnz,
        });
    }
}

#[cfg(test)]
#[path = "../../tests/src_inline/input/mtx.rs"]
mod tests;
