pub mod barcodes;
pub mod cache;
pub mod detect;
pub mod features;
pub mod meta;
pub mod mtx;

use std::path::{Path, PathBuf};
use std::{fmt, io};

use thiserror::Error;

#[derive(Error, Debug)]
pub enum InputError {
    #[error("missing required file: {0}")]
    MissingFile(String),
    #[error("invalid matrix market header: {0}")]
    InvalidMtxHeader(String),
    #[error("invalid matrix dimensions: {0}")]
    InvalidMtxDimensions(String),
    #[error("invalid TSV row at line {line}: {reason}")]
    InvalidTsvRow { line: usize, reason: String },
    #[error("empty barcode at line {0}")]
    EmptyBarcode(usize),
    #[error("meta file missing required column: {0}")]
    MissingMetaColumn(String),
    #[error("meta row missing cell_id at line {0}")]
    MissingMetaCellId(usize),
    #[error("unsupported gzip input without feature enabled: {0}")]
    GzipNotEnabled(PathBuf),
    #[error("io error: {0}")]
    Io(#[from] io::Error),
}

pub fn open_reader(path: &Path) -> Result<Box<dyn io::BufRead>, InputError> {
    let file = std::fs::File::open(path)?;
    if path.extension().and_then(|s| s.to_str()) == Some("gz") {
        #[cfg(feature = "gz")]
        {
            let decoder = flate2::read::GzDecoder::new(file);
            return Ok(Box::new(io::BufReader::new(decoder)));
        }
        #[cfg(not(feature = "gz"))]
        {
            return Err(InputError::GzipNotEnabled(path.to_path_buf()));
        }
    }
    Ok(Box::new(io::BufReader::new(file)))
}

pub fn path_display(path: &Path) -> impl fmt::Display + '_ {
    path.to_string_lossy()
}
