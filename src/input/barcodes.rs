use std::path::Path;

use kira_scio::api::{Reader, ReaderOptions};
use kira_scio::detect::DetectedFormat;

use crate::input::InputError;

pub fn read_barcodes(path: &Path) -> Result<Vec<String>, InputError> {
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

    if md.barcodes.is_empty() {
        return Err(InputError::InvalidTsvRow {
            line: 0,
            reason: "no barcodes found".to_string(),
        });
    }

    Ok(md.barcodes)
}
