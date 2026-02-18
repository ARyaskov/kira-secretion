use std::path::Path;

use crate::input::{InputError, open_reader};

pub fn read_barcodes(path: &Path) -> Result<Vec<String>, InputError> {
    let mut reader = open_reader(path)?;
    let mut line = String::new();
    let mut barcodes = Vec::new();
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
            return Err(InputError::EmptyBarcode(line_no));
        }
        barcodes.push(value.to_string());
    }
    if barcodes.is_empty() {
        return Err(InputError::InvalidTsvRow {
            line: 0,
            reason: "no barcodes found".to_string(),
        });
    }
    Ok(barcodes)
}
