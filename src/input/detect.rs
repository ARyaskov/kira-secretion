use std::path::{Path, PathBuf};

use tracing::warn;

use crate::input::InputError;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TenXFormat {
    TenXv2,
    TenXv3,
    Unknown,
}

impl std::fmt::Display for TenXFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TenXFormat::TenXv2 => write!(f, "tenx_v2"),
            TenXFormat::TenXv3 => write!(f, "tenx_v3"),
            TenXFormat::Unknown => write!(f, "unknown"),
        }
    }
}

#[derive(Debug, Clone)]
pub struct TenXLayout {
    pub format: TenXFormat,
    pub matrix_path: PathBuf,
    pub features_path: PathBuf,
    pub barcodes_path: PathBuf,
    pub prefix: Option<String>,
}

pub fn detect_10x_dir(dir: &Path) -> Result<TenXLayout, InputError> {
    let ds = kira_scio::discover(dir).map_err(|e| InputError::MissingFile(e.message))?;
    let barcodes = ds
        .barcodes
        .ok_or_else(|| InputError::MissingFile("barcodes.tsv[.gz]".to_string()))?;

    let (features_path, format) = if let Some(features) = ds.features {
        (features, TenXFormat::TenXv3)
    } else if let Some(genes) = ds.genes {
        (genes, TenXFormat::TenXv2)
    } else {
        return Err(InputError::MissingFile(
            "features.tsv/genes.tsv[.gz]".to_string(),
        ));
    };

    Ok(TenXLayout {
        format,
        matrix_path: ds.matrix,
        features_path,
        barcodes_path: barcodes,
        prefix: ds.prefix,
    })
}

pub fn detect_prefix(dir: &Path) -> Result<Option<String>, InputError> {
    kira_scio::detect_prefix(dir).map_err(|e| InputError::InvalidTsvRow {
        line: 0,
        reason: e.to_string(),
    })
}

pub fn resolve_shared_cache_file_name(prefix: Option<&str>) -> String {
    kira_scio::resolve_shared_cache_filename(prefix)
}

pub fn find_shared_cache_file(
    dir: &Path,
    prefix: Option<&str>,
) -> Result<Option<PathBuf>, InputError> {
    let expected = dir.join(resolve_shared_cache_file_name(prefix));
    if expected.is_file() {
        return Ok(Some(expected));
    }

    let mut candidates = Vec::new();
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        if !entry.file_type()?.is_file() {
            continue;
        }
        let name = entry.file_name();
        let name = name.to_string_lossy();
        if name.ends_with("kira-organelle.bin") {
            candidates.push(entry.path());
        }
    }

    candidates.sort();
    if candidates.len() > 1 {
        warn!(
            count = candidates.len(),
            selected = %candidates[0].to_string_lossy(),
            "multiple shared cache files found; selecting lexicographically first"
        );
    }

    Ok(candidates.into_iter().next())
}

#[cfg(test)]
#[path = "../../tests/src_inline/input/detect.rs"]
mod tests;
