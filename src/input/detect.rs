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
    let prefix = detect_prefix(dir)?;

    let barcodes = pick_file_with_prefix(dir, &prefix, "barcodes.tsv")
        .ok_or_else(|| InputError::MissingFile("barcodes.tsv[.gz]".to_string()))?;
    let matrix = pick_file_with_prefix(dir, &prefix, "matrix.mtx")
        .ok_or_else(|| InputError::MissingFile("matrix.mtx[.gz]".to_string()))?;

    let features = pick_file_with_prefix(dir, &prefix, "features.tsv");
    let genes = pick_file_with_prefix(dir, &prefix, "genes.tsv");

    match (features, genes) {
        (Some(features_path), None) => Ok(TenXLayout {
            format: TenXFormat::TenXv3,
            matrix_path: matrix,
            features_path,
            barcodes_path: barcodes,
            prefix,
        }),
        (None, Some(genes_path)) => Ok(TenXLayout {
            format: TenXFormat::TenXv2,
            matrix_path: matrix,
            features_path: genes_path,
            barcodes_path: barcodes,
            prefix,
        }),
        (Some(_), Some(_)) => Err(InputError::InvalidTsvRow {
            line: 0,
            reason: "both features.tsv and genes.tsv present".to_string(),
        }),
        (None, None) => Err(InputError::MissingFile(
            "features.tsv/genes.tsv[.gz]".to_string(),
        )),
    }
}

pub fn detect_prefix(dir: &Path) -> Result<Option<String>, InputError> {
    let mut prefixes = Vec::new();
    let entries = std::fs::read_dir(dir)?;
    for entry in entries {
        let entry = entry?;
        let file_name = entry.file_name();
        let file_name = file_name.to_string_lossy();
        if let Some(prefix) = prefix_from_name(&file_name, "matrix.mtx") {
            prefixes.push(prefix);
        } else if let Some(prefix) = prefix_from_name(&file_name, "features.tsv") {
            prefixes.push(prefix);
        } else if let Some(prefix) = prefix_from_name(&file_name, "genes.tsv") {
            prefixes.push(prefix);
        } else if let Some(prefix) = prefix_from_name(&file_name, "barcodes.tsv") {
            prefixes.push(prefix);
        }
    }

    prefixes.sort();
    prefixes.dedup();
    if prefixes.is_empty() {
        Ok(None)
    } else if prefixes.len() == 1 {
        Ok(Some(prefixes.remove(0)))
    } else {
        Err(InputError::InvalidTsvRow {
            line: 0,
            reason: "multiple dataset prefixes detected".to_string(),
        })
    }
}

pub fn resolve_shared_cache_file_name(prefix: Option<&str>) -> String {
    match prefix {
        Some(p) if !p.is_empty() => format!("{}.kira-organelle.bin", p),
        _ => "kira-organelle.bin".to_string(),
    }
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

fn prefix_from_name(name: &str, base: &str) -> Option<String> {
    let plain = format!("_{}", base);
    let gz = format!("_{}.gz", base);
    if name.ends_with(&plain) {
        let prefix = name.trim_end_matches(&plain).to_string();
        if !prefix.is_empty() {
            return Some(prefix);
        }
    }
    if name.ends_with(&gz) {
        let prefix = name.trim_end_matches(&gz).to_string();
        if !prefix.is_empty() {
            return Some(prefix);
        }
    }
    None
}

fn pick_file_with_prefix(dir: &Path, prefix: &Option<String>, base: &str) -> Option<PathBuf> {
    if let Some(prefix) = prefix {
        return pick_file(dir, &format!("{}_{}", prefix, base));
    }
    pick_file(dir, base)
}

fn pick_file(dir: &Path, base: &str) -> Option<PathBuf> {
    let plain = dir.join(base);
    let gz = dir.join(format!("{}.gz", base));

    let plain_exists = plain.is_file();
    let gz_exists = gz.is_file();

    match (plain_exists, gz_exists) {
        (true, true) => {
            warn!(file = base, "both plain and .gz present; choosing plain");
            Some(plain)
        }
        (true, false) => Some(plain),
        (false, true) => Some(gz),
        (false, false) => None,
    }
}

#[cfg(test)]
#[path = "../../tests/src_inline/input/detect.rs"]
mod tests;
