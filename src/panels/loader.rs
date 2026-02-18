use std::fs;
use std::path::{Path, PathBuf};

use thiserror::Error;

use crate::panels::defs::PanelSet;

#[derive(Debug, Error)]
pub enum PanelLoadError {
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
    #[error("toml parse error: {0}")]
    Toml(#[from] toml::de::Error),
    #[error("no panels found in {0}")]
    Empty(String),
}

#[derive(serde::Deserialize)]
struct PanelFile {
    #[serde(default)]
    panel: Vec<crate::panels::defs::PanelDef>,
}

pub fn load_panels_from_dir(dir: &Path) -> Result<PanelSet, PanelLoadError> {
    let mut files = list_toml_files(dir)?;
    files.sort();

    let mut panels = Vec::new();
    for file in files {
        let text = fs::read_to_string(&file)?;
        let parsed: PanelFile = toml::from_str(&text)?;
        panels.extend(parsed.panel.into_iter());
    }

    if panels.is_empty() {
        return Err(PanelLoadError::Empty(dir.to_string_lossy().to_string()));
    }

    Ok(PanelSet { panels })
}

pub fn default_panels_dir() -> PathBuf {
    let relative = Path::new("assets").join("panels");
    if relative.is_dir() {
        return relative;
    }

    let manifest = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("assets")
        .join("panels");
    if manifest.is_dir() {
        return manifest;
    }

    if let Ok(exe) = std::env::current_exe()
        && let Some(dir) = exe.parent()
    {
        let sibling = dir.join("assets").join("panels");
        if sibling.is_dir() {
            return sibling;
        }
        let parent = dir.join("..").join("assets").join("panels");
        if parent.is_dir() {
            return parent;
        }
    }

    relative
}

fn list_toml_files(dir: &Path) -> Result<Vec<PathBuf>, std::io::Error> {
    let mut files = Vec::new();
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.extension().and_then(|s| s.to_str()) == Some("toml") {
            files.push(path);
        }
    }
    Ok(files)
}

#[cfg(test)]
#[path = "../../tests/src_inline/panels/loader.rs"]
mod tests;
