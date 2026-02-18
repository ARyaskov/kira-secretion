use std::path::{Path, PathBuf};

use thiserror::Error;
use tracing::warn;

use crate::input::InputError;
use crate::input::barcodes::read_barcodes;
use crate::input::cache::read_shared_cache_metadata;
use crate::input::detect::{
    TenXFormat, TenXLayout, detect_10x_dir, detect_prefix, find_shared_cache_file,
    resolve_shared_cache_file_name,
};
use crate::input::features::{DuplicateGene, FeatureRow, build_gene_index, read_features};
use crate::input::meta::read_meta;
use crate::input::mtx::{count_nnz_lines, read_header};

#[derive(Debug, Error)]
pub enum Stage1Error {
    #[error("input error: {0}")]
    Input(#[from] InputError),
    #[error("cache error: {0}")]
    Cache(#[from] crate::input::cache::CacheError),
    #[error("matrix dimensions do not match features/barcodes")]
    DimensionMismatch {
        expected_rows: usize,
        expected_cols: usize,
        found_rows: usize,
        found_cols: usize,
    },
    #[error("nnz line count mismatch: expected {expected}, found {found}")]
    NnzMismatch { expected: usize, found: usize },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunMode {
    Standalone,
    Pipeline,
}

#[derive(Debug, Clone)]
pub struct DatasetCtx {
    pub format: TenXFormat,
    pub matrix_path: PathBuf,
    pub features_path: PathBuf,
    pub barcodes_path: PathBuf,
    pub shared_cache_path: Option<PathBuf>,
    pub resolved_shared_cache_path: Option<PathBuf>,
    pub gene_index: crate::input::features::GeneIndex,
    pub barcodes: Vec<String>,
    pub n_genes: usize,
    pub n_cells: usize,
    pub nnz: usize,
    pub duplicate_gene_symbols_count: usize,
    pub duplicate_gene_symbols: Vec<DuplicateGene>,
    pub meta_present: bool,
    pub meta_cells_matched: usize,
    pub meta_cells_missing: usize,
}

pub fn run_stage1(
    input_dir: &Path,
    meta_path: Option<&Path>,
    out_dir: &Path,
    fast: bool,
    run_mode: RunMode,
    cache_override: Option<&Path>,
) -> Result<DatasetCtx, Stage1Error> {
    let _ = out_dir;

    if run_mode == RunMode::Pipeline {
        if let Some(cache_path) = cache_override {
            return run_stage1_shared_cache(input_dir, cache_path.to_path_buf(), meta_path);
        }
        let prefix = detect_prefix(input_dir)?;
        let cache_name = resolve_shared_cache_file_name(prefix.as_deref());
        let expected_cache = input_dir.join(cache_name);
        if let Some(cache_path) = find_shared_cache_file(input_dir, prefix.as_deref())? {
            return run_stage1_shared_cache(input_dir, cache_path, meta_path);
        }
        warn!(
            expected_cache = %expected_cache.to_string_lossy(),
            "shared cache not found, falling back to MTX input"
        );
        let layout = detect_10x_dir(input_dir)?;
        let mut ctx = run_stage1_layout(input_dir, layout, meta_path, fast)?;
        ctx.resolved_shared_cache_path = Some(expected_cache);
        return Ok(ctx);
    }

    let layout = detect_10x_dir(input_dir)?;
    run_stage1_layout(input_dir, layout, meta_path, fast)
}

fn run_stage1_shared_cache(
    input_dir: &Path,
    shared_cache_path: PathBuf,
    meta_path: Option<&Path>,
) -> Result<DatasetCtx, Stage1Error> {
    let metadata = read_shared_cache_metadata(&shared_cache_path)?;

    let rows: Vec<FeatureRow> = metadata
        .genes
        .iter()
        .map(|g| FeatureRow {
            id: g.clone(),
            symbol: g.clone(),
        })
        .collect();
    let gene_index = build_gene_index(rows);
    let duplicate_gene_symbols_count = gene_index.duplicates.len();
    let duplicate_gene_symbols = gene_index.duplicates.clone();

    let mut meta_present = false;
    let mut meta_cells_matched = 0usize;
    let mut meta_cells_missing = 0usize;
    if let Some(meta) = meta_path {
        meta_present = true;
        let stats = read_meta(meta, &metadata.barcodes)?;
        meta_cells_matched = stats.matched;
        meta_cells_missing = stats.missing;
    }

    Ok(DatasetCtx {
        format: TenXFormat::Unknown,
        matrix_path: input_dir.join("matrix.mtx"),
        features_path: input_dir.join("features.tsv"),
        barcodes_path: input_dir.join("barcodes.tsv"),
        shared_cache_path: Some(shared_cache_path.clone()),
        resolved_shared_cache_path: Some(shared_cache_path),
        gene_index,
        barcodes: metadata.barcodes,
        n_genes: metadata.n_genes,
        n_cells: metadata.n_cells,
        nnz: metadata.nnz,
        duplicate_gene_symbols_count,
        duplicate_gene_symbols,
        meta_present,
        meta_cells_matched,
        meta_cells_missing,
    })
}

fn run_stage1_layout(
    input_dir: &Path,
    layout: TenXLayout,
    meta_path: Option<&Path>,
    fast: bool,
) -> Result<DatasetCtx, Stage1Error> {
    let barcodes = read_barcodes(&layout.barcodes_path)?;
    let gene_index = read_features(&layout.features_path)?;
    let n_genes = gene_index.rows.len();
    let duplicate_gene_symbols_count = gene_index.duplicates.len();
    let duplicate_gene_symbols = gene_index.duplicates.clone();
    let header = read_header(&layout.matrix_path)?;

    if header.n_rows != gene_index.rows.len() || header.n_cols != barcodes.len() {
        return Err(Stage1Error::DimensionMismatch {
            expected_rows: gene_index.rows.len(),
            expected_cols: barcodes.len(),
            found_rows: header.n_rows,
            found_cols: header.n_cols,
        });
    }

    if !fast {
        let counted = count_nnz_lines(&layout.matrix_path)?;
        if counted != header.nnz {
            return Err(Stage1Error::NnzMismatch {
                expected: header.nnz,
                found: counted,
            });
        }
    }

    let mut meta_present = false;
    let mut meta_cells_matched = 0usize;
    let mut meta_cells_missing = 0usize;

    if let Some(meta) = meta_path {
        meta_present = true;
        let stats = read_meta(meta, &barcodes)?;
        meta_cells_matched = stats.matched;
        meta_cells_missing = stats.missing;
    }

    Ok(DatasetCtx {
        format: layout.format,
        matrix_path: layout.matrix_path,
        features_path: layout.features_path,
        barcodes_path: layout.barcodes_path,
        shared_cache_path: None,
        resolved_shared_cache_path: layout
            .prefix
            .as_deref()
            .map(|p| input_dir.join(resolve_shared_cache_file_name(Some(p)))),
        gene_index,
        barcodes,
        n_genes,
        n_cells: header.n_cols,
        nnz: header.nnz,
        duplicate_gene_symbols_count,
        duplicate_gene_symbols,
        meta_present,
        meta_cells_matched,
        meta_cells_missing,
    })
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage1_load.rs"]
mod tests;
