use std::path::Path;

use thiserror::Error;

use crate::expr::csc::{CellStats, ExprCsc};
use crate::expr::normalize::Normalization;
use crate::input::InputError;
use crate::input::cache::{SharedCacheMapped, mmap_shared_cache, mmap_shared_cache_unchecked};
use crate::pipeline::stage1_load::DatasetCtx;

#[derive(Debug, Error)]
pub enum Stage2Error {
    #[error("input error: {0}")]
    Input(#[from] InputError),
    #[error("cache error: {0}")]
    Cache(#[from] crate::input::cache::CacheError),
}

#[derive(Debug, Clone)]
pub enum ExprMatrix {
    Owned(ExprCsc),
    Shared(SharedCacheMapped),
}

impl ExprMatrix {
    pub fn n_genes(&self) -> usize {
        match self {
            ExprMatrix::Owned(e) => e.n_genes,
            ExprMatrix::Shared(e) => e.n_genes,
        }
    }

    pub fn n_cells(&self) -> usize {
        match self {
            ExprMatrix::Owned(e) => e.n_cells,
            ExprMatrix::Shared(e) => e.n_cells,
        }
    }

    pub fn nnz(&self) -> usize {
        match self {
            ExprMatrix::Owned(e) => e.nnz,
            ExprMatrix::Shared(e) => e.nnz,
        }
    }

    pub fn for_each_cell_norm<F>(
        &self,
        cell_idx: usize,
        norm: &Normalization,
        cell_stats: &CellStats,
        mut f: F,
    ) where
        F: FnMut(u32, f32),
    {
        match self {
            ExprMatrix::Owned(e) => {
                for (row, value) in e.iter_cell_norm(cell_idx, norm, cell_stats) {
                    f(row, value);
                }
            }
            ExprMatrix::Shared(e) => e.for_each_cell_norm(cell_idx, norm, cell_stats, f),
        }
    }

    pub fn for_each_cell_raw<F>(&self, cell_idx: usize, mut f: F)
    where
        F: FnMut(u32, u32),
    {
        match self {
            ExprMatrix::Owned(e) => {
                for (row, value) in e.iter_cell_raw(cell_idx) {
                    f(row, value);
                }
            }
            ExprMatrix::Shared(e) => e.for_each_cell_raw(cell_idx, f),
        }
    }
}

#[derive(Debug, Clone)]
pub struct ExprContext {
    pub expr: ExprMatrix,
    pub cell_stats: Vec<CellStats>,
    pub normalization: Normalization,
}

pub fn run_stage2(
    ctx: &DatasetCtx,
    _out_dir: &Path,
    normalization: Normalization,
    fast: bool,
) -> Result<ExprContext, Stage2Error> {
    if let Some(shared_cache_path) = &ctx.shared_cache_path {
        // Stage 1 already performed strict validation in pipeline mode.
        let shared = mmap_shared_cache_unchecked(shared_cache_path)
            .or_else(|_| mmap_shared_cache(shared_cache_path))?;
        let cell_stats = shared.compute_cell_stats();
        return Ok(ExprContext {
            expr: ExprMatrix::Shared(shared),
            cell_stats,
            normalization,
        });
    }

    let (expr, cell_stats) = ExprCsc::from_mtx(&ctx.matrix_path, ctx.n_genes, ctx.n_cells, fast)?;

    Ok(ExprContext {
        expr: ExprMatrix::Owned(expr),
        cell_stats,
        normalization,
    })
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage2_normalize.rs"]
mod tests;
