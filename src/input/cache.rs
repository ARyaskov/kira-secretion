use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::sync::Arc;

use crc::{CRC_64_ECMA_182, Crc};
use memmap2::Mmap;
use thiserror::Error;

use crate::expr::csc::{CellStats, ExprCsc};
use crate::expr::normalize::Normalization;
use crate::simd;

const MAGIC_EXPR: &[u8; 8] = b"KIRAEXPR";
const VERSION_EXPR: u32 = 1;

const SHARED_MAGIC: &[u8; 4] = b"KORG";
const SHARED_ENDIAN_TAG: u32 = 0x1234_5678;
const SHARED_HEADER_SIZE: usize = 256;
const CRC64: Crc<u64> = Crc::<u64>::new(&CRC_64_ECMA_182);

#[derive(Debug, Error)]
pub enum CacheError {
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
    #[error("invalid cache magic")]
    InvalidMagic,
    #[error("unsupported cache version: {0}")]
    UnsupportedVersion(u32),
    #[error("invalid cache format: {0}")]
    InvalidFormat(String),
}

#[derive(Debug, Clone)]
pub struct SharedCacheMetadata {
    pub n_genes: usize,
    pub n_cells: usize,
    pub nnz: usize,
    pub genes: Vec<String>,
    pub barcodes: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct SharedCacheMapped {
    mmap: Arc<Mmap>,
    pub n_genes: usize,
    pub n_cells: usize,
    pub nnz: usize,
    pub genes: Vec<String>,
    pub barcodes: Vec<String>,
    col_ptr_offset: usize,
    row_idx_offset: usize,
    values_offset: usize,
}

impl SharedCacheMapped {
    pub fn metadata(&self) -> SharedCacheMetadata {
        SharedCacheMetadata {
            n_genes: self.n_genes,
            n_cells: self.n_cells,
            nnz: self.nnz,
            genes: self.genes.clone(),
            barcodes: self.barcodes.clone(),
        }
    }

    pub fn col_ptr_at(&self, i: usize) -> u64 {
        let base = self.col_ptr_offset + i * 8;
        read_u64_slice(&self.mmap[base..base + 8])
    }

    pub fn row_idx_at(&self, i: usize) -> u32 {
        let base = self.row_idx_offset + i * 4;
        read_u32_slice(&self.mmap[base..base + 4])
    }

    pub fn value_at(&self, i: usize) -> u32 {
        let base = self.values_offset + i * 4;
        read_u32_slice(&self.mmap[base..base + 4])
    }

    pub fn compute_cell_stats(&self) -> Vec<CellStats> {
        let mut stats = vec![CellStats::default(); self.n_cells];
        for (cell, stat) in stats.iter_mut().enumerate().take(self.n_cells) {
            let start = self.col_ptr_at(cell) as usize;
            let end = self.col_ptr_at(cell + 1) as usize;
            stat.detected = (end - start) as u32;
            stat.libsize = self.sum_values_range(start, end);
        }
        stats
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
        self.for_each_cell_raw(cell_idx, |row, raw_count| {
            let raw = raw_count as f32;
            let out = if norm.enabled {
                let denom = cell_stats.libsize as f32 + norm.epsilon;
                let scaled = raw * (norm.scale / denom);
                scaled.ln_1p()
            } else {
                raw
            };
            f(row, out);
        });
    }

    pub fn for_each_cell_raw<F>(&self, cell_idx: usize, mut f: F)
    where
        F: FnMut(u32, u32),
    {
        let start = self.col_ptr_at(cell_idx) as usize;
        let end = self.col_ptr_at(cell_idx + 1) as usize;

        #[cfg(target_endian = "little")]
        {
            // SAFETY: row/value sections are validated and range comes from validated col_ptr.
            unsafe {
                let rows_ptr =
                    self.mmap.as_ptr().add(self.row_idx_offset + start * 4) as *const u32;
                let vals_ptr = self.mmap.as_ptr().add(self.values_offset + start * 4) as *const u32;
                let len = end - start;
                let rows = std::slice::from_raw_parts(rows_ptr, len);
                let vals = std::slice::from_raw_parts(vals_ptr, len);
                for i in 0..len {
                    f(rows[i], vals[i]);
                }
            }
            return;
        }

        #[cfg(not(target_endian = "little"))]
        {
            for i in start..end {
                f(self.row_idx_at(i), self.value_at(i));
            }
        }
    }

    fn sum_values_range(&self, start: usize, end: usize) -> u64 {
        #[cfg(target_endian = "little")]
        {
            // SAFETY: values_u32 section is validated; bounds are constrained by caller using col_ptr.
            unsafe {
                let ptr = self.mmap.as_ptr().add(self.values_offset + start * 4) as *const u32;
                let slice = std::slice::from_raw_parts(ptr, end - start);
                simd::sum_u32(slice)
            }
        }
        #[cfg(not(target_endian = "little"))]
        {
            let mut sum = 0u64;
            for i in start..end {
                sum += self.value_at(i) as u64;
            }
            sum
        }
    }
}

pub fn read_shared_cache_metadata(path: &Path) -> Result<SharedCacheMetadata, CacheError> {
    let mapped = mmap_shared_cache(path)?;
    Ok(mapped.metadata())
}

pub fn mmap_shared_cache(path: &Path) -> Result<SharedCacheMapped, CacheError> {
    let file = File::open(path)?;
    let mmap = {
        // SAFETY: mapping file read-only and holding Arc<Mmap> for lifetime of view.
        unsafe { Mmap::map(&file)? }
    };
    parse_shared_cache(Arc::new(mmap), true)
}

pub fn mmap_shared_cache_unchecked(path: &Path) -> Result<SharedCacheMapped, CacheError> {
    let file = File::open(path)?;
    let mmap = {
        // SAFETY: mapping file read-only and holding Arc<Mmap> for lifetime of view.
        unsafe { Mmap::map(&file)? }
    };
    parse_shared_cache(Arc::new(mmap), false)
}

fn parse_shared_cache(
    mmap: Arc<Mmap>,
    validate_csc_strict: bool,
) -> Result<SharedCacheMapped, CacheError> {
    if mmap.len() < SHARED_HEADER_SIZE {
        return Err(CacheError::InvalidFormat(
            "file smaller than header".to_string(),
        ));
    }

    let header = &mmap[..SHARED_HEADER_SIZE];
    if &header[0..4] != SHARED_MAGIC {
        return Err(CacheError::InvalidMagic);
    }
    let version_major = read_u16_slice(&header[4..6]);
    let version_minor = read_u16_slice(&header[6..8]);
    if version_major != 1 {
        return Err(CacheError::UnsupportedVersion(version_major as u32));
    }
    if version_minor != 0 {
        return Err(CacheError::InvalidFormat(
            "unsupported minor version".to_string(),
        ));
    }
    let endian_tag = read_u32_slice(&header[8..12]);
    if endian_tag != SHARED_ENDIAN_TAG {
        return Err(CacheError::InvalidFormat("invalid endian tag".to_string()));
    }
    let header_size = read_u32_slice(&header[12..16]) as usize;
    if header_size != SHARED_HEADER_SIZE {
        return Err(CacheError::InvalidFormat("invalid header size".to_string()));
    }

    let n_genes = read_u64_slice(&header[16..24]) as usize;
    let n_cells = read_u64_slice(&header[24..32]) as usize;
    let nnz = read_u64_slice(&header[32..40]) as usize;

    let genes_table_offset = read_u64_slice(&header[40..48]) as usize;
    let genes_table_bytes = read_u64_slice(&header[48..56]) as usize;
    let barcodes_table_offset = read_u64_slice(&header[56..64]) as usize;
    let barcodes_table_bytes = read_u64_slice(&header[64..72]) as usize;

    let col_ptr_offset = read_u64_slice(&header[72..80]) as usize;
    let row_idx_offset = read_u64_slice(&header[80..88]) as usize;
    let values_offset = read_u64_slice(&header[88..96]) as usize;

    let n_blocks = read_u64_slice(&header[96..104]);
    let blocks_offset = read_u64_slice(&header[104..112]);
    let file_bytes = read_u64_slice(&header[112..120]) as usize;
    let header_crc64 = read_u64_slice(&header[120..128]);

    if n_blocks != 0 || blocks_offset != 0 {
        return Err(CacheError::InvalidFormat(
            "unsupported optional blocks in v1".to_string(),
        ));
    }
    if file_bytes != mmap.len() {
        return Err(CacheError::InvalidFormat(
            "file_bytes does not match file length".to_string(),
        ));
    }

    let mut header_for_crc = header.to_vec();
    header_for_crc[120..128].fill(0);
    let crc = CRC64.checksum(&header_for_crc);
    if crc != header_crc64 {
        return Err(CacheError::InvalidFormat(
            "header CRC64 mismatch".to_string(),
        ));
    }

    check_bounds(
        mmap.len(),
        genes_table_offset,
        genes_table_bytes,
        "genes table",
    )?;
    check_bounds(
        mmap.len(),
        barcodes_table_offset,
        barcodes_table_bytes,
        "barcodes table",
    )?;

    let col_ptr_bytes = (n_cells + 1)
        .checked_mul(8)
        .ok_or_else(|| CacheError::InvalidFormat("col_ptr size overflow".to_string()))?;
    let row_idx_bytes = nnz
        .checked_mul(4)
        .ok_or_else(|| CacheError::InvalidFormat("row_idx size overflow".to_string()))?;
    let values_bytes = nnz
        .checked_mul(4)
        .ok_or_else(|| CacheError::InvalidFormat("values size overflow".to_string()))?;

    check_bounds(mmap.len(), col_ptr_offset, col_ptr_bytes, "col_ptr")?;
    check_bounds(mmap.len(), row_idx_offset, row_idx_bytes, "row_idx")?;
    check_bounds(mmap.len(), values_offset, values_bytes, "values")?;

    let genes = parse_string_table(
        &mmap,
        genes_table_offset,
        genes_table_bytes,
        n_genes,
        "genes",
    )?;
    let barcodes = parse_string_table(
        &mmap,
        barcodes_table_offset,
        barcodes_table_bytes,
        n_cells,
        "barcodes",
    )?;

    if validate_csc_strict {
        validate_csc(&mmap, n_genes, n_cells, nnz, col_ptr_offset, row_idx_offset)?;
    }

    Ok(SharedCacheMapped {
        mmap,
        n_genes,
        n_cells,
        nnz,
        genes,
        barcodes,
        col_ptr_offset,
        row_idx_offset,
        values_offset,
    })
}

fn check_bounds(
    file_len: usize,
    offset: usize,
    bytes: usize,
    label: &str,
) -> Result<(), CacheError> {
    if offset < SHARED_HEADER_SIZE {
        return Err(CacheError::InvalidFormat(format!(
            "{} offset before header",
            label
        )));
    }
    let end = offset
        .checked_add(bytes)
        .ok_or_else(|| CacheError::InvalidFormat(format!("{} bounds overflow", label)))?;
    if end > file_len {
        return Err(CacheError::InvalidFormat(format!(
            "{} out of file bounds",
            label
        )));
    }
    Ok(())
}

fn parse_string_table(
    mmap: &[u8],
    offset: usize,
    bytes: usize,
    expected_count: usize,
    label: &str,
) -> Result<Vec<String>, CacheError> {
    if bytes < 4 {
        return Err(CacheError::InvalidFormat(format!(
            "{} table too small",
            label
        )));
    }
    let table = &mmap[offset..offset + bytes];
    let count = read_u32_slice(&table[0..4]) as usize;
    if count != expected_count {
        return Err(CacheError::InvalidFormat(format!(
            "{} table count mismatch",
            label
        )));
    }

    let offsets_bytes = (count + 1)
        .checked_mul(4)
        .ok_or_else(|| CacheError::InvalidFormat(format!("{} offsets overflow", label)))?;
    if bytes < 4 + offsets_bytes {
        return Err(CacheError::InvalidFormat(format!(
            "{} table missing offsets",
            label
        )));
    }

    let mut offsets = Vec::with_capacity(count + 1);
    for i in 0..=count {
        let start = 4 + i * 4;
        offsets.push(read_u32_slice(&table[start..start + 4]) as usize);
    }

    for i in 0..count {
        if offsets[i + 1] < offsets[i] {
            return Err(CacheError::InvalidFormat(format!(
                "{} offsets not monotonic",
                label
            )));
        }
    }

    let blob = &table[4 + offsets_bytes..];
    if offsets[count] != blob.len() {
        return Err(CacheError::InvalidFormat(format!(
            "{} terminal offset mismatch",
            label
        )));
    }

    let mut out = Vec::with_capacity(count);
    for i in 0..count {
        let start = offsets[i];
        let end = offsets[i + 1];
        let s = std::str::from_utf8(&blob[start..end]).map_err(|_| {
            CacheError::InvalidFormat(format!("{} contains invalid UTF-8 string", label))
        })?;
        out.push(s.to_string());
    }
    Ok(out)
}

fn validate_csc(
    mmap: &[u8],
    n_genes: usize,
    n_cells: usize,
    nnz: usize,
    col_ptr_offset: usize,
    row_idx_offset: usize,
) -> Result<(), CacheError> {
    let mut prev_col_ptr = 0u64;
    for i in 0..=n_cells {
        let base = col_ptr_offset + i * 8;
        let v = read_u64_slice(&mmap[base..base + 8]);
        if i == 0 && v != 0 {
            return Err(CacheError::InvalidFormat(
                "col_ptr[0] must be 0".to_string(),
            ));
        }
        if i > 0 && v < prev_col_ptr {
            return Err(CacheError::InvalidFormat(
                "col_ptr must be monotonic".to_string(),
            ));
        }
        prev_col_ptr = v;
    }
    if prev_col_ptr as usize != nnz {
        return Err(CacheError::InvalidFormat(
            "col_ptr[n_cells] must equal nnz".to_string(),
        ));
    }

    for cell in 0..n_cells {
        let start = read_u64_slice(&mmap[col_ptr_offset + cell * 8..col_ptr_offset + cell * 8 + 8])
            as usize;
        let end = read_u64_slice(
            &mmap[col_ptr_offset + (cell + 1) * 8..col_ptr_offset + (cell + 1) * 8 + 8],
        ) as usize;
        let mut prev_row: Option<u32> = None;
        for i in start..end {
            let row = read_u32_slice(&mmap[row_idx_offset + i * 4..row_idx_offset + i * 4 + 4]);
            if row as usize >= n_genes {
                return Err(CacheError::InvalidFormat(
                    "row_idx out of bounds".to_string(),
                ));
            }
            if let Some(prev) = prev_row {
                if row <= prev {
                    return Err(CacheError::InvalidFormat(
                        "row_idx must be strictly increasing per column".to_string(),
                    ));
                }
            }
            prev_row = Some(row);
        }
    }

    Ok(())
}

fn read_u16_slice(slice: &[u8]) -> u16 {
    let mut buf = [0u8; 2];
    buf.copy_from_slice(slice);
    u16::from_le_bytes(buf)
}

fn read_u32_slice(slice: &[u8]) -> u32 {
    let mut buf = [0u8; 4];
    buf.copy_from_slice(slice);
    u32::from_le_bytes(buf)
}

fn read_u64_slice(slice: &[u8]) -> u64 {
    let mut buf = [0u8; 8];
    buf.copy_from_slice(slice);
    u64::from_le_bytes(buf)
}

pub fn write_expr_cache(
    path: &Path,
    expr: &ExprCsc,
    stats: &[CellStats],
) -> Result<(), CacheError> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    writer.write_all(MAGIC_EXPR)?;
    writer.write_all(&VERSION_EXPR.to_le_bytes())?;
    writer.write_all(&(expr.n_genes as u64).to_le_bytes())?;
    writer.write_all(&(expr.n_cells as u64).to_le_bytes())?;
    writer.write_all(&(expr.nnz as u64).to_le_bytes())?;
    writer.write_all(&(expr.col_ptr.len() as u64).to_le_bytes())?;
    writer.write_all(&(expr.row_idx.len() as u64).to_le_bytes())?;
    writer.write_all(&(expr.values.len() as u64).to_le_bytes())?;
    writer.write_all(&(stats.len() as u64).to_le_bytes())?;

    for v in &expr.col_ptr {
        writer.write_all(&v.to_le_bytes())?;
    }
    for v in &expr.row_idx {
        writer.write_all(&v.to_le_bytes())?;
    }
    for v in &expr.values {
        writer.write_all(&v.to_le_bytes())?;
    }
    for cell in stats {
        writer.write_all(&cell.libsize.to_le_bytes())?;
        writer.write_all(&cell.detected.to_le_bytes())?;
        writer.write_all(&0u32.to_le_bytes())?;
    }

    writer.flush()?;
    Ok(())
}

pub fn read_expr_cache(path: &Path) -> Result<(ExprCsc, Vec<CellStats>), CacheError> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);

    let mut magic = [0u8; 8];
    reader.read_exact(&mut magic)?;
    if &magic != MAGIC_EXPR {
        return Err(CacheError::InvalidMagic);
    }

    let version = read_u32(&mut reader)?;
    if version != VERSION_EXPR {
        return Err(CacheError::UnsupportedVersion(version));
    }

    let n_genes = read_u64(&mut reader)? as usize;
    let n_cells = read_u64(&mut reader)? as usize;
    let nnz = read_u64(&mut reader)? as usize;
    let col_len = read_u64(&mut reader)? as usize;
    let row_len = read_u64(&mut reader)? as usize;
    let val_len = read_u64(&mut reader)? as usize;
    let stats_len = read_u64(&mut reader)? as usize;

    if col_len != n_cells + 1 || row_len != nnz || val_len != nnz || stats_len != n_cells {
        return Err(CacheError::InvalidFormat(
            "lengths do not match header".to_string(),
        ));
    }

    let mut col_ptr = vec![0u64; col_len];
    for v in &mut col_ptr {
        *v = read_u64(&mut reader)?;
    }
    let mut row_idx = vec![0u32; row_len];
    for v in &mut row_idx {
        *v = read_u32(&mut reader)?;
    }
    let mut values = vec![0u32; val_len];
    for v in &mut values {
        *v = read_u32(&mut reader)?;
    }

    let mut stats = vec![CellStats::default(); stats_len];
    for cell in &mut stats {
        cell.libsize = read_u64(&mut reader)?;
        cell.detected = read_u32(&mut reader)?;
        let _ = read_u32(&mut reader)?;
    }

    Ok((
        ExprCsc {
            n_genes,
            n_cells,
            nnz,
            col_ptr,
            row_idx,
            values,
        },
        stats,
    ))
}

fn read_u32(reader: &mut dyn Read) -> Result<u32, std::io::Error> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

fn read_u64(reader: &mut dyn Read) -> Result<u64, std::io::Error> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

#[cfg(test)]
#[path = "../../tests/src_inline/input/cache.rs"]
mod tests;
