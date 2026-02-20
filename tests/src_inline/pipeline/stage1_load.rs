use super::*;
use crc::{CRC_64_ECMA_182, Crc};
use std::fs;
use tempfile::tempdir;

const HEADER_SIZE: usize = 256;
const CRC64: Crc<u64> = Crc::<u64>::new(&CRC_64_ECMA_182);

fn write_file(path: &Path, contents: &str) {
    fs::write(path, contents).expect("write file");
}

fn align64(x: usize) -> usize {
    (x + 63) & !63
}

fn encode_string_table(values: &[&str]) -> Vec<u8> {
    let mut blob = Vec::new();
    let mut offsets = Vec::with_capacity(values.len() + 1);
    offsets.push(0u32);
    for s in values {
        blob.extend_from_slice(s.as_bytes());
        offsets.push(blob.len() as u32);
    }
    let mut out = Vec::new();
    out.extend_from_slice(&(values.len() as u32).to_le_bytes());
    for off in offsets {
        out.extend_from_slice(&off.to_le_bytes());
    }
    out.extend_from_slice(&blob);
    out
}

fn write_shared_cache(path: &Path) {
    let genes = ["G1", "G2"];
    let barcodes = ["c1", "c2"];
    let col_ptr = [0u64, 1, 2];
    let row_idx = [0u32, 1];
    let values = [3u32, 4];

    let genes_table = encode_string_table(&genes);
    let barcodes_table = encode_string_table(&barcodes);

    let mut offset = HEADER_SIZE;
    let genes_off = align64(offset);
    offset = genes_off + genes_table.len();
    let barcodes_off = align64(offset);
    offset = barcodes_off + barcodes_table.len();
    let col_ptr_off = align64(offset);
    offset = col_ptr_off + col_ptr.len() * 8;
    let row_idx_off = align64(offset);
    offset = row_idx_off + row_idx.len() * 4;
    let values_off = align64(offset);
    offset = values_off + values.len() * 4;
    let file_bytes = offset;

    let mut out = vec![0u8; file_bytes];
    out[0..4].copy_from_slice(b"KORG");
    out[4..6].copy_from_slice(&1u16.to_le_bytes());
    out[6..8].copy_from_slice(&0u16.to_le_bytes());
    out[8..12].copy_from_slice(&0x1234_5678u32.to_le_bytes());
    out[12..16].copy_from_slice(&(HEADER_SIZE as u32).to_le_bytes());
    out[16..24].copy_from_slice(&(genes.len() as u64).to_le_bytes());
    out[24..32].copy_from_slice(&(barcodes.len() as u64).to_le_bytes());
    out[32..40].copy_from_slice(&(values.len() as u64).to_le_bytes());
    out[40..48].copy_from_slice(&(genes_off as u64).to_le_bytes());
    out[48..56].copy_from_slice(&(genes_table.len() as u64).to_le_bytes());
    out[56..64].copy_from_slice(&(barcodes_off as u64).to_le_bytes());
    out[64..72].copy_from_slice(&(barcodes_table.len() as u64).to_le_bytes());
    out[72..80].copy_from_slice(&(col_ptr_off as u64).to_le_bytes());
    out[80..88].copy_from_slice(&(row_idx_off as u64).to_le_bytes());
    out[88..96].copy_from_slice(&(values_off as u64).to_le_bytes());
    out[96..104].copy_from_slice(&0u64.to_le_bytes());
    out[104..112].copy_from_slice(&0u64.to_le_bytes());
    out[112..120].copy_from_slice(&(file_bytes as u64).to_le_bytes());
    out[128..136].copy_from_slice(&0u64.to_le_bytes());

    let mut hdr = out[0..HEADER_SIZE].to_vec();
    hdr[120..128].fill(0);
    let crc = CRC64.checksum(&hdr);
    out[120..128].copy_from_slice(&crc.to_le_bytes());

    out[genes_off..genes_off + genes_table.len()].copy_from_slice(&genes_table);
    out[barcodes_off..barcodes_off + barcodes_table.len()].copy_from_slice(&barcodes_table);
    for (i, v) in col_ptr.iter().enumerate() {
        let base = col_ptr_off + i * 8;
        out[base..base + 8].copy_from_slice(&v.to_le_bytes());
    }
    for (i, v) in row_idx.iter().enumerate() {
        let base = row_idx_off + i * 4;
        out[base..base + 4].copy_from_slice(&v.to_le_bytes());
    }
    for (i, v) in values.iter().enumerate() {
        let base = values_off + i * 4;
        out[base..base + 4].copy_from_slice(&v.to_le_bytes());
    }
    fs::write(path, out).expect("write shared cache");
}

#[test]
fn stage1_v3_ok() {
    let dir = tempdir().expect("tempdir");
    write_file(&dir.path().join("features.tsv"), "f1\tG1\nf2\tG2\n");
    write_file(&dir.path().join("barcodes.tsv"), "c1\nc2\nc3\n");
    write_file(
        &dir.path().join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n2 3 3\n1 1 1\n1 2 1\n2 3 1\n",
    );

    let ctx = run_stage1(
        dir.path(),
        None,
        dir.path(),
        false,
        RunMode::Standalone,
        None,
    )
    .expect("stage1 ok");
    assert_eq!(ctx.format, TenXFormat::TenXv3);
    assert_eq!(ctx.n_genes, 2);
    assert_eq!(ctx.n_cells, 3);
    assert_eq!(ctx.nnz, 3);
}

#[test]
fn stage1_dim_mismatch() {
    let dir = tempdir().expect("tempdir");
    write_file(&dir.path().join("features.tsv"), "f1\tG1\nf2\tG2\n");
    write_file(&dir.path().join("barcodes.tsv"), "c1\nc2\nc3\n");
    write_file(
        &dir.path().join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n2 4 3\n1 1 1\n1 2 1\n2 3 1\n",
    );

    let err = run_stage1(
        dir.path(),
        None,
        dir.path(),
        true,
        RunMode::Standalone,
        None,
    )
    .unwrap_err();
    match err {
        Stage1Error::DimensionMismatch { .. } => {}
        other => panic!("unexpected error: {other:?}"),
    }
}

#[test]
fn stage1_duplicate_gene_symbols() {
    let dir = tempdir().expect("tempdir");
    write_file(&dir.path().join("features.tsv"), "f1\tG1\nf2\tG1\nf3\tG2\n");
    write_file(&dir.path().join("barcodes.tsv"), "c1\nc2\n");
    write_file(
        &dir.path().join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n3 2 1\n1 1 1\n",
    );

    let ctx = run_stage1(
        dir.path(),
        None,
        dir.path(),
        true,
        RunMode::Standalone,
        None,
    )
    .expect("stage1 ok");
    assert_eq!(ctx.duplicate_gene_symbols_count, 1);
    assert_eq!(ctx.duplicate_gene_symbols[0].symbol, "G1");
    assert_eq!(ctx.duplicate_gene_symbols[0].first_row, 1);
    assert_eq!(ctx.duplicate_gene_symbols[0].dup_row, 2);
}

#[test]
fn stage1_meta_validation() {
    let dir = tempdir().expect("tempdir");
    write_file(&dir.path().join("features.tsv"), "f1\tG1\n");
    write_file(&dir.path().join("barcodes.tsv"), "c1\nc2\n");
    write_file(
        &dir.path().join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n1 2 1\n1 1 1\n",
    );
    write_file(
        &dir.path().join("meta.tsv"),
        "cell_id\tsample_id\nc1\ts1\nc1\ts1\nmissing\ts2\n",
    );

    let ctx = run_stage1(
        dir.path(),
        Some(&dir.path().join("meta.tsv")),
        dir.path(),
        true,
        RunMode::Standalone,
        None,
    )
    .expect("stage1 ok");
    assert!(ctx.meta_present);
    assert_eq!(ctx.meta_cells_matched, 1);
    assert_eq!(ctx.meta_cells_missing, 1);
}

#[test]
fn pipeline_mode_uses_cache_when_present() {
    let dir = tempdir().expect("tempdir");
    let cache = dir.path().join("kira-organelle.bin");
    write_shared_cache(&cache);

    let ctx = run_stage1(dir.path(), None, dir.path(), true, RunMode::Pipeline, None).expect("ctx");
    assert_eq!(ctx.shared_cache_path, Some(cache.clone()));
    assert_eq!(ctx.resolved_shared_cache_path, Some(cache));
    assert_eq!(ctx.n_genes, 2);
    assert_eq!(ctx.n_cells, 2);
    assert_eq!(ctx.nnz, 2);
}

#[test]
fn pipeline_mode_uses_suffix_cache_when_present() {
    let dir = tempdir().expect("tempdir");
    let cache = dir.path().join("GSM1.kira-organelle.bin");
    write_shared_cache(&cache);

    let ctx = run_stage1(dir.path(), None, dir.path(), true, RunMode::Pipeline, None).expect("ctx");
    assert_eq!(ctx.shared_cache_path, Some(cache.clone()));
    assert_eq!(ctx.resolved_shared_cache_path, Some(cache));
    assert_eq!(ctx.n_genes, 2);
    assert_eq!(ctx.n_cells, 2);
    assert_eq!(ctx.nnz, 2);
}

#[test]
fn pipeline_mode_falls_back_when_cache_missing() {
    let dir = tempdir().expect("tempdir");
    write_file(&dir.path().join("features.tsv"), "f1\tG1\n");
    write_file(&dir.path().join("barcodes.tsv"), "c1\n");
    write_file(
        &dir.path().join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n1 1 1\n1 1 1\n",
    );

    let ctx = run_stage1(dir.path(), None, dir.path(), true, RunMode::Pipeline, None).expect("ctx");
    assert!(ctx.shared_cache_path.is_none());
    assert_eq!(
        ctx.resolved_shared_cache_path,
        Some(dir.path().join("kira-organelle.bin"))
    );
}

#[test]
fn pipeline_mode_invalid_cache_hard_fails() {
    let dir = tempdir().expect("tempdir");
    fs::write(dir.path().join("kira-organelle.bin"), b"bad").expect("write");
    let err = run_stage1(dir.path(), None, dir.path(), true, RunMode::Pipeline, None).unwrap_err();
    match err {
        Stage1Error::Cache(_) => {}
        other => panic!("unexpected error: {other:?}"),
    }
}
