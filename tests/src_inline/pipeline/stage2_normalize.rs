
use super::*;
use crc::{CRC_64_ECMA_182, Crc};
use std::collections::HashMap;
use std::fs;
use tempfile::tempdir;

const HEADER_SIZE: usize = 256;
const CRC64: Crc<u64> = Crc::<u64>::new(&CRC_64_ECMA_182);

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
    let genes = ["G1"];
    let barcodes = ["c1"];
    let col_ptr = [0u64, 1];
    let row_idx = [0u32];
    let values = [5u32];

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
fn stage2_uses_shared_cache_when_present() {
    let dir = tempdir().expect("tempdir");
    let cache = dir.path().join("kira-organelle.bin");
    write_shared_cache(&cache);
    let ctx = DatasetCtx {
        format: crate::input::detect::TenXFormat::Unknown,
        matrix_path: dir.path().join("matrix.mtx"),
        features_path: dir.path().join("features.tsv"),
        barcodes_path: dir.path().join("barcodes.tsv"),
        shared_cache_path: Some(cache.clone()),
        resolved_shared_cache_path: Some(cache),
        gene_index: crate::input::features::GeneIndex {
            rows: Vec::new(),
            duplicates: Vec::new(),
            first_index_by_symbol: HashMap::new(),
        },
        barcodes: vec!["c1".to_string()],
        n_genes: 1,
        n_cells: 1,
        nnz: 1,
        duplicate_gene_symbols_count: 0,
        duplicate_gene_symbols: Vec::new(),
        meta_present: false,
        meta_cells_matched: 0,
        meta_cells_missing: 0,
    };

    let expr = run_stage2(&ctx, dir.path(), Normalization::default(), true).expect("stage2");
    match expr.expr {
        ExprMatrix::Shared(ref shared) => {
            assert_eq!(shared.n_genes, 1);
            assert_eq!(shared.n_cells, 1);
            assert_eq!(shared.nnz, 1);
            assert_eq!(shared.genes, vec!["G1"]);
            assert_eq!(shared.barcodes, vec!["c1"]);
            assert_eq!(shared.col_ptr_at(0), 0);
            assert_eq!(shared.col_ptr_at(1), 1);
            assert_eq!(shared.row_idx_at(0), 0);
            assert_eq!(shared.value_at(0), 5);
        }
        ExprMatrix::Owned(_) => panic!("expected shared cache expression"),
    }
}
