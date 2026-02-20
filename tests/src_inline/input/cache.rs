use super::*;
use std::fs;
use tempfile::tempdir;

fn write_shared_cache(path: &Path, tamper_crc: bool) {
    let genes = ["G1", "G2", "G3"];
    let barcodes = ["C1", "C2"];
    let col_ptr = [0u64, 2, 3];
    let row_idx = [0u32, 2, 1];
    let values = [5u32, 1, 7];

    let genes_table = encode_string_table(&genes);
    let barcodes_table = encode_string_table(&barcodes);

    let mut offset = SHARED_HEADER_SIZE;
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

    out[0..4].copy_from_slice(SHARED_MAGIC);
    out[4..6].copy_from_slice(&1u16.to_le_bytes());
    out[6..8].copy_from_slice(&0u16.to_le_bytes());
    out[8..12].copy_from_slice(&SHARED_ENDIAN_TAG.to_le_bytes());
    out[12..16].copy_from_slice(&(SHARED_HEADER_SIZE as u32).to_le_bytes());
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

    let mut header_for_crc = out[0..SHARED_HEADER_SIZE].to_vec();
    header_for_crc[120..128].fill(0);
    let mut crc = CRC64.checksum(&header_for_crc);
    if tamper_crc {
        crc ^= 0xFFFF;
    }
    out[120..128].copy_from_slice(&crc.to_le_bytes());
    out[128..136].copy_from_slice(&0u64.to_le_bytes());

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

#[test]
fn shared_cache_valid() {
    let dir = tempdir().expect("tempdir");
    let path = dir.path().join("kira-organelle.bin");
    write_shared_cache(&path, false);
    let mapped = mmap_shared_cache(&path).expect("shared cache");
    assert_eq!(mapped.n_genes, 3);
    assert_eq!(mapped.n_cells, 2);
    assert_eq!(mapped.nnz, 3);
    assert_eq!(mapped.genes, vec!["G1", "G2", "G3"]);
    assert_eq!(mapped.barcodes, vec!["C1", "C2"]);
    assert_eq!(mapped.col_ptr_at(2), 3);
    assert_eq!(mapped.row_idx_at(2), 1);
    assert_eq!(mapped.value_at(2), 7);
}

#[test]
fn shared_cache_bad_crc_rejected() {
    let dir = tempdir().expect("tempdir");
    let path = dir.path().join("kira-organelle.bin");
    write_shared_cache(&path, true);
    let err = mmap_shared_cache(&path).expect_err("expected error");
    assert!(format!("{err}").contains("CRC64"));
}

#[test]
fn cache_roundtrip_deterministic() {
    let dir = tempdir().expect("tempdir");
    let path = dir.path().join("expr.bin");

    let expr = ExprCsc {
        n_genes: 3,
        n_cells: 2,
        nnz: 3,
        col_ptr: vec![0, 2, 3],
        row_idx: vec![0, 2, 1],
        values: vec![1, 2, 3],
    };
    let stats = vec![
        CellStats {
            libsize: 3,
            detected: 2,
        },
        CellStats {
            libsize: 3,
            detected: 1,
        },
    ];

    write_expr_cache(&path, &expr, &stats).expect("write cache");
    let bytes_first = fs::read(&path).expect("read cache");

    write_expr_cache(&path, &expr, &stats).expect("write cache");
    let bytes_second = fs::read(&path).expect("read cache");

    assert_eq!(bytes_first, bytes_second);

    let (expr2, stats2) = read_expr_cache(&path).expect("read back");
    assert_eq!(expr2.col_ptr, expr.col_ptr);
    assert_eq!(expr2.row_idx, expr.row_idx);
    assert_eq!(expr2.values, expr.values);
    assert_eq!(stats2.len(), stats.len());
    assert_eq!(stats2[0].libsize, stats[0].libsize);
    assert_eq!(stats2[0].detected, stats[0].detected);
}
