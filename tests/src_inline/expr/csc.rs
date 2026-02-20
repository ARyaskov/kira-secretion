use super::*;
use std::fs;
use tempfile::tempdir;

#[test]
fn build_csc_tiny() {
    let dir = tempdir().expect("tempdir");
    let path = dir.path().join("matrix.mtx");
    fs::write(
        &path,
        "%%MatrixMarket matrix coordinate integer general\n3 2 4\n1 1 1\n2 1 2\n3 2 3\n1 2 4\n",
    )
    .expect("write file");

    let (csc, stats) = ExprCsc::from_mtx(&path, 3, 2, false).expect("csc");
    assert_eq!(csc.col_ptr, vec![0, 2, 4]);
    assert_eq!(csc.row_idx, vec![0, 1, 0, 2]);
    assert_eq!(csc.values, vec![1, 2, 4, 3]);
    assert_eq!(stats[0].libsize, 3);
    assert_eq!(stats[0].detected, 2);
    assert_eq!(stats[1].libsize, 7);
    assert_eq!(stats[1].detected, 2);
}

#[test]
fn normalization_values() {
    let dir = tempdir().expect("tempdir");
    let path = dir.path().join("matrix.mtx");
    fs::write(
        &path,
        "%%MatrixMarket matrix coordinate integer general\n2 1 2\n1 1 1\n2 1 3\n",
    )
    .expect("write file");

    let (csc, stats) = ExprCsc::from_mtx(&path, 2, 1, false).expect("csc");
    let norm = Normalization::default();
    let values: Vec<(u32, f32)> = csc.iter_cell_norm(0, &norm, &stats[0]).collect();
    let denom = stats[0].libsize as f32 + norm.epsilon;
    let v0 = (1.0 * (norm.scale / denom)).ln_1p();
    let v1 = (3.0 * (norm.scale / denom)).ln_1p();
    assert_eq!(values.len(), 2);
    assert!((values[0].1 - v0).abs() < 1e-6);
    assert!((values[1].1 - v1).abs() < 1e-6);
}

#[test]
fn determinism_repeat_build() {
    let dir = tempdir().expect("tempdir");
    let path = dir.path().join("matrix.mtx");
    fs::write(
        &path,
        "%%MatrixMarket matrix coordinate integer general\n3 2 3\n1 1 2\n2 2 3\n3 1 4\n",
    )
    .expect("write file");

    let (csc1, stats1) = ExprCsc::from_mtx(&path, 3, 2, false).expect("csc1");
    let (csc2, stats2) = ExprCsc::from_mtx(&path, 3, 2, false).expect("csc2");

    assert_eq!(csc1.col_ptr, csc2.col_ptr);
    assert_eq!(csc1.row_idx, csc2.row_idx);
    assert_eq!(csc1.values, csc2.values);
    assert_eq!(stats1[0].libsize, stats2[0].libsize);
    assert_eq!(stats1[1].detected, stats2[1].detected);
}
