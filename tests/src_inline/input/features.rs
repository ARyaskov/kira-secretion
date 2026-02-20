use super::*;
use std::fs;
use tempfile::tempdir;

#[test]
fn duplicate_gene_symbols_first_wins() {
    let dir = tempdir().expect("tempdir");
    let path = dir.path().join("features.tsv");
    fs::write(&path, "f1\tG1\nf2\tG1\nf3\tG2\n").expect("write file");

    let index = read_features(&path).expect("read features");
    assert_eq!(index.rows.len(), 3);
    assert_eq!(index.duplicates.len(), 1);
    assert_eq!(index.duplicates[0].symbol, "G1");
    assert_eq!(index.duplicates[0].first_row, 1);
    assert_eq!(index.duplicates[0].dup_row, 2);
}
