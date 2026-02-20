use super::*;
use std::fs;
use tempfile::tempdir;

#[test]
fn parse_header_and_dims() {
    let dir = tempdir().expect("tempdir");
    let path = dir.path().join("matrix.mtx");
    fs::write(
        &path,
        "%%MatrixMarket matrix coordinate integer general\n% comment\n2 3 4\n1 1 1\n",
    )
    .expect("write file");

    let header = read_header(&path).expect("read header");
    assert_eq!(header.n_rows, 2);
    assert_eq!(header.n_cols, 3);
    assert_eq!(header.nnz, 4);
}
