
use super::*;
use tempfile::tempdir;

#[test]
fn detects_prefix_none() {
    let dir = tempdir().expect("tempdir");
    std::fs::write(dir.path().join("matrix.mtx"), "x").expect("write");
    let got = detect_prefix(dir.path()).expect("prefix");
    assert_eq!(got, None);
}

#[test]
fn detects_prefix_present() {
    let dir = tempdir().expect("tempdir");
    std::fs::write(dir.path().join("ABC_matrix.mtx"), "x").expect("write");
    let got = detect_prefix(dir.path()).expect("prefix");
    assert_eq!(got, Some("ABC".to_string()));
}

#[test]
fn resolves_cache_name() {
    assert_eq!(resolve_shared_cache_file_name(None), "kira-organelle.bin");
    assert_eq!(
        resolve_shared_cache_file_name(Some("XYZ")),
        "XYZ.kira-organelle.bin"
    );
}

#[test]
fn finds_cache_by_suffix_when_exact_missing() {
    let dir = tempdir().expect("tempdir");
    std::fs::write(dir.path().join("ABC.kira-organelle.bin"), "x").expect("write");
    let got = find_shared_cache_file(dir.path(), None).expect("find");
    assert_eq!(got, Some(dir.path().join("ABC.kira-organelle.bin")));
}

#[test]
fn finds_exact_cache_preferred() {
    let dir = tempdir().expect("tempdir");
    std::fs::write(dir.path().join("ABC.kira-organelle.bin"), "x").expect("write");
    std::fs::write(dir.path().join("kira-organelle.bin"), "x").expect("write");
    let got = find_shared_cache_file(dir.path(), None).expect("find");
    assert_eq!(got, Some(dir.path().join("kira-organelle.bin")));
}
