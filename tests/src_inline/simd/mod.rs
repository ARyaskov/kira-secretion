use super::*;

#[test]
fn backend_name_is_supported() {
    assert!(matches!(backend_name(), "scalar" | "avx2" | "neon"));
}

#[test]
fn sum_u32_matches_scalar() {
    let data = vec![1u32, 3, 7, 11, 13, 17, 23, 31, 37, 41, 43];
    let expected: u64 = data.iter().map(|v| *v as u64).sum();
    assert_eq!(sum_u32(&data), expected);
}
