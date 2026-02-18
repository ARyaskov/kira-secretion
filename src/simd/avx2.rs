#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
pub fn sum_u32(values: &[u32]) -> u64 {
    // SAFETY: this function is compiled only when target includes AVX2.
    unsafe { sum_u32_avx2(values) }
}

#[cfg(not(all(target_arch = "x86_64", target_feature = "avx2")))]
pub fn sum_u32(values: &[u32]) -> u64 {
    values.iter().map(|v| *v as u64).sum()
}

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
#[target_feature(enable = "avx2")]
unsafe fn sum_u32_avx2(values: &[u32]) -> u64 {
    let mut i = 0usize;
    let len = values.len();
    let mut acc_lo = _mm256_setzero_si256();
    let mut acc_hi = _mm256_setzero_si256();

    while i + 8 <= len {
        let ptr = values.as_ptr().add(i) as *const __m256i;
        let v = _mm256_loadu_si256(ptr);

        let lo128 = _mm256_castsi256_si128(v);
        let hi128 = _mm256_extracti128_si256(v, 1);

        let lo64 = _mm256_cvtepu32_epi64(lo128);
        let hi64 = _mm256_cvtepu32_epi64(hi128);

        acc_lo = _mm256_add_epi64(acc_lo, lo64);
        acc_hi = _mm256_add_epi64(acc_hi, hi64);

        i += 8;
    }

    let mut buf_lo = [0u64; 4];
    let mut buf_hi = [0u64; 4];
    _mm256_storeu_si256(buf_lo.as_mut_ptr() as *mut __m256i, acc_lo);
    _mm256_storeu_si256(buf_hi.as_mut_ptr() as *mut __m256i, acc_hi);

    let mut sum = buf_lo.iter().copied().sum::<u64>() + buf_hi.iter().copied().sum::<u64>();

    while i < len {
        sum += values[i] as u64;
        i += 1;
    }
    sum
}
