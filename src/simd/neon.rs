#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
use std::arch::aarch64::*;

#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
pub fn sum_u32(values: &[u32]) -> u64 {
    // SAFETY: function uses NEON intrinsics and is compiled for aarch64+neon.
    unsafe { sum_u32_neon(values) }
}

#[cfg(not(all(target_arch = "aarch64", target_feature = "neon")))]
pub fn sum_u32(values: &[u32]) -> u64 {
    values.iter().map(|v| *v as u64).sum()
}

#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
unsafe fn sum_u32_neon(values: &[u32]) -> u64 {
    let mut i = 0usize;
    let len = values.len();
    let mut acc = unsafe { vdupq_n_u64(0) };

    while i + 4 <= len {
        let v = unsafe { vld1q_u32(values.as_ptr().add(i)) };
        let wide = unsafe { vpaddlq_u32(v) };
        acc = unsafe { vaddq_u64(acc, wide) };
        i += 4;
    }

    let mut buf = [0u64; 2];
    unsafe { vst1q_u64(buf.as_mut_ptr(), acc) };
    let mut sum = buf[0] + buf[1];

    while i < len {
        sum += values[i] as u64;
        i += 1;
    }
    sum
}
