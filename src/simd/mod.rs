pub mod avx2;
pub mod neon;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Backend {
    Scalar,
    Avx2,
    Neon,
}

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
pub const BACKEND: Backend = Backend::Avx2;

#[cfg(all(
    not(all(target_arch = "x86_64", target_feature = "avx2")),
    target_arch = "aarch64",
    target_feature = "neon"
))]
pub const BACKEND: Backend = Backend::Neon;

#[cfg(not(any(
    all(target_arch = "x86_64", target_feature = "avx2"),
    all(target_arch = "aarch64", target_feature = "neon")
)))]
pub const BACKEND: Backend = Backend::Scalar;

pub fn backend_name() -> &'static str {
    match BACKEND {
        Backend::Scalar => "scalar",
        Backend::Avx2 => "avx2",
        Backend::Neon => "neon",
    }
}

pub fn sum_u32(values: &[u32]) -> u64 {
    #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
    {
        return avx2::sum_u32(values);
    }

    #[cfg(all(
        not(all(target_arch = "x86_64", target_feature = "avx2")),
        target_arch = "aarch64",
        target_feature = "neon"
    ))]
    {
        return neon::sum_u32(values);
    }

    #[cfg(not(any(
        all(target_arch = "x86_64", target_feature = "avx2"),
        all(target_arch = "aarch64", target_feature = "neon")
    )))]
    {
        values.iter().map(|v| *v as u64).sum()
    }
}

#[cfg(test)]
#[path = "../../tests/src_inline/simd/mod.rs"]
mod tests;
