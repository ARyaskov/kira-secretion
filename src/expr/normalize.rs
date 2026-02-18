#[derive(Debug, Clone)]
pub struct Normalization {
    pub enabled: bool,
    pub scale: f32,
    pub epsilon: f32,
}

impl Default for Normalization {
    fn default() -> Self {
        Self {
            enabled: true,
            scale: 10_000.0,
            epsilon: 1e-8,
        }
    }
}
