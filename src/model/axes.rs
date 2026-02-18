#[derive(Debug, Clone, Copy)]
pub struct AxisConfig {
    pub k: f32,
    pub epsilon: f32,
}

impl Default for AxisConfig {
    fn default() -> Self {
        Self {
            k: 1.0,
            epsilon: 1e-8,
        }
    }
}

pub fn saturating_map(x: f32, k: f32) -> f32 {
    if x <= 0.0 { 0.0 } else { x / (x + k) }
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct AxisValues {
    pub sia: f32,
    pub eeb: f32,
    pub sli: f32,
    pub mei: f32,
    pub ecmi: f32,
    pub apci: f32,
    pub gdi: f32,
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct AxisCoverage {
    pub sia: f32,
    pub eeb: f32,
    pub sli: f32,
    pub mei: f32,
    pub ecmi: f32,
    pub apci: f32,
    pub gdi: f32,
}
