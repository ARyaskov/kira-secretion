#[derive(Debug, Clone, Copy)]
pub struct Thresholds {
    pub low_counts: u64,
    pub few_detected: u32,
    pub cov_min: f32,
    pub oii_hi: f32,
    pub esi_hi: f32,
    pub esi_very: f32,
    pub sia_low: f32,
    pub sia_mid: f32,
    pub sia_hi: f32,
    pub pos_eeb_hi: f32,
    pub pos_eeb_mid: f32,
    pub pos_eeb_low: f32,
    pub sli_hi: f32,
    pub mei_hi: f32,
    pub ecmi_hi: f32,
    pub gdi_hi: f32,
    pub apci_hi: f32,
    pub ambient_gdi: f32,
    pub ambient_sia: f32,
}

impl Default for Thresholds {
    fn default() -> Self {
        Self {
            low_counts: 500,
            few_detected: 300,
            cov_min: 0.60,
            oii_hi: 0.65,
            esi_hi: 0.65,
            esi_very: 0.75,
            sia_low: 0.35,
            sia_mid: 0.40,
            sia_hi: 0.55,
            pos_eeb_hi: 0.70,
            pos_eeb_mid: 0.55,
            pos_eeb_low: 0.45,
            sli_hi: 0.70,
            mei_hi: 0.70,
            ecmi_hi: 0.70,
            gdi_hi: 0.75,
            apci_hi: 0.70,
            ambient_gdi: 0.75,
            ambient_sia: 0.45,
        }
    }
}
