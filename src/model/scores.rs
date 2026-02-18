#[derive(Debug, Clone, Copy)]
pub struct WeightsDefault {
    pub oii: OiiWeights,
    pub iai_with_apci: IaiWeights,
    pub iai_no_apci: IaiNoApciWeights,
    pub esi: EsiWeights,
}

#[derive(Debug, Clone, Copy)]
pub struct OiiWeights {
    pub sia: f32,
    pub pos_eeb: f32,
    pub sli: f32,
    pub mei: f32,
    pub ecmi: f32,
    pub gdi: f32,
}

#[derive(Debug, Clone, Copy)]
pub struct IaiWeights {
    pub mei: f32,
    pub gdi: f32,
    pub apci: f32,
    pub sia: f32,
    pub pos_eeb: f32,
}

#[derive(Debug, Clone, Copy)]
pub struct IaiNoApciWeights {
    pub mei: f32,
    pub gdi: f32,
    pub sia: f32,
    pub pos_eeb: f32,
}

#[derive(Debug, Clone, Copy)]
pub struct EsiWeights {
    pub ecmi: f32,
    pub mei: f32,
    pub pos_eeb: f32,
    pub sli: f32,
}

impl Default for WeightsDefault {
    fn default() -> Self {
        Self {
            oii: OiiWeights {
                sia: 0.22,
                pos_eeb: 0.18,
                sli: 0.12,
                mei: 0.16,
                ecmi: 0.16,
                gdi: 0.16,
            },
            iai_with_apci: IaiWeights {
                mei: 0.22,
                gdi: 0.22,
                apci: 0.26,
                sia: 0.18,
                pos_eeb: 0.12,
            },
            iai_no_apci: IaiNoApciWeights {
                mei: 0.30,
                gdi: 0.30,
                sia: 0.25,
                pos_eeb: 0.15,
            },
            esi: EsiWeights {
                ecmi: 0.34,
                mei: 0.26,
                pos_eeb: 0.20,
                sli: 0.20,
            },
        }
    }
}

pub fn clamp01(x: f32) -> f32 {
    if x.is_nan() {
        0.0
    } else if x < 0.0 {
        0.0
    } else if x > 1.0 {
        1.0
    } else {
        x
    }
}

pub fn pos_eeb(eeb: f32) -> f32 {
    (eeb + 1.0) * 0.5
}
