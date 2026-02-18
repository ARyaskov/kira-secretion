use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Flags {
    bits: u8,
}

impl Flags {
    pub const LOW_CONFIDENCE: u8 = 0b0001;
    pub const FEW_DETECTED_GENES: u8 = 0b0010;
    pub const LOW_COUNTS: u8 = 0b0100;
    pub const HIGH_AMBIENT_RISK: u8 = 0b1000;

    pub fn empty() -> Self {
        Self { bits: 0 }
    }

    pub fn set(&mut self, bit: u8) {
        self.bits |= bit;
    }

    pub fn contains(&self, bit: u8) -> bool {
        self.bits & bit != 0
    }

    pub fn to_csv(&self) -> String {
        if self.bits == 0 {
            return ".".to_string();
        }
        let mut parts = Vec::new();
        if self.contains(Self::LOW_CONFIDENCE) {
            parts.push("LOW_CONFIDENCE");
        }
        if self.contains(Self::FEW_DETECTED_GENES) {
            parts.push("FEW_DETECTED_GENES");
        }
        if self.contains(Self::LOW_COUNTS) {
            parts.push("LOW_COUNTS");
        }
        if self.contains(Self::HIGH_AMBIENT_RISK) {
            parts.push("HIGH_AMBIENT_RISK");
        }
        parts.join(",")
    }
}
