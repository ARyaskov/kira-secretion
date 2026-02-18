use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum Regime {
    SelfPreserving,
    EnvironmentShaping,
    ExportDominant,
    SecretoryLysosomeActive,
    MetabolicSuppressive,
    InflammatorySignaler,
    PresentationHigh,
    Unclassified,
}

impl Regime {
    pub fn as_str(&self) -> &'static str {
        match self {
            Regime::SelfPreserving => "SelfPreserving",
            Regime::EnvironmentShaping => "EnvironmentShaping",
            Regime::ExportDominant => "ExportDominant",
            Regime::SecretoryLysosomeActive => "SecretoryLysosomeActive",
            Regime::MetabolicSuppressive => "MetabolicSuppressive",
            Regime::InflammatorySignaler => "InflammatorySignaler",
            Regime::PresentationHigh => "PresentationHigh",
            Regime::Unclassified => "Unclassified",
        }
    }

    pub fn ordered() -> &'static [Regime] {
        &[
            Regime::SelfPreserving,
            Regime::EnvironmentShaping,
            Regime::ExportDominant,
            Regime::SecretoryLysosomeActive,
            Regime::MetabolicSuppressive,
            Regime::InflammatorySignaler,
            Regime::PresentationHigh,
            Regime::Unclassified,
        ]
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum RuleId {
    R1SelfPreserving,
    R2SecretoryLysosomeActive,
    R3ExportDominant,
    R4MetabolicSuppressive,
    R5InflammatorySignaler,
    R6PresentationHigh,
    R7EnvironmentShaping,
    R0Unclassified,
}

impl RuleId {
    pub fn as_str(&self) -> &'static str {
        match self {
            RuleId::R1SelfPreserving => "R1_SELF_PRESERVING",
            RuleId::R2SecretoryLysosomeActive => "R2_SECRETORY_LYSOSOME_ACTIVE",
            RuleId::R3ExportDominant => "R3_EXPORT_DOMINANT",
            RuleId::R4MetabolicSuppressive => "R4_METABOLIC_SUPPRESSIVE",
            RuleId::R5InflammatorySignaler => "R5_INFLAMMATORY_SIGNALER",
            RuleId::R6PresentationHigh => "R6_PRESENTATION_HIGH",
            RuleId::R7EnvironmentShaping => "R7_ENVIRONMENT_SHAPING",
            RuleId::R0Unclassified => "R0_UNCLASSIFIED",
        }
    }
}
