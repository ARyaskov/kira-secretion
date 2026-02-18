use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(transparent)]
pub struct PanelGene {
    pub symbol: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PanelDef {
    pub id: String,
    pub description: String,
    pub axis: String,
    #[serde(default)]
    pub genes: Vec<PanelGene>,
    #[serde(default)]
    pub required: Vec<String>,
    #[serde(default)]
    pub weights: Option<Vec<f32>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PanelSet {
    #[serde(default)]
    pub panels: Vec<PanelDef>,
}

impl PanelDef {
    pub fn gene_symbols(&self) -> impl Iterator<Item = &str> {
        self.genes.iter().map(|g| g.symbol.as_str())
    }
}
