
use super::*;
use std::collections::HashMap;

#[test]
fn mapping_missing_required() {
    let mut index = GeneIndex {
        rows: Vec::new(),
        duplicates: Vec::new(),
        first_index_by_symbol: HashMap::new(),
    };
    index.first_index_by_symbol.insert("A".to_string(), 1);
    index.first_index_by_symbol.insert("B".to_string(), 2);

    let panel = PanelDef {
        id: "P1".to_string(),
        description: "".to_string(),
        axis: "X".to_string(),
        genes: vec![
            crate::panels::defs::PanelGene {
                symbol: "A".to_string(),
            },
            crate::panels::defs::PanelGene {
                symbol: "C".to_string(),
            },
        ],
        required: vec!["A".to_string(), "C".to_string()],
        weights: None,
    };

    let (mapping, warning) = map_panel(&panel, &index);
    assert_eq!(mapping.mapped.len(), 2);
    assert_eq!(mapping.required_hits, 1);
    assert_eq!(mapping.required_total, 2);
    assert!(warning.is_some());
    assert_eq!(warning.unwrap().missing_required, vec!["C".to_string()]);
}
