use super::*;
use crate::panels::defs::{PanelDef, PanelGene, PanelSet};
use crate::pipeline::stage3_panels::PanelsContext;
use std::collections::HashMap;
use std::fs;
use tempfile::tempdir;

fn make_panels_ctx() -> PanelsContext {
    let panels = PanelSet {
        panels: vec![
            PanelDef {
                id: "P_SIA".to_string(),
                description: "".to_string(),
                axis: "SIA".to_string(),
                genes: vec![PanelGene {
                    symbol: "A".to_string(),
                }],
                required: vec!["A".to_string()],
                weights: None,
            },
            PanelDef {
                id: "P_EXP".to_string(),
                description: "".to_string(),
                axis: "EEB_EXPORT".to_string(),
                genes: vec![PanelGene {
                    symbol: "B".to_string(),
                }],
                required: vec!["B".to_string()],
                weights: None,
            },
            PanelDef {
                id: "P_DEG".to_string(),
                description: "".to_string(),
                axis: "EEB_DEGRADE".to_string(),
                genes: vec![PanelGene {
                    symbol: "C".to_string(),
                }],
                required: vec!["C".to_string()],
                weights: None,
            },
        ],
    };
    let mut mappings = Vec::new();
    for panel in &panels.panels {
        mappings.push(crate::panels::mapping::GeneMapping {
            panel_id: panel.id.clone(),
            mapped: vec![Some(0)],
            required_hits: panel.required.len(),
            required_total: panel.required.len(),
        });
    }
    PanelsContext {
        panels,
        mappings,
        warnings: Vec::new(),
        cell_ids: vec!["c1".to_string()],
        per_cell: vec![PanelCellPacked {
            sums: vec![2.0, 3.0, 1.0],
            hits: vec![1, 1, 1],
            required_missing: vec![0, 0, 0],
        }],
    }
}

#[test]
fn axis_correctness() {
    let ctx = make_panels_ctx();
    let dir = tempdir().expect("tempdir");
    fs::create_dir_all(dir.path()).expect("mkdir");
    let dummy = DatasetCtx {
        format: crate::input::detect::TenXFormat::TenXv3,
        matrix_path: dir.path().join("matrix.mtx"),
        features_path: dir.path().join("features.tsv"),
        barcodes_path: dir.path().join("barcodes.tsv"),
        shared_cache_path: None,
        resolved_shared_cache_path: None,
        gene_index: crate::input::features::GeneIndex {
            rows: Vec::new(),
            duplicates: Vec::new(),
            first_index_by_symbol: HashMap::new(),
        },
        barcodes: vec!["c1".to_string()],
        n_genes: 3,
        n_cells: 1,
        nnz: 3,
        duplicate_gene_symbols_count: 0,
        duplicate_gene_symbols: Vec::new(),
        meta_present: false,
        meta_cells_matched: 0,
        meta_cells_missing: 0,
    };
    let axes = run_stage4_axes(&dummy, &ctx, dir.path()).expect("axes");
    let sia = axes.values[0].sia;
    let eeb = axes.values[0].eeb;
    let sia_expected = 2.0 / (2.0 + 1.0);
    assert!((sia - sia_expected).abs() < 1e-6);
    let eeb_expected = (3.0 - 1.0) / (3.0 + 1.0 + 1e-8);
    assert!((eeb - eeb_expected).abs() < 1e-6);
}

#[test]
fn driver_determinism() {
    let ids = vec!["B".to_string(), "A".to_string()];
    let vals = vec![1.0, 1.0];
    let drivers = top_k_panels(&ids, &vals, 2);
    assert_eq!(drivers[0].panel_id, "A");
}

#[test]
fn axes_determinism_report() {
    let ctx = make_panels_ctx();
    let dir = tempdir().expect("tempdir");
    fs::create_dir_all(dir.path()).expect("mkdir");
    let dummy = DatasetCtx {
        format: crate::input::detect::TenXFormat::TenXv3,
        matrix_path: dir.path().join("matrix.mtx"),
        features_path: dir.path().join("features.tsv"),
        barcodes_path: dir.path().join("barcodes.tsv"),
        shared_cache_path: None,
        resolved_shared_cache_path: None,
        gene_index: crate::input::features::GeneIndex {
            rows: Vec::new(),
            duplicates: Vec::new(),
            first_index_by_symbol: HashMap::new(),
        },
        barcodes: vec!["c1".to_string()],
        n_genes: 3,
        n_cells: 1,
        nnz: 3,
        duplicate_gene_symbols_count: 0,
        duplicate_gene_symbols: Vec::new(),
        meta_present: false,
        meta_cells_matched: 0,
        meta_cells_missing: 0,
    };
    let out1 = dir.path().join("out1");
    let out2 = dir.path().join("out2");
    fs::create_dir_all(&out1).expect("mkdir");
    fs::create_dir_all(&out2).expect("mkdir");
    run_stage4_axes(&dummy, &ctx, &out1).expect("axes1");
    run_stage4_axes(&dummy, &ctx, &out2).expect("axes2");
    let a = fs::read(out1.join("axes.tsv")).expect("read1");
    let b = fs::read(out2.join("axes.tsv")).expect("read2");
    assert_eq!(a, b);
}

#[test]
fn coverage_correctness() {
    let panels = PanelSet {
        panels: vec![PanelDef {
            id: "P1".to_string(),
            description: "".to_string(),
            axis: "SIA".to_string(),
            genes: vec![PanelGene {
                symbol: "A".to_string(),
            }],
            required: vec!["A".to_string(), "B".to_string()],
            weights: None,
        }],
    };
    let mappings = vec![crate::panels::mapping::GeneMapping {
        panel_id: "P1".to_string(),
        mapped: vec![Some(0)],
        required_hits: 1,
        required_total: 2,
    }];
    let ctx = PanelsContext {
        panels,
        mappings,
        warnings: Vec::new(),
        cell_ids: vec!["c1".to_string()],
        per_cell: vec![PanelCellPacked {
            sums: vec![1.0],
            hits: vec![1],
            required_missing: vec![1],
        }],
    };
    let indices = build_axis_indices(&ctx.panels);
    let (vals, cov, _) =
        compute_cell_axes(&indices, &ctx, &ctx.per_cell[0], &AxisConfig::default());
    assert!((vals.sia - 0.5).abs() < 1e-6);
    assert!((cov.sia - 0.5).abs() < 1e-6);
}
