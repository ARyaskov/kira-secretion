use super::*;
use crate::expr::csc::ExprCsc;
use crate::expr::normalize::Normalization;
use crate::pipeline::stage2_normalize::ExprMatrix;
use std::collections::HashMap;
use std::fs;
use tempfile::tempdir;

fn build_gene_index() -> GeneIndex {
    let mut idx = GeneIndex {
        rows: Vec::new(),
        duplicates: Vec::new(),
        first_index_by_symbol: HashMap::new(),
    };
    idx.first_index_by_symbol.insert("A".to_string(), 1);
    idx.first_index_by_symbol.insert("B".to_string(), 2);
    idx.first_index_by_symbol.insert("C".to_string(), 3);
    idx
}

#[test]
fn panel_accumulation_correctness() {
    let dir = tempdir().expect("tempdir");
    let mtx = dir.path().join("matrix.mtx");
    fs::write(
        &mtx,
        "%%MatrixMarket matrix coordinate integer general\n3 2 3\n1 1 1\n2 1 2\n3 2 3\n",
    )
    .expect("write file");

    let (expr, stats) = ExprCsc::from_mtx(&mtx, 3, 2, false).expect("csc");
    let expr_ctx = ExprContext {
        expr: ExprMatrix::Owned(expr),
        cell_stats: stats,
        normalization: Normalization {
            enabled: false,
            scale: 10_000.0,
            epsilon: 1e-8,
        },
    };

    let panels = PanelSet {
        panels: vec![crate::panels::defs::PanelDef {
            id: "P1".to_string(),
            description: "".to_string(),
            axis: "X".to_string(),
            genes: vec![
                crate::panels::defs::PanelGene {
                    symbol: "A".to_string(),
                },
                crate::panels::defs::PanelGene {
                    symbol: "B".to_string(),
                },
                crate::panels::defs::PanelGene {
                    symbol: "C".to_string(),
                },
            ],
            required: vec!["A".to_string()],
            weights: None,
        }],
    };

    let cell_ids = vec!["c1".to_string(), "c2".to_string()];

    let out_dir = dir.path().join("out");
    fs::create_dir_all(&out_dir).expect("mkdir");

    let ctx = run_stage3_panels(&expr_ctx, &panels, &build_gene_index(), &cell_ids, &out_dir)
        .expect("stage3");
    assert_eq!(ctx.mappings.len(), 1);

    let report = fs::read_to_string(out_dir.join("panels_report.tsv")).expect("report");
    assert!(report.contains("c1\tP1\tX\t3.000000\t2\t1.000000\t0"));
    assert!(report.contains("c2\tP1\tX\t3.000000\t1\t1.000000\t0"));
}

#[test]
fn determinism_report_bytes() {
    let dir = tempdir().expect("tempdir");
    let mtx = dir.path().join("matrix.mtx");
    fs::write(
        &mtx,
        "%%MatrixMarket matrix coordinate integer general\n2 1 2\n1 1 1\n2 1 2\n",
    )
    .expect("write file");

    let (expr, stats) = ExprCsc::from_mtx(&mtx, 2, 1, false).expect("csc");
    let expr_ctx = ExprContext {
        expr: ExprMatrix::Owned(expr),
        cell_stats: stats,
        normalization: Normalization::default(),
    };
    let panels = PanelSet {
        panels: vec![crate::panels::defs::PanelDef {
            id: "P1".to_string(),
            description: "".to_string(),
            axis: "X".to_string(),
            genes: vec![
                crate::panels::defs::PanelGene {
                    symbol: "A".to_string(),
                },
                crate::panels::defs::PanelGene {
                    symbol: "B".to_string(),
                },
            ],
            required: vec!["A".to_string()],
            weights: None,
        }],
    };
    let mut idx = GeneIndex {
        rows: Vec::new(),
        duplicates: Vec::new(),
        first_index_by_symbol: HashMap::new(),
    };
    idx.first_index_by_symbol.insert("A".to_string(), 1);
    idx.first_index_by_symbol.insert("B".to_string(), 2);

    let cell_ids = vec!["c1".to_string()];

    let out1 = dir.path().join("out1");
    let out2 = dir.path().join("out2");
    fs::create_dir_all(&out1).expect("mkdir");
    fs::create_dir_all(&out2).expect("mkdir");

    run_stage3_panels(&expr_ctx, &panels, &idx, &cell_ids, &out1).expect("stage3-1");
    run_stage3_panels(&expr_ctx, &panels, &idx, &cell_ids, &out2).expect("stage3-2");

    let bytes1 = fs::read(out1.join("panels_report.tsv")).expect("read1");
    let bytes2 = fs::read(out2.join("panels_report.tsv")).expect("read2");
    assert_eq!(bytes1, bytes2);
}
