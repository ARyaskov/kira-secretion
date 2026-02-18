
use super::*;
use crate::model::axes::{AxisCoverage, AxisValues};
use crate::pipeline::stage2_normalize::ExprMatrix;
use crate::pipeline::stage4_axes::{
    AxesContext, AxesSummary, AxisDrivers, AxisStats, AxisSummaryEntry,
};
use crate::pipeline::stage5_scores::{CompositeStats, CompositesSummary, ScoresContext};
use std::collections::HashMap;
use tempfile::tempdir;

fn dummy_axes(values: AxisValues) -> AxesContext {
    AxesContext {
        cell_ids: vec!["c1".to_string()],
        values: vec![values],
        coverage: vec![AxisCoverage {
            sia: 1.0,
            eeb: 1.0,
            sli: 1.0,
            mei: 1.0,
            ecmi: 1.0,
            apci: 1.0,
            gdi: 1.0,
        }],
        drivers: vec![AxisDrivers {
            sia: "".to_string(),
            eeb: "".to_string(),
            sli: "".to_string(),
            mei: "".to_string(),
            ecmi: "".to_string(),
            apci: "".to_string(),
            gdi: "".to_string(),
        }],
        stats: AxesSummary {
            sia: AxisSummaryEntry {
                present: true,
                value: AxisStats {
                    median: 0.0,
                    p90: 0.0,
                    p99: 0.0,
                    frac_ge_0_65: 0.0,
                    frac_ge_0_80: 0.0,
                },
                coverage: AxisStats {
                    median: 0.0,
                    p90: 0.0,
                    p99: 0.0,
                    frac_ge_0_65: 0.0,
                    frac_ge_0_80: 0.0,
                },
            },
            eeb: AxisSummaryEntry {
                present: true,
                value: AxisStats {
                    median: 0.0,
                    p90: 0.0,
                    p99: 0.0,
                    frac_ge_0_65: 0.0,
                    frac_ge_0_80: 0.0,
                },
                coverage: AxisStats {
                    median: 0.0,
                    p90: 0.0,
                    p99: 0.0,
                    frac_ge_0_65: 0.0,
                    frac_ge_0_80: 0.0,
                },
            },
            sli: AxisSummaryEntry {
                present: true,
                value: AxisStats {
                    median: 0.0,
                    p90: 0.0,
                    p99: 0.0,
                    frac_ge_0_65: 0.0,
                    frac_ge_0_80: 0.0,
                },
                coverage: AxisStats {
                    median: 0.0,
                    p90: 0.0,
                    p99: 0.0,
                    frac_ge_0_65: 0.0,
                    frac_ge_0_80: 0.0,
                },
            },
            mei: AxisSummaryEntry {
                present: true,
                value: AxisStats {
                    median: 0.0,
                    p90: 0.0,
                    p99: 0.0,
                    frac_ge_0_65: 0.0,
                    frac_ge_0_80: 0.0,
                },
                coverage: AxisStats {
                    median: 0.0,
                    p90: 0.0,
                    p99: 0.0,
                    frac_ge_0_65: 0.0,
                    frac_ge_0_80: 0.0,
                },
            },
            ecmi: AxisSummaryEntry {
                present: true,
                value: AxisStats {
                    median: 0.0,
                    p90: 0.0,
                    p99: 0.0,
                    frac_ge_0_65: 0.0,
                    frac_ge_0_80: 0.0,
                },
                coverage: AxisStats {
                    median: 0.0,
                    p90: 0.0,
                    p99: 0.0,
                    frac_ge_0_65: 0.0,
                    frac_ge_0_80: 0.0,
                },
            },
            apci: AxisSummaryEntry {
                present: true,
                value: AxisStats {
                    median: 0.0,
                    p90: 0.0,
                    p99: 0.0,
                    frac_ge_0_65: 0.0,
                    frac_ge_0_80: 0.0,
                },
                coverage: AxisStats {
                    median: 0.0,
                    p90: 0.0,
                    p99: 0.0,
                    frac_ge_0_65: 0.0,
                    frac_ge_0_80: 0.0,
                },
            },
            gdi: AxisSummaryEntry {
                present: true,
                value: AxisStats {
                    median: 0.0,
                    p90: 0.0,
                    p99: 0.0,
                    frac_ge_0_65: 0.0,
                    frac_ge_0_80: 0.0,
                },
                coverage: AxisStats {
                    median: 0.0,
                    p90: 0.0,
                    p99: 0.0,
                    frac_ge_0_65: 0.0,
                    frac_ge_0_80: 0.0,
                },
            },
        },
    }
}

fn dummy_scores(oii: f32, esi: f32) -> ScoresContext {
    ScoresContext {
        oii: vec![oii],
        iai: vec![0.0],
        esi: vec![esi],
        cov_oii: vec![1.0],
        cov_iai: vec![1.0],
        cov_esi: vec![1.0],
        drivers_oii: vec!["".to_string()],
        drivers_iai: vec!["".to_string()],
        drivers_esi: vec!["".to_string()],
        summary: CompositesSummary {
            oii: CompositeStats {
                median: 0.0,
                p90: 0.0,
                p99: 0.0,
                frac_ge_0_65: 0.0,
                frac_ge_0_80: 0.0,
            },
            iai: CompositeStats {
                median: 0.0,
                p90: 0.0,
                p99: 0.0,
                frac_ge_0_65: 0.0,
                frac_ge_0_80: 0.0,
            },
            esi: CompositeStats {
                median: 0.0,
                p90: 0.0,
                p99: 0.0,
                frac_ge_0_65: 0.0,
                frac_ge_0_80: 0.0,
            },
        },
    }
}

fn dummy_dataset(n: usize) -> DatasetCtx {
    let dir = tempdir().expect("tempdir");
    let barcodes = dir.path().join("barcodes.tsv");
    std::fs::write(&barcodes, "c1\n").expect("write");
    DatasetCtx {
        format: crate::input::detect::TenXFormat::TenXv3,
        matrix_path: dir.path().join("matrix.mtx"),
        features_path: dir.path().join("features.tsv"),
        barcodes_path: barcodes,
        shared_cache_path: None,
        resolved_shared_cache_path: None,
        gene_index: crate::input::features::GeneIndex {
            rows: Vec::new(),
            duplicates: Vec::new(),
            first_index_by_symbol: HashMap::new(),
        },
        barcodes: (0..n).map(|i| format!("c{}", i + 1)).collect(),
        n_genes: 1,
        n_cells: n,
        nnz: 0,
        duplicate_gene_symbols_count: 0,
        duplicate_gene_symbols: Vec::new(),
        meta_present: false,
        meta_cells_matched: 0,
        meta_cells_missing: 0,
    }
}

#[test]
fn rule_boundary_self_preserving() {
    let t = Thresholds::default();
    let axis = AxisValues {
        sia: 0.2,
        eeb: -0.2,
        sli: 0.1,
        mei: 0.2,
        ecmi: 0.2,
        apci: 0.0,
        gdi: 0.2,
    };
    let (regime, rule) = classify_cell(&axis, pos_eeb(axis.eeb), 0.0, 0.0, &t);
    assert_eq!(regime, Regime::SelfPreserving);
    assert_eq!(rule, RuleId::R1SelfPreserving);
}

#[test]
fn rule_boundary_secretory_lysosome_active() {
    let t = Thresholds::default();
    let axis = AxisValues {
        sia: 0.5,
        eeb: 0.0,
        sli: 0.8,
        mei: 0.2,
        ecmi: 0.2,
        apci: 0.0,
        gdi: 0.2,
    };
    let (regime, rule) = classify_cell(&axis, pos_eeb(axis.eeb), 0.0, 0.0, &t);
    assert_eq!(regime, Regime::SecretoryLysosomeActive);
    assert_eq!(rule, RuleId::R2SecretoryLysosomeActive);
}

#[test]
fn rule_boundary_export_dominant() {
    let t = Thresholds::default();
    let axis = AxisValues {
        sia: 0.6,
        eeb: 0.6,
        sli: 0.1,
        mei: 0.2,
        ecmi: 0.2,
        apci: 0.0,
        gdi: 0.2,
    };
    let (regime, rule) = classify_cell(&axis, pos_eeb(axis.eeb), 0.61, 0.0, &t);
    assert_eq!(regime, Regime::ExportDominant);
    assert_eq!(rule, RuleId::R3ExportDominant);
}

#[test]
fn rule_boundary_metabolic_suppressive() {
    let t = Thresholds::default();
    let axis = AxisValues {
        sia: 0.4,
        eeb: 0.1,
        sli: 0.1,
        mei: 0.8,
        ecmi: 0.2,
        apci: 0.0,
        gdi: 0.2,
    };
    let (regime, rule) = classify_cell(&axis, pos_eeb(axis.eeb), 0.0, 0.0, &t);
    assert_eq!(regime, Regime::MetabolicSuppressive);
    assert_eq!(rule, RuleId::R4MetabolicSuppressive);
}

#[test]
fn rule_boundary_inflammatory_signaler() {
    let t = Thresholds::default();
    let axis = AxisValues {
        sia: 0.4,
        eeb: 0.0,
        sli: 0.1,
        mei: 0.2,
        ecmi: 0.2,
        apci: 0.0,
        gdi: 0.8,
    };
    let (regime, rule) = classify_cell(&axis, pos_eeb(axis.eeb), 0.0, 0.0, &t);
    assert_eq!(regime, Regime::InflammatorySignaler);
    assert_eq!(rule, RuleId::R5InflammatorySignaler);
}

#[test]
fn rule_boundary_presentation_high() {
    let t = Thresholds::default();
    let axis = AxisValues {
        sia: 0.5,
        eeb: 0.0,
        sli: 0.1,
        mei: 0.2,
        ecmi: 0.2,
        apci: 0.8,
        gdi: 0.2,
    };
    let (regime, rule) = classify_cell(&axis, pos_eeb(axis.eeb), 0.0, 0.0, &t);
    assert_eq!(regime, Regime::PresentationHigh);
    assert_eq!(rule, RuleId::R6PresentationHigh);
}

#[test]
fn rule_boundary_environment_shaping() {
    let t = Thresholds::default();
    let axis = AxisValues {
        sia: 0.6,
        eeb: 0.0,
        sli: 0.1,
        mei: 0.6,
        ecmi: 0.6,
        apci: 0.0,
        gdi: 0.2,
    };
    let (regime, rule) = classify_cell(&axis, pos_eeb(axis.eeb), 0.70, 0.70, &t);
    assert_eq!(regime, Regime::EnvironmentShaping);
    assert_eq!(rule, RuleId::R7EnvironmentShaping);
}

#[test]
fn flags_low_counts_and_detected() {
    let axes = dummy_axes(AxisValues {
        sia: 0.5,
        eeb: 0.0,
        sli: 0.1,
        mei: 0.1,
        ecmi: 0.1,
        apci: 0.0,
        gdi: 0.1,
    });
    let scores = dummy_scores(0.0, 0.0);
    let dataset = dummy_dataset(1);
    let expr = ExprContext {
        expr: ExprMatrix::Owned(crate::expr::csc::ExprCsc {
            n_genes: 0,
            n_cells: 1,
            nnz: 0,
            col_ptr: vec![0, 0],
            row_idx: vec![],
            values: vec![],
        }),
        cell_stats: vec![crate::expr::csc::CellStats {
            libsize: 100,
            detected: 10,
        }],
        normalization: crate::expr::normalize::Normalization::default(),
    };
    let dir = tempdir().expect("tempdir");
    let ctx = run_stage6_classify(&dataset, &expr, &axes, &scores, dir.path()).expect("classify");
    let f = ctx.flags[0];
    assert!(f.contains(Flags::LOW_COUNTS));
    assert!(f.contains(Flags::FEW_DETECTED_GENES));
}

#[test]
fn flags_low_confidence_and_ambient() {
    let mut axes = dummy_axes(AxisValues {
        sia: 0.4,
        eeb: 0.0,
        sli: 0.1,
        mei: 0.1,
        ecmi: 0.1,
        apci: 0.0,
        gdi: 0.8,
    });
    axes.coverage[0].sia = 0.5;
    let scores = dummy_scores(0.0, 0.0);
    let dataset = dummy_dataset(1);
    let expr = ExprContext {
        expr: ExprMatrix::Owned(crate::expr::csc::ExprCsc {
            n_genes: 0,
            n_cells: 1,
            nnz: 0,
            col_ptr: vec![0, 0],
            row_idx: vec![],
            values: vec![],
        }),
        cell_stats: vec![crate::expr::csc::CellStats {
            libsize: 1000,
            detected: 10,
        }],
        normalization: crate::expr::normalize::Normalization::default(),
    };
    let dir = tempdir().expect("tempdir");
    let ctx = run_stage6_classify(&dataset, &expr, &axes, &scores, dir.path()).expect("classify");
    let f = ctx.flags[0];
    assert!(f.contains(Flags::LOW_CONFIDENCE));
    assert!(f.contains(Flags::FEW_DETECTED_GENES));
    assert!(f.contains(Flags::HIGH_AMBIENT_RISK));
}

#[test]
fn determinism_classify_tsv() {
    let axes = dummy_axes(AxisValues {
        sia: 0.5,
        eeb: 0.0,
        sli: 0.1,
        mei: 0.1,
        ecmi: 0.1,
        apci: 0.0,
        gdi: 0.1,
    });
    let scores = dummy_scores(0.0, 0.0);
    let dataset = dummy_dataset(1);
    let expr = ExprContext {
        expr: ExprMatrix::Owned(crate::expr::csc::ExprCsc {
            n_genes: 0,
            n_cells: 1,
            nnz: 0,
            col_ptr: vec![0, 0],
            row_idx: vec![],
            values: vec![],
        }),
        cell_stats: vec![crate::expr::csc::CellStats {
            libsize: 1000,
            detected: 1000,
        }],
        normalization: crate::expr::normalize::Normalization::default(),
    };
    let dir = tempdir().expect("tempdir");
    let out1 = dir.path().join("out1");
    let out2 = dir.path().join("out2");
    std::fs::create_dir_all(&out1).expect("mkdir");
    std::fs::create_dir_all(&out2).expect("mkdir");
    run_stage6_classify(&dataset, &expr, &axes, &scores, &out1).expect("c1");
    run_stage6_classify(&dataset, &expr, &axes, &scores, &out2).expect("c2");
    let a = std::fs::read(out1.join("classify.tsv")).expect("read1");
    let b = std::fs::read(out2.join("classify.tsv")).expect("read2");
    assert_eq!(a, b);
}

#[test]
fn summary_counts_and_fractions() {
    let regimes = vec![
        Regime::SelfPreserving,
        Regime::Unclassified,
        Regime::Unclassified,
    ];
    let flags = vec![Flags::empty(), Flags::empty(), Flags::empty()];
    let summary = summarize(&regimes, &flags);
    let count_self = summary
        .counts
        .iter()
        .find(|(r, _)| *r == Regime::SelfPreserving)
        .unwrap()
        .1;
    let count_un = summary
        .counts
        .iter()
        .find(|(r, _)| *r == Regime::Unclassified)
        .unwrap()
        .1;
    assert_eq!(count_self, 1);
    assert_eq!(count_un, 2);
    let frac_un = summary
        .fractions
        .iter()
        .find(|(r, _)| *r == Regime::Unclassified)
        .unwrap()
        .1;
    assert!((frac_un - 2.0 / 3.0).abs() < 1e-6);
}
