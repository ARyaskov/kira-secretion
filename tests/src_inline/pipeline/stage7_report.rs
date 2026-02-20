use super::*;
use crate::expr::csc::{CellStats, ExprCsc};
use crate::expr::normalize::Normalization;
use crate::input::detect::TenXFormat;
use crate::input::features::GeneIndex;
use crate::model::axes::{AxisCoverage, AxisValues};
use crate::model::regimes::RuleId;
use crate::panels::defs::{PanelDef, PanelGene, PanelSet};
use crate::panels::mapping::GeneMapping;
use crate::pipeline::stage2_normalize::ExprMatrix;
use crate::pipeline::stage3_panels::{PanelCellPacked, PanelsContext};
use crate::pipeline::stage4_axes::{
    AxesContext, AxesSummary, AxisDrivers, AxisStats, AxisSummaryEntry,
};
use crate::pipeline::stage5_scores::{CompositeStats, CompositesSummary, ScoresContext};
use crate::pipeline::stage6_classify::{ClassifyContext, RegimeSummary as Stage6RegimeSummary};
use tempfile::tempdir;

fn dummy_dataset() -> DatasetCtx {
    DatasetCtx {
        format: TenXFormat::TenXv3,
        matrix_path: "matrix.mtx".into(),
        features_path: "features.tsv".into(),
        barcodes_path: "barcodes.tsv".into(),
        shared_cache_path: None,
        resolved_shared_cache_path: None,
        gene_index: GeneIndex {
            rows: vec![],
            duplicates: vec![],
            first_index_by_symbol: HashMap::new(),
        },
        barcodes: vec!["c1".to_string(), "c2".to_string()],
        n_genes: 2,
        n_cells: 2,
        nnz: 2,
        duplicate_gene_symbols_count: 0,
        duplicate_gene_symbols: vec![],
        meta_present: false,
        meta_cells_matched: 0,
        meta_cells_missing: 0,
    }
}

fn dummy_expr() -> ExprContext {
    ExprContext {
        expr: ExprMatrix::Owned(ExprCsc {
            n_genes: 2,
            n_cells: 2,
            nnz: 2,
            col_ptr: vec![0, 1, 2],
            row_idx: vec![0, 1],
            values: vec![10, 20],
        }),
        cell_stats: vec![
            CellStats {
                libsize: 1000,
                detected: 10,
            },
            CellStats {
                libsize: 2000,
                detected: 20,
            },
        ],
        normalization: Normalization::default(),
    }
}

fn dummy_axes() -> AxesContext {
    AxesContext {
        cell_ids: vec!["c1".to_string(), "c2".to_string()],
        values: vec![
            AxisValues {
                sia: 0.6,
                eeb: 0.2,
                sli: 0.5,
                mei: 0.4,
                ecmi: 0.4,
                apci: 0.3,
                gdi: 0.2,
            },
            AxisValues {
                sia: 0.2,
                eeb: -0.2,
                sli: 0.1,
                mei: 0.1,
                ecmi: 0.1,
                apci: 0.1,
                gdi: 0.8,
            },
        ],
        coverage: vec![
            AxisCoverage {
                sia: 0.9,
                eeb: 0.9,
                sli: 0.9,
                mei: 0.9,
                ecmi: 0.9,
                apci: 0.9,
                gdi: 0.9,
            },
            AxisCoverage {
                sia: 0.5,
                eeb: 0.5,
                sli: 0.5,
                mei: 0.5,
                ecmi: 0.5,
                apci: 0.5,
                gdi: 0.5,
            },
        ],
        drivers: vec![
            AxisDrivers {
                sia: "".to_string(),
                eeb: "".to_string(),
                sli: "".to_string(),
                mei: "".to_string(),
                ecmi: "".to_string(),
                apci: "".to_string(),
                gdi: "".to_string(),
            },
            AxisDrivers {
                sia: "".to_string(),
                eeb: "".to_string(),
                sli: "".to_string(),
                mei: "".to_string(),
                ecmi: "".to_string(),
                apci: "".to_string(),
                gdi: "".to_string(),
            },
        ],
        stats: AxesSummary {
            sia: zero_axis_summary(),
            eeb: zero_axis_summary(),
            sli: zero_axis_summary(),
            mei: zero_axis_summary(),
            ecmi: zero_axis_summary(),
            apci: zero_axis_summary(),
            gdi: zero_axis_summary(),
        },
    }
}

fn zero_axis_summary() -> AxisSummaryEntry {
    AxisSummaryEntry {
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
    }
}

fn dummy_scores() -> ScoresContext {
    ScoresContext {
        oii: vec![0.7, 0.1],
        iai: vec![0.6, 0.2],
        esi: vec![0.65, 0.15],
        cov_oii: vec![0.9, 0.5],
        cov_iai: vec![0.9, 0.5],
        cov_esi: vec![0.9, 0.5],
        drivers_oii: vec!["".to_string(), "".to_string()],
        drivers_iai: vec!["".to_string(), "".to_string()],
        drivers_esi: vec!["".to_string(), "".to_string()],
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

fn dummy_classify() -> ClassifyContext {
    ClassifyContext {
        regimes: vec![Regime::EnvironmentShaping, Regime::SelfPreserving],
        rule_ids: vec![RuleId::R7EnvironmentShaping, RuleId::R1SelfPreserving],
        flags: vec![Flags::empty(), Flags::empty()],
        summary: Stage6RegimeSummary {
            counts: vec![],
            fractions: vec![],
            flagged_fractions: vec![],
        },
    }
}

fn dummy_panels() -> PanelsContext {
    PanelsContext {
        panels: PanelSet {
            panels: vec![PanelDef {
                id: "P1".to_string(),
                description: "Panel One".to_string(),
                axis: "SIA".to_string(),
                genes: vec![PanelGene {
                    symbol: "G1".to_string(),
                }],
                required: vec!["G1".to_string()],
                weights: None,
            }],
        },
        mappings: vec![GeneMapping {
            panel_id: "P1".to_string(),
            mapped: vec![Some(0)],
            required_hits: 1,
            required_total: 1,
        }],
        warnings: vec![],
        cell_ids: vec!["c1".to_string(), "c2".to_string()],
        per_cell: vec![
            PanelCellPacked {
                sums: vec![1.0],
                hits: vec![1],
                required_missing: vec![0],
            },
            PanelCellPacked {
                sums: vec![2.0],
                hits: vec![1],
                required_missing: vec![0],
            },
        ],
    }
}

#[test]
fn secretion_tsv_header_exact_order() {
    let dir = tempdir().expect("tempdir");
    run_stage7_report(
        &dummy_dataset(),
        &dummy_expr(),
        &dummy_axes(),
        &dummy_scores(),
        &dummy_classify(),
        &dummy_panels(),
        dir.path(),
        "cell",
        RunMode::Standalone,
        None,
    )
    .expect("stage7");

    let txt = std::fs::read_to_string(dir.path().join("secretion.tsv")).expect("read");
    let header = txt.lines().next().unwrap_or("");
    assert_eq!(
        header,
        "barcode\tsample\tcondition\tspecies\tlibsize\tnnz\texpressed_genes\tsecretory_load\texocytosis_bias\tvesicle_traffic_intensity\ter_golgi_pressure\tparacrine_signal_potential\tstress_secretion_index\tregime\tflags\tconfidence"
    );
}

#[test]
fn summary_json_schema() {
    let dir = tempdir().expect("tempdir");
    run_stage7_report(
        &dummy_dataset(),
        &dummy_expr(),
        &dummy_axes(),
        &dummy_scores(),
        &dummy_classify(),
        &dummy_panels(),
        dir.path(),
        "cell",
        RunMode::Standalone,
        None,
    )
    .expect("stage7");

    let v: serde_json::Value =
        serde_json::from_slice(&std::fs::read(dir.path().join("summary.json")).expect("read"))
            .expect("json");
    assert!(v.get("tool").is_some());
    assert!(v.get("input").is_some());
    assert!(v.get("distributions").is_some());
    assert!(v.get("regimes").is_some());
    assert!(v.get("qc").is_some());
    assert!(v["distributions"]["secretory_load"]["median"].is_number());
}

#[test]
fn pipeline_step_json_schema() {
    let dir = tempdir().expect("tempdir");
    run_stage7_report(
        &dummy_dataset(),
        &dummy_expr(),
        &dummy_axes(),
        &dummy_scores(),
        &dummy_classify(),
        &dummy_panels(),
        dir.path(),
        "cell",
        RunMode::Pipeline,
        None,
    )
    .expect("stage7");

    let v: serde_json::Value = serde_json::from_slice(
        &std::fs::read(dir.path().join("pipeline_step.json")).expect("read"),
    )
    .expect("json");
    assert_eq!(v["tool"]["name"], "kira-secretion");
    assert_eq!(v["tool"]["stage"], "secretion");
    assert_eq!(v["artifacts"]["primary_metrics"], "secretion.tsv");
    assert_eq!(v["cell_metrics"]["file"], "secretion.tsv");
    assert_eq!(v["cell_metrics"]["id_column"], "barcode");
    assert_eq!(v["cell_metrics"]["regime_column"], "regime");
    assert_eq!(v["cell_metrics"]["confidence_column"], "confidence");
    assert_eq!(v["cell_metrics"]["flag_column"], "flags");
    assert!(v["regimes"].is_array());
}

#[test]
fn deterministic_outputs() {
    let dir = tempdir().expect("tempdir");

    run_stage7_report(
        &dummy_dataset(),
        &dummy_expr(),
        &dummy_axes(),
        &dummy_scores(),
        &dummy_classify(),
        &dummy_panels(),
        dir.path(),
        "cell",
        RunMode::Pipeline,
        None,
    )
    .expect("stage7-1");

    let s1 = std::fs::read(dir.path().join("secretion.tsv")).expect("read1");
    let p1 = std::fs::read(dir.path().join("panels_report.tsv")).expect("read1");
    let j1 = std::fs::read(dir.path().join("summary.json")).expect("read1");
    let m1 = std::fs::read(dir.path().join("pipeline_step.json")).expect("read1");

    run_stage7_report(
        &dummy_dataset(),
        &dummy_expr(),
        &dummy_axes(),
        &dummy_scores(),
        &dummy_classify(),
        &dummy_panels(),
        dir.path(),
        "cell",
        RunMode::Pipeline,
        None,
    )
    .expect("stage7-2");

    let s2 = std::fs::read(dir.path().join("secretion.tsv")).expect("read2");
    let p2 = std::fs::read(dir.path().join("panels_report.tsv")).expect("read2");
    let j2 = std::fs::read(dir.path().join("summary.json")).expect("read2");
    let m2 = std::fs::read(dir.path().join("pipeline_step.json")).expect("read2");

    assert_eq!(s1, s2);
    assert_eq!(p1, p2);
    assert_eq!(j1, j2);
    assert_eq!(m1, m2);
}

#[test]
fn pipeline_step_written_only_in_pipeline_mode() {
    let dir = tempdir().expect("tempdir");
    run_stage7_report(
        &dummy_dataset(),
        &dummy_expr(),
        &dummy_axes(),
        &dummy_scores(),
        &dummy_classify(),
        &dummy_panels(),
        dir.path(),
        "cell",
        RunMode::Standalone,
        None,
    )
    .expect("stage7");
    assert!(!dir.path().join("pipeline_step.json").exists());

    run_stage7_report(
        &dummy_dataset(),
        &dummy_expr(),
        &dummy_axes(),
        &dummy_scores(),
        &dummy_classify(),
        &dummy_panels(),
        dir.path(),
        "cell",
        RunMode::Pipeline,
        None,
    )
    .expect("stage7");
    assert!(dir.path().join("pipeline_step.json").exists());
}
