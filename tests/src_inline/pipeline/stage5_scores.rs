
use super::*;
use crate::model::axes::{AxisCoverage, AxisValues};
use crate::pipeline::stage4_axes::{
    AxesContext, AxesSummary, AxisDrivers, AxisStats, AxisSummaryEntry,
};
use tempfile::tempdir;

fn dummy_axes(values: AxisValues, coverage: AxisCoverage) -> AxesContext {
    AxesContext {
        cell_ids: vec!["c1".to_string()],
        values: vec![values],
        coverage: vec![coverage],
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

#[test]
fn composite_correctness_apci_present() {
    let axes = dummy_axes(
        AxisValues {
            sia: 0.5,
            eeb: 0.0,
            sli: 0.2,
            mei: 0.4,
            ecmi: 0.3,
            apci: 0.6,
            gdi: 0.1,
        },
        AxisCoverage {
            sia: 1.0,
            eeb: 1.0,
            sli: 1.0,
            mei: 1.0,
            ecmi: 1.0,
            apci: 1.0,
            gdi: 1.0,
        },
    );
    let dir = tempdir().expect("tempdir");
    let scores = run_stage5_scores(&axes, dir.path()).expect("scores");
    let eeb_pos = 0.5;
    let expected =
        clamp01(0.22 * 0.5 + 0.18 * eeb_pos + 0.12 * 0.2 + 0.16 * 0.4 + 0.16 * 0.3 + 0.16 * 0.1);
    assert!((scores.oii[0] - expected).abs() < 1e-6);
}

#[test]
fn composite_correctness_apci_absent() {
    let axes = dummy_axes(
        AxisValues {
            sia: 0.2,
            eeb: -0.2,
            sli: 0.1,
            mei: 0.4,
            ecmi: 0.3,
            apci: f32::NAN,
            gdi: 0.5,
        },
        AxisCoverage {
            sia: 1.0,
            eeb: 1.0,
            sli: 1.0,
            mei: 1.0,
            ecmi: 1.0,
            apci: 0.0,
            gdi: 1.0,
        },
    );
    let dir = tempdir().expect("tempdir");
    let scores = run_stage5_scores(&axes, dir.path()).expect("scores");
    let eeb_pos = pos_eeb(-0.2);
    let expected = clamp01(0.30 * 0.4 + 0.30 * 0.5 + 0.25 * 0.2 + 0.15 * eeb_pos);
    assert!((scores.iai[0] - expected).abs() < 1e-6);
}

#[test]
fn coverage_weighted_mean() {
    let axes = dummy_axes(
        AxisValues {
            sia: 0.0,
            eeb: 0.0,
            sli: 0.0,
            mei: 0.0,
            ecmi: 0.0,
            apci: 0.0,
            gdi: 0.0,
        },
        AxisCoverage {
            sia: 0.5,
            eeb: 1.0,
            sli: 0.0,
            mei: 1.0,
            ecmi: 0.5,
            apci: 0.0,
            gdi: 1.0,
        },
    );
    let dir = tempdir().expect("tempdir");
    let scores = run_stage5_scores(&axes, dir.path()).expect("scores");
    let w = WeightsDefault::default();
    let expected = weighted_cov_oii(&axes.coverage[0], &w);
    assert!((scores.cov_oii[0] - expected).abs() < 1e-6);
}

#[test]
fn determinism_composites_tsv() {
    let axes = dummy_axes(
        AxisValues {
            sia: 0.1,
            eeb: 0.2,
            sli: 0.3,
            mei: 0.4,
            ecmi: 0.5,
            apci: 0.6,
            gdi: 0.7,
        },
        AxisCoverage {
            sia: 1.0,
            eeb: 1.0,
            sli: 1.0,
            mei: 1.0,
            ecmi: 1.0,
            apci: 1.0,
            gdi: 1.0,
        },
    );
    let dir = tempdir().expect("tempdir");
    let out1 = dir.path().join("out1");
    let out2 = dir.path().join("out2");
    std::fs::create_dir_all(&out1).expect("mkdir");
    std::fs::create_dir_all(&out2).expect("mkdir");
    run_stage5_scores(&axes, &out1).expect("scores1");
    run_stage5_scores(&axes, &out2).expect("scores2");
    let a = std::fs::read(out1.join("composites.tsv")).expect("read1");
    let b = std::fs::read(out2.join("composites.tsv")).expect("read2");
    assert_eq!(a, b);
}
