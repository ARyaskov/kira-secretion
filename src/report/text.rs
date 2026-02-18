use crate::pipeline::stage7_report::FinalSummary;

pub fn render_report(summary: &FinalSummary) -> String {
    let mut out = String::new();
    out.push_str("Kira Secretion Report\n");
    out.push_str("======================\n\n");
    out.push_str("This report summarizes transcript-derived proxy signals. ");
    out.push_str("It does not measure proteins, does not establish causality, and should be interpreted conservatively.\n\n");

    out.push_str("Dataset overview:\n");
    out.push_str(&format!("- Cells: {}\n", summary.input.n_cells));
    out.push_str(&format!("- Species: {}\n\n", summary.input.species));

    out.push_str("Dominant regimes:\n");
    let top = top_regimes(&summary.regimes.fractions, 2);
    for (name, frac) in top {
        out.push_str(&format!("- {}: {:.2}%\n", name, frac * 100.0));
    }
    out.push_str("\n");

    out.push_str("Distribution tails:\n");
    out.push_str(&format!(
        "- Secretory load p99: {:.4}\n",
        summary.distributions.secretory_load.p99
    ));
    out.push_str(&format!(
        "- ER-Golgi pressure p99: {:.4}\n",
        summary.distributions.er_golgi_pressure.p99
    ));
    out.push_str(&format!(
        "- Stress secretion index p99: {:.4}\n",
        summary.distributions.stress_secretion_index.p99
    ));
    out.push_str("\n");

    out.push_str("Confidence and QC flags:\n");
    out.push_str(&format!(
        "- LOW_CONFIDENCE: {:.2}%\n",
        summary.qc.low_confidence_fraction * 100.0
    ));
    out.push_str(&format!(
        "- LOW_SECRETORY_SIGNAL: {:.2}%\n",
        summary.qc.low_secretory_signal_fraction * 100.0
    ));
    out.push_str("\n");

    out
}

fn top_regimes(regimes: &std::collections::BTreeMap<String, f32>, k: usize) -> Vec<(String, f32)> {
    let mut pairs: Vec<(String, f32)> = regimes.iter().map(|(r, f)| (r.clone(), *f)).collect();
    pairs.sort_by(
        |a, b| match b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal) {
            std::cmp::Ordering::Equal => a.0.cmp(&b.0),
            other => other,
        },
    );
    pairs.truncate(k);
    pairs
}
