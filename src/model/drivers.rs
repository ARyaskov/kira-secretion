#[derive(Debug, Clone)]
pub struct PanelDriver {
    pub panel_id: String,
    pub score: f32,
}

pub fn top_k_panels(panel_ids: &[String], contributions: &[f32], k: usize) -> Vec<PanelDriver> {
    let mut pairs: Vec<PanelDriver> = panel_ids
        .iter()
        .zip(contributions.iter())
        .map(|(id, v)| PanelDriver {
            panel_id: id.clone(),
            score: *v,
        })
        .collect();

    pairs.sort_by(|a, b| {
        match b
            .score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
        {
            std::cmp::Ordering::Equal => a.panel_id.cmp(&b.panel_id),
            other => other,
        }
    });

    pairs.truncate(k);
    pairs
}

pub fn top_k_eeb_drivers(
    export_ids: &[String],
    export_vals: &[f32],
    degrade_ids: &[String],
    degrade_vals: &[f32],
    k: usize,
) -> (Vec<PanelDriver>, Vec<PanelDriver>) {
    let export = top_k_panels(export_ids, export_vals, k);
    let mut degrade = top_k_panels(degrade_ids, degrade_vals, k);
    for d in &mut degrade {
        d.score = -d.score.abs();
    }
    (export, degrade)
}

pub fn format_drivers(drivers: &[PanelDriver]) -> String {
    if drivers.is_empty() {
        return ".".to_string();
    }
    let mut parts = Vec::with_capacity(drivers.len());
    for d in drivers {
        parts.push(format!("{}={:.4}", d.panel_id, d.score));
    }
    parts.join(",")
}

pub fn format_eeb_drivers(export: &[PanelDriver], degrade: &[PanelDriver]) -> String {
    let export_str = format_drivers(export);
    let degrade_str = format_drivers(degrade);
    format!("EXPORT:{};DEGRADE:{}", export_str, degrade_str)
}

pub fn top_k_components(names: &[&str], contribs: &[f32], k: usize) -> String {
    if names.is_empty() || contribs.is_empty() {
        return ".".to_string();
    }
    let mut pairs: Vec<(String, f32)> = names
        .iter()
        .zip(contribs.iter())
        .map(|(n, v)| ((*n).to_string(), *v))
        .collect();
    pairs.sort_by(
        |a, b| match b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal) {
            std::cmp::Ordering::Equal => a.0.cmp(&b.0),
            other => other,
        },
    );
    if pairs.len() > k {
        pairs.truncate(k);
    }
    let mut out = Vec::with_capacity(pairs.len());
    for (name, value) in pairs {
        out.push(format!("{}={:.4}", name, value));
    }
    out.join(",")
}

#[cfg(test)]
#[path = "../../tests/src_inline/model/drivers.rs"]
mod tests;
