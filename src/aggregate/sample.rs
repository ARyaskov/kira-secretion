use std::collections::BTreeMap;

use crate::model::flags::Flags;
use crate::model::regimes::Regime;

pub fn median_f32(values: &mut [f32]) -> f32 {
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    percentile(values, 0.5)
}

pub fn median_ignore_nan(values: &mut Vec<f32>) -> f32 {
    values.retain(|v| !v.is_nan());
    if values.is_empty() {
        return f32::NAN;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    percentile(values, 0.5)
}

fn percentile(values: &[f32], p: f32) -> f32 {
    if values.is_empty() {
        return f32::NAN;
    }
    let n = values.len();
    let idx = ((p * (n as f32 - 1.0)).floor() as usize).min(n - 1);
    values[idx]
}

pub fn majority_regime(regimes: &[Regime]) -> Regime {
    let mut counts: BTreeMap<Regime, usize> = BTreeMap::new();
    for r in regimes {
        *counts.entry(*r).or_insert(0) += 1;
    }
    let mut best = Regime::Unclassified;
    let mut best_count = 0usize;
    for r in Regime::ordered() {
        let c = *counts.get(r).unwrap_or(&0);
        if c > best_count {
            best_count = c;
            best = *r;
        }
    }
    best
}

pub fn flags_from_fraction(flags: &[Flags], threshold: f32) -> Flags {
    let mut out = Flags::empty();
    let n = flags.len() as f32;
    if n == 0.0 {
        return out;
    }
    let bits = [
        Flags::LOW_CONFIDENCE,
        Flags::FEW_DETECTED_GENES,
        Flags::LOW_COUNTS,
        Flags::HIGH_AMBIENT_RISK,
    ];
    for bit in bits {
        let c = flags.iter().filter(|f| f.contains(bit)).count() as f32;
        if c / n >= threshold {
            out.set(bit);
        }
    }
    out
}
