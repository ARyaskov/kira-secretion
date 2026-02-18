
use super::*;

#[test]
fn drivers_tie_break() {
    let ids = vec!["B".to_string(), "A".to_string(), "C".to_string()];
    let vals = vec![1.0, 1.0, 0.5];
    let drivers = top_k_panels(&ids, &vals, 2);
    assert_eq!(drivers[0].panel_id, "A");
    assert_eq!(drivers[1].panel_id, "B");
}

#[test]
fn components_tie_break() {
    let names = vec!["B", "A", "C"];
    let vals = vec![0.5, 0.5, 0.4];
    let out = top_k_components(&names, &vals, 2);
    assert_eq!(out, "A=0.5000,B=0.5000");
}
