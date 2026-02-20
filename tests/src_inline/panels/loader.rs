use super::*;

#[test]
fn load_core_panels() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("assets/panels");
    let set = load_panels_from_dir(&dir).expect("load panels");
    assert!(set.panels.len() >= 4);
    assert_eq!(set.panels[0].id, "ER_GOLGI_TRAFFICKING");
    assert_eq!(set.panels[0].genes[0].symbol, "SEC23A");
}
