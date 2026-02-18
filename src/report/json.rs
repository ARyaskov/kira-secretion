use crate::pipeline::stage7_report::FinalSummary;

pub type Summary = FinalSummary;

pub fn write_summary(out_dir: &std::path::Path, summary: &Summary) -> anyhow::Result<()> {
    let json = serde_json::to_string_pretty(summary)?;
    let path = out_dir.join("summary.json");
    std::fs::write(path, json)?;
    Ok(())
}
