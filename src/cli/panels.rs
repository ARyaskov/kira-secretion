use std::path::PathBuf;

use clap::{Args, Subcommand};

use crate::panels::loader::{default_panels_dir, load_panels_from_dir};

#[derive(Args, Debug)]
pub struct PanelsArgs {
    #[command(subcommand)]
    command: PanelsCommand,
}

#[derive(Subcommand, Debug)]
enum PanelsCommand {
    List,
    Dump(PanelsDumpArgs),
}

#[derive(Args, Debug)]
pub struct PanelsDumpArgs {
    /// Output directory
    #[arg(long)]
    out: PathBuf,
}

pub fn handle(args: PanelsArgs) -> anyhow::Result<()> {
    match args.command {
        PanelsCommand::List => list_panels(),
        PanelsCommand::Dump(args) => dump_panels(args),
    }
}

fn list_panels() -> anyhow::Result<()> {
    let dir = default_panels_dir();
    let panels = load_panels_from_dir(&dir)?;
    println!("panel_id\taxis\tn_genes\tn_required");
    for panel in panels.panels {
        println!(
            "{}\t{}\t{}\t{}",
            panel.id,
            panel.axis,
            panel.genes.len(),
            panel.required.len()
        );
    }
    Ok(())
}

fn dump_panels(args: PanelsDumpArgs) -> anyhow::Result<()> {
    std::fs::create_dir_all(&args.out)?;
    let dir = default_panels_dir();
    let panels = load_panels_from_dir(&dir)?;
    let json = serde_json::to_string_pretty(&panels)?;
    let path = args.out.join("panels_manifest.json");
    std::fs::write(path, json)?;
    Ok(())
}
