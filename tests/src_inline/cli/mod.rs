
use super::*;
use clap::Parser;

#[test]
fn run_mode_default_is_standalone() {
    let cli = Cli::parse_from(["kira-secretion", "run", "--input", "in", "--out", "out"]);
    match cli.command {
        Command::Run(args) => {
            assert_eq!(args.run_mode, run::RunModeArg::Standalone);
        }
        _ => panic!("expected run command"),
    }
}

#[test]
fn run_mode_pipeline_parses() {
    let cli = Cli::parse_from([
        "kira-secretion",
        "run",
        "--input",
        "in",
        "--out",
        "out",
        "--run-mode",
        "pipeline",
    ]);
    match cli.command {
        Command::Run(args) => {
            assert_eq!(args.run_mode, run::RunModeArg::Pipeline);
        }
        _ => panic!("expected run command"),
    }
}
