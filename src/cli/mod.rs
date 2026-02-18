use clap::{Parser, Subcommand};

mod panels;
mod run;
mod validate;

#[derive(Parser, Debug)]
#[command(name = "kira-secretion", version, about = "Kira Secretion CLI")]
pub struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    Run(run::RunArgs),
    Validate(validate::ValidateArgs),
    Panels(panels::PanelsArgs),
}

impl Cli {
    pub fn dispatch(self) -> anyhow::Result<()> {
        match self.command {
            Command::Run(args) => run::handle(args),
            Command::Validate(args) => validate::handle(args),
            Command::Panels(args) => panels::handle(args),
        }
    }
}

#[cfg(test)]
#[path = "../../tests/src_inline/cli/mod.rs"]
mod tests;
