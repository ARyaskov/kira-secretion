use clap::Parser;
use kira_secretion::cli::Cli;
use kira_secretion::simd;
use tracing_subscriber::EnvFilter;
use tracing_subscriber::fmt::time::UtcTime;

fn main() -> anyhow::Result<()> {
    let filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info"));
    tracing_subscriber::fmt()
        .with_env_filter(filter)
        .with_timer(UtcTime::rfc_3339())
        .with_target(false)
        .init();

    tracing::info!(
        simd_backend = simd::backend_name(),
        "simd backend selected at build time"
    );

    let cli = Cli::parse();
    cli.dispatch()
}
