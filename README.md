# kira-secretion

Deterministic, explainable secretion-state QC for single-cell expression data.

## Build requirements

- Rust >= 1.95

## Install

Install from crates.io:

```bash
cargo install kira-secretion
```

## Usage examples

Standalone run (cell mode):

```bash
kira-secretion run \
  --input ./data/inf \
  --out ./out/inf \
  --mode cell \
  --run-mode standalone
```

Standalone run (sample mode):

```bash
kira-secretion run \
  --input ./data/inf \
  --out ./out/inf-sample \
  --mode sample \
  --run-mode standalone
```

Pipeline run (shared cache lookup + pipeline artifacts):

```bash
kira-secretion run \
  --input ./data/scc \
  --out ./out/scc-cell \
  --mode cell \
  --run-mode pipeline
```

Validation command:

```bash
kira-secretion validate \
  --input ./data/inf \
  --run-mode pipeline
```

Panels listing:

```bash
kira-secretion panels list
```

Panels manifest dump:

```bash
kira-secretion panels dump --out ./out/panels
```

## Modes

- `--run-mode standalone` (default): standard MTX/TSV input flow.
- `--run-mode pipeline`: pipeline contract mode for `kira-organelle`.

## Pipeline cache lookup

In pipeline mode, `kira-secretion` first searches for shared cache in the input directory:

- exact expected file by prefix:
  - no prefix: `kira-organelle.bin`
  - prefixed dataset: `<PREFIX>.kira-organelle.bin`
- if exact file is missing: any file ending with `kira-organelle.bin` (deterministic lexicographic choice when multiple match)

Behavior:

- cache exists and valid: use shared cache path.
- cache missing: warn once and fall back to MTX input.
- cache exists but invalid: hard error (no silent fallback).

## Pipeline output contract

In pipeline mode, outputs are written to:

- `--out <DIR>` -> `<DIR>/kira-secretion/`

Required artifacts:

- `secretion.tsv` (per-cell contract table)
- `summary.json` (run-level aggregates)
- `panels_report.tsv` (panel audit)
- `pipeline_step.json` (ingestion manifest for `kira-organelle`)

All TSV float values are fixed `%.6f`.

## Shared cache specification

- Canonical format: `CACHE_FILE.md`
