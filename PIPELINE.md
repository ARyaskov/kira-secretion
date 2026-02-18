# Pipeline

`kira-secretion run` executes a fixed 7-stage pipeline:

1. `stage1_load`
2. `stage2_normalize`
3. `stage3_panels`
4. `stage4_axes`
5. `stage5_scores`
6. `stage6_classify`
7. `stage7_report`

## Run modes

- `--run-mode standalone` (default): outputs go directly to `--out`.
- `--run-mode pipeline`: outputs go to `--out/kira-secretion/` and include pipeline manifest for downstream ingestion.

## Stage-by-stage outputs

1. `stage1_load`
- Discovers input source (shared cache vs MTX/TSV), validates dimensions/metadata, builds `DatasetCtx`.
- No direct artifact file.

2. `stage2_normalize`
- Loads expression matrix (from shared cache or MTX) and computes per-cell stats.
- No direct artifact file.

3. `stage3_panels`
- Computes per-cell panel accumulations and mapping coverage.
- Writes `panels_report.tsv` (per-cell panel diagnostics; intermediate).

4. `stage4_axes`
- Builds secretion axes + coverage + axis drivers.
- Writes `axes.tsv`.

5. `stage5_scores`
- Computes composite scores (OII/IAI/ESI), coverage, and score drivers.
- Writes `composites.tsv`.

6. `stage6_classify`
- Assigns regime/rule/flags from axes + composites + QC thresholds.
- Writes `classify.tsv`.

7. `stage7_report`
- Produces final contract-facing tables and aggregates.
- Writes:
  - `secretion.tsv` (primary per-cell contract table; barcode-sorted)
  - `summary.json` (deterministic aggregated summary)
  - `panels_report.tsv` (final panel-level aggregate report; replaces stage3 intermediate file)
  - `report.txt`
  - `pipeline_step.json` (only in `--run-mode pipeline`)

## Shared cache resolution (pipeline mode)

In `--run-mode pipeline`, Stage 1 resolves shared cache in this order:

1. If `--cache <PATH>` is provided: use it directly.
2. Try expected cache name by detected prefix:
  - no prefix: `kira-organelle.bin`
  - prefixed dataset: `<PREFIX>.kira-organelle.bin`
3. If exact expected name is absent: pick any file ending with `kira-organelle.bin` (deterministic lexicographic choice if multiple).
4. If no cache found: warn once and fallback to MTX/TSV (`matrix.mtx[.gz]`, `features.tsv/genes.tsv[.gz]`, `barcodes.tsv[.gz]`).

Failure policy:

- cache exists and valid: use shared cache path
- cache missing: fallback to MTX/TSV
- cache exists but invalid: hard error (no silent fallback)

Validation rules follow `CACHE_FILE.md` (header/magic/version/endian/header-size/file-bytes, CRC64-ECMA, section bounds, string tables, CSC invariants).

## Pipeline aggregator contract (`kira-organelle` ingestion)

In `--run-mode pipeline`, downstream `kira-organelle` must read only:

- `summary.json`
- `secretion.tsv`
- `pipeline_step.json`

`pipeline_step.json` defines contract-relevant artifact mapping:

- `artifacts.summary = "summary.json"`
- `artifacts.primary_metrics = "secretion.tsv"`
- `artifacts.panels = "panels_report.tsv"`
- `cell_metrics.file = "secretion.tsv"`
- `cell_metrics.id_column = "barcode"`
- `cell_metrics.regime_column = "regime"`
- `cell_metrics.confidence_column = "confidence"`
- `cell_metrics.flag_column = "flags"`
