# Unit Tests Overview

These are the current pytest suites under `tests/`. Most tests are fast “real code path” checks with minimal/no mocking.

## Unit tests (`tests/unit/*`)

`tests/unit/test_utils.py`
- `parse_coordinate()` happy-path + validation of invalid formats/values.
- `validate_output_dir()` creates temp dirs and rejects non-directories.

`tests/unit/test_bam_parser.py`
- `validate_bam_file()` raises on missing BAMs.
- Fixture **HG002** @ **chr16 hba**: minimum read counts for paraphase non-split / `only_split`; segments on the query chrom must **overlap** the region; synthetic indexed BAM with `SA` tag for deterministic split parsing.

`tests/unit/test_vcf_parser.py`
- Uses a tiny synthetic uncompressed VCF to validate:
  - SNP/INS/DEL classification
  - region overlap filtering
  - genotype-to-haplotype extraction (`haplotypes`)

`tests/unit/test_gtf_parser.py`
- Uses a tiny synthetic GTF written as `bgzip + tabix` indexed `.gtf.gz` to validate:
  - `parse_annotation_file()` returns overlapping gene annotations
  - missing `.tbi` triggers a `ValueError`

`tests/unit/test_cli_validation.py`
- CLI rejects invalid `--prefix`, too many `--other-bam` / `--other-vcf`, mismatched `--other-sample-label`, and bad `--coord` format.

`tests/unit/test_cli_workflow_edge_cases.py`
- Multi-`--coord` plot succeeds; `--prefix` + two regions yields expected multi-region filenames; `deploy` with missing `--outdir` fails; invalid second `--coord` fails without writing plot HTML.

## Integration tests (`tests/integration/*`)

`tests/unit/test_plot_json_schema.py`
- Validates plot `*_bokeh.json.gz` structure (`docs_json` + `render_items`) and that genomic x-ranges appear inside `docs_json`.

`tests/integration/test_cli_plot_smoke.py`
- Runs `orographer plot` via the CLI and asserts that at least one `*_bokeh.html` and `*_bokeh.json.gz` are produced in `--outdir`.
- Asserts deeper JSON schema and genomic range bounds in `docs_json`.

`tests/integration/test_deploy_smoke.py`
- Starts `orographer deploy` and fetches an expected HTML page over HTTP.
- Marked `@pytest.mark.slow`.

## How to run

- Fast: `mamba run -n paraviewer -- pytest -m "not slow" -q`
- Slow (deploy): `mamba run -n paraviewer -- pytest -m slow -q`

