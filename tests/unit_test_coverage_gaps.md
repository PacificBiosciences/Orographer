# Unit Test Coverage Gaps Report

This file summarizes gaps identified in the current unit/integration test coverage with a goal of **production-quality** confidence (minimal permissiveness, meaningful assertions, and deterministic behavior).

## 1) What tests should exist but don’t?

1. **Output filename generation is untested at the unit level**
   - Functions like `generate_output_filename()` and `generate_multi_region_filename()` are used heavily but are not directly unit-tested for:
     - prefix behavior
     - sanitization rules (especially multi-region colon/hyphen replacement)
     - exact expected output naming for single vs multi-region.
   - Status: **COMPLETED** (added `tests/unit/test_plot_filenames.py` with exact filename assertions).

2. **BAM parser correctness beyond smoke/invariants**
   - Current tests focus on structural invariants and ordering, but don’t validate meaningful correctness such as:
     - `extract_variants()` producing expected SNV/INS/DEL calls for known reads and a known reference slice
     - SA parsing deterministically influencing split segment content (e.g., verifying specific SA-derived positions appear)
     - `order_split_alignments()` tie-breaking and robustness for unusual SA patterns (missing/odd entries).
   - Status: **COMPLETED (partial)** — `test_extract_variants.py` for variant tuples; **`tests/unit/test_order_split_alignments.py`** for SA ordering, malformed-SA fallback to query-span sort, and fallback tie-breaks (`order_split_alignments` skips bad SA fields; empty valid keys → interval sort).

3. **`only_split=True` is not deterministically tested**
   - Unit tests allow that `only_split=True` may return empty results, without asserting a known fixture region where SA-tagged reads must exist.
   - Production confidence benefits from at least one deterministic “SA-present” scenario asserting non-empty results and expected properties derived from SA.
   - Status: **COMPLETED** — (a) synthetic indexed BAM + `SA` tag: `test_fetch_all_alignments_only_split_sa_present`; (b) **fixture** path: chr16 hba + minimum **55** SA-tagged reads in `test_fetch_all_alignments_only_split_smoke` (see §2 “Recently completed” item 1).

4. **VCF tabix/indexed branch**
   - `parse_vcf_file()` includes a `.tbi`/tabix code path, but tests only exercise the uncompressed branch.
   - Production confidence should cover at least one indexed VCF scenario to validate:
     - contig naming matching
     - tabix coordinate transforms and filtering behavior.
   - Status: **COMPLETED** (`test_vcf_parser_tabix.py`: tabix branch, **`chr1`↔`1` aliasing** via `get_tabix_chromosome_name` + `parse_vcf_file` on both contig styles)

5. **GTF/GFF3 branch coverage**
   - Current GTF tests cover basic parsing and overlap, but not:
     - the GFF3 parsing attribute branch
     - edge cases in `detect_annotation_format()`
     - attribute fallback logic (e.g., gene name fallback behaviors).
   - Status: **COMPLETED** (added `tests/unit/test_gtf_gff3_branch.py` + `write_bgzip_tabix_gff3` helper: GFF3 tabix parse, `detect_annotation_format`, `parse_gff3_attributes`, `extract_gff3_gene_info`, name fallback when `Name` absent).

6. **CLI workflow edge cases**
   - `tests/unit/test_cli_validation.py` covers some validation errors, but missing coverage for:
     - multi-`--coord` behavior at the unit level (it’s only indirectly checked in integration)
     - deploy failure modes (wrong `--outdir`, port issues) outside of slow smoke
     - interactions like prefix + multi-region naming consistency as a direct assertion.
   - Status: **COMPLETED** (added `tests/unit/test_cli_workflow_edge_cases.py`: multi-`--coord` plot success, `--prefix` + two coords output basename, `deploy` with missing `--outdir` exits non-zero, second `--coord` invalid fails before writing HTML).

7. **Plot JSON schema deeper validation**
   - Integration tests assert top-level JSON keys (`docs_json`, `render_items`) and coordinate token presence.
   - Production-quality confidence would benefit from schema-derived assertions on specific stable internal JSON fields/structure (not “anywhere in the payload” substring matching).
   - Status: **COMPLETED** (`tests/helpers/plot_json_schema.py`: `validate_plot_json_schema`, genomic x-range extraction; `tests/unit/test_plot_json_schema.py`; plot/complex_sv/multisample integration tests assert schema + expected `(start, end)` in `docs_json`).

## 2) What tests are more permissive than they should be? *(updated after gap fixes)*

Many earlier “permissive” gaps are **narrowed** by new tests (`extract_variants`, deterministic `only_split` via synthetic SA BAM, `plot_json_schema`, CLI workflow cases). Below is what is **still** looser than ideal for production-grade rigor, plus a few **new** caveats introduced by the current approach.

### Recently completed (former “still permissive” items 1–3)

1. **Fixture-based `only_split` smoke** — **COMPLETED**
   - **Improvement:** `test_fetch_all_alignments_only_split_smoke` now uses the same **chr16 hba** window as the paraphase non-split test and **requires** at least `_MIN_READS_CHR16_PARAPHASE_ONLY_SPLIT` (55) SA-tagged reads, so an empty split result fails the suite. This complements the synthetic SA BAM test with **real** SA/CIGAR data on a fixed fixture region.

2. **BAM segment geometry vs query region** — **COMPLETED**
   - **Improvement:** `_assert_segment_overlaps_query_region()` checks that every segment whose `chrom` matches the queried chromosome has `[pos, end)` overlapping the BED-derived 1-based region (converted to 0-based half-open). Applied to fixture paraphase tests and to segments in the synthetic `only_split` SA test. Segments on other chromosomes (e.g. SA off-target) are skipped.

3. **Minimum counts on fixture BAM smoke paths** — **COMPLETED**
   - **Improvement:** For **HG002.paraphase.small** @ **chr16:171801–175500**, non-split mode must return at least **60** read keys and at least **60** total segment lists; `only_split` must return at least **55** each. Constants live in `tests/unit/test_bam_parser.py` so bounds can be raised if the fixture is refreshed.

### Recently completed (former §2 “still permissive” items 1–3 — CLI / SA / VCF alias)

4. **CLI exit codes vs stage** — **COMPLETED**
   - **Improvement:** `parser.error()` paths assert **`returncode == 2`**; coordinate `ValueError` and deploy missing dir assert **`== 1`**. (`test_cli_validation.py`, `test_cli_workflow_edge_cases.py`.)

5. **`order_split_alignments` / odd SA** — **COMPLETED**
   - **Improvement:** `bam_parser.order_split_alignments` skips malformed SA entries; if none parse, falls back to query-span ordering. **`tests/unit/test_order_split_alignments.py`:** SA order vs tag, malformed skip, full-malformed fallback, tie-breaks (longer span, primary vs supplementary).

6. **VCF tabix `chr` vs numeric contig** — **COMPLETED**
   - **Improvement:** `test_vcf_parser_tabix.py` asserts `get_tabix_chromosome_name` maps **chr1→1** on numeric-CHROM VCF and **1→chr1** on chr-prefixed VCF; `parse_vcf_file` returns expected positions for both query styles.

### Recently completed (former §2 “still permissive” plot JSON + sample labels)

7. **Plot JSON genomic / region bounds** — **COMPLETED**
   - **Improvement:** `collect_orig_region_bounds_from_docs_json` + `assert_orig_region_bounds_in_docs_json` read **CustomJS** `orig_start` / `orig_end` args (orographer’s stable 1-based bounds). Integration tests assert these **before** heuristic x-range checks. `collect_genomic_like_start_end_pairs(..., always_include_ranges=...)` lets expected regions below the usual magnitude floor still match when explicitly requested.

8. **Sample labels in JSON** — **COMPLETED**
   - **Improvement:** `assert_sample_labels_in_docs_json` requires each label on a **Bokeh Div** (`attributes.text` exact match) **and** in **ColumnDataSource** `sample_label` (plain dict or Bokeh 3 **`type: map` / `entries`** encoding). Multisample integration uses this instead of raw JSON substring checks.

### Still permissive (residual / optional)

- **Heuristic x-ranges:** `assert_genomic_ranges_in_docs_json` still supplements orig-bound checks using start/end pairs (with expected ranges always included); unrelated tiny numeric pairs elsewhere remain filtered by the default floor.

### Resolved (no longer “overly permissive” in the way originally described)

- **Plot coordinate checks:** integration tests now use **`validate_plot_json_schema` + `assert_genomic_ranges_in_docs_json`** so region bounds must appear as plausible x-ranges inside `docs_json`, not arbitrary substring hits in the full payload.
- **Deterministic split alignment path:** synthetic SA BAM test **plus** fixture chr16 hba minimum-count `only_split` test (§2 items 1 & 3).
- **BAM fixture strictness:** segment overlap vs query region + minimum read counts on chr16 hba (§2 items 2 & 3).
- **CLI / SA / VCF:** exit code staging, `order_split_alignments` robustness + unit tests, tabix contig aliasing (§2 items 4–6).
- **Plot JSON / labels:** orig_start/orig_end assertions, always-include genomic pairs, structured sample-label checks (§2 items 7–8).

