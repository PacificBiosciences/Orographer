# User Guide

## Installation

Orographer requires Python 3.10+ and can be installed from source. We recommend a dedicated conda/mamba environment:

```bash
mamba create -n orographer_env "python>=3.10" pip
mamba activate orographer_env
cd orographer
pip install .
```

## How to run

### Command-line arguments

After installation, run `orographer plot -h` and `orographer deploy -h` to see options:

```text
usage: orographer [-h] [-v] {plot,deploy} ...

Commands:
  plot    Generate alignment plot
  deploy  Start HTTP server to serve generated plots
```

#### Basic usage
**Plot command** (`orographer plot -h`):

- `--bam`: Path to the input BAM file (required)
- `--coord`: Genomic coordinate(s) in format `chrom:start-end` (e.g. `chr1:1000-2000`). Can be given multiple times for multiple regions, tested for up to 3.
- `--region-type`: Currently prefer `paraphase`, additional option `complex_sv` coming soon (required)
- `--ref`: Path to reference FASTA file (required for mismatch visualization)
- `--outdir`: Directory for output HTML and JSON files (required)
- `--prefix`: Alphanumeric prefix (with underscores) for output filenames
- `--gtf`: _Optional_ path to bgzip compressed GTF/GFF3 annotation file with `.tbi` index for gene track (optional). Index with: `tabix -p gff file.gtf.gz`
- `--vcf`: _Optional_ path to VCF file for variant track (optional). If a `.tbi` index exists, tabix is used for region access.
- `--sample-label`: Display label for the primary BAM (--bam).
- `--verbose`: Verbose logging to stderr

#### Additional options for multiple samples:
- `--other-bam`: Path to an additional BAM file. May be specified up to two times. Order in plot: first other-bam at top, second below it, primary --bam at bottom.
- `--other-vcf`: Path to VCF for the corresponding --other-bam (order matches: first --other-vcf for first --other-bam). May be specified up to two times.
- `--other-sample-label`: Display label for --other-bam. Specify once per --other-bam (order matches). May be specified up to two times.

**Deploy command** (`orographer deploy -h`):

- `--outdir`: Directory containing the generated HTML and JSON files to serve (required)
- `--port`: _Optional_ port to serve on (default: 8000)

### Basic usage

Single region:

```bash
orographer plot \
  --bam input.bam \
  --coord chr1:1000-2000 \
  --region-type complex_sv \
  --ref reference.fa \
  --outdir ./output
```

With optional gene and variant tracks:

```bash
orographer plot \
  --bam input.bam \
  --coord chr1:1000-2000 \
  --region-type complex_sv \
  --ref reference.fa \
  --gtf genes.gtf.gz \
  --vcf variants.vcf.gz \
  --outdir ./output
```

Multiple regions:

```bash
orographer plot \
  --bam input.bam \
  --coord chr1:1000-2000 \
  --coord chr1:5000-6000 \
  --region-type paraphase \
  --ref reference.fa \
  --outdir ./output
```

### Coordinate format

Coordinates must be `chromosome:start-end`:

- **chromosome**: Any string (e.g. `chr1`, `1`, `chrX`)
- **start**, **end**: Integer positions (1-based, inclusive)
- Start must be less than end

Examples: `chr1:1000-2000`, `1:150000-160000`, `chrX:50000-60000`.

### Results and viewing

The plot command writes HTML and JSON under `--outdir`. To view locally, serve the directory and open in a browser:

```bash
orographer deploy --outdir ./output --port 8000
```

Then open `http://localhost:8000/` and navigate to the generated HTML file(s). File names include the region (and optional prefix), e.g. `chr1_1000_2000_bokeh.html`.

Plots must be served over HTTP(S); opening the HTML file directly via `file://` will not load the external JSON data.