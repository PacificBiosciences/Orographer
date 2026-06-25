# User Guide

## Installation

Orographer requires Python 3.10+ and can be installed from source or 
[Bioconda](https://bioconda.github.io/). We recommend a dedicated 
conda/mamba environment:

```bash
mamba create -n orographer_env
mamba activate orographer_env
mamba install -c bioconda orographer
```

## How to run

### Command-line arguments

After installation, run `orographer plot -h` and `orographer deploy -h` to see
options.

```text
Orographer v0.3.0
usage: orographer [-h] [-v] COMMAND ...

Parse BAM file and genomic coordinates, then log the information and analyze split alignments.

positional arguments:
  COMMAND        Command to run
    plot         Generate alignment plot
    deploy       Start HTTP server to serve generated plots

options:
  -h, --help     show this help message and exit
  -v, --version  Installed version (0.3.0)
```

#### Basic usage
**Plot command** (`orographer plot -h`):

- `--bam`: Path to the input BAM file (required)
- `--coord`: Genomic coordinate(s) in format `chrom:start-end` (e.g.
  `chr1:1000-2000`). Can be given multiple times for multiple regions.
- `--region-type`: Region type to plot: `paraphase` (required).
- `--ref`: Path to reference FASTA file (required for mismatch visualization)
- `--outdir`: Directory for output HTML and JSON files (required)
- `--prefix`: Alphanumeric prefix (with underscores) for output filenames
- `--gtf`: _Optional_ path to bgzip compressed GTF/GFF3 annotation file with
  `.tbi` index for gene tracks. Index with: `tabix -p gff file.gtf.gz`
- `--vcf`: _Optional_ path to VCF file for variant track. If a `.tbi` index
  exists, tabix is used for region access.
- `--sample-label`: Display label for the primary BAM (--bam).
- `--verbose`: Verbose logging to stderr

#### Additional options for multiple samples:
- `--other-bam`: Path to an additional BAM file. May be specified up to two
  times. Order in plot: first other-bam at top, second below it, primary --bam
  at bottom.
- `--other-vcf`: Path to VCF for the corresponding --other-bam. May be
  specified up to two times.
- `--other-sample-label`: Display label for --other-bam. Specify once per
  --other-bam. May be specified up to two times.

Multiple samples are supported for all region types.

**Deploy command** (`orographer deploy -h`):

- `--outdir`: Directory containing the generated HTML and JSON files to serve (required)
- `--port`: _Optional_ port to serve on (default: 8000)

### Basic usage by region type

#### Paraphase region type

Paraphase renders requested region(s) directly with phased read, coverage,
gene, and VCF context.

```bash
orographer plot \
  --bam input.bam \
  --coord chr1:1000-2000 \
  --region-type paraphase \
  --ref reference.fa \
  --outdir ./output
```

With optional gene and variant tracks:

```bash
orographer plot \
  --bam input.bam \
  --coord chr1:1000-2000 \
  --region-type paraphase \
  --ref reference.fa \
  --gtf genes.gtf.gz \
  --vcf variants.vcf.gz \
  --outdir ./output
```

Multiple regions can be supplied with repeated `--coord` arguments:

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

The plot command writes HTML and JSON under `--outdir`. To view locally, serve
the directory and open in a browser:

```bash
orographer deploy --outdir ./output --port 8000
```

Then open `http://localhost:8000/` and navigate to the generated HTML file(s).
File names include the region (and optional prefix), e.g.
`chr1_1000_2000_bokeh.html`.

Plots must be served over HTTP(S); opening the HTML file directly via `file://`
will not load the external JSON data.

## Reference dotplot

Each plot includes a **reference self-identity dotplot** — a square dotplot 
that compares the reference sequence against itself to reveal repetitive regions.
Click the "Show ref identity" button in the coordinate toolbar to open it. See 
[Reference Self-Identity Dotplot](dotplot.md) for details.

### Detailed visualization guides

- [Paraphase Visualization Guide](paraphase.md)
