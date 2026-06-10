# Paraphase Visualization Guide

## Table Of Contents

- [Region Context](#region-context)
  - [Input requirements](#input-requirements)
  - [Reference dotplot](#reference-dotplot)
  - [Gene track](#gene-track)
  - [Multi-BAM visualization](#multi-bam-visualization)
- [Read Evidence](#read-evidence)
  - [Which reads are included](#which-reads-are-included)
  - [Read representation](#read-representation)
  - [Alignment numbering](#alignment-numbering)
  - [Read names and haplotypes](#read-names-and-haplotypes)
- [Haplotype And Phase](#haplotype-and-phase)
  - [Haplotype blocks](#haplotype-blocks)
  - [Phase set delimiters](#phase-set-delimiters)
  - [Coverage by haplotype](#coverage-by-haplotype)
- [Coverage And Variants](#coverage-and-variants)
  - [Coverage plots](#coverage-plots)
  - [Read-level variant annotations](#read-level-variant-annotations)
  - [VCF track](#vcf-track)
- [Interaction And Navigation](#interaction-and-navigation)
  - [Navigation features](#navigation-features)
  - [Read selection options](#read-selection-options)
  - [Display controls](#display-controls)

Paraphase plots are designed for phased, region-focused alignment review. The
display shows reads from the requested coordinate interval, organizes them by
haplotype when HP tags are available, overlays read-level sequence differences,
and links coverage, variants, genes, and navigation controls to a shared genomic
x-axis.

## Region Context

### Input requirements

Paraphase mode is selected with `--region-type paraphase`. It renders the
explicit coordinate intervals supplied on the command line and does not perform
connected-region discovery. Coordinates use the same `chrom:start-end` format
as other Orographer modes, including comma-formatted positions such as
`chr5:70,917,001-70,961,200`.

A BAM file and reference FASTA are required. The reference FASTA is used for
read-level mismatch, insertion, and deletion visualization. When a reference is
not available to the parser, the plot can still render alignments, but
reference-base comparisons are not available for read-level variant markers.

Paraphase BAMs commonly include HP, PS, and YC auxiliary tags. HP tags identify
haplotype assignment, PS tags identify phase-set starts, and YC tags provide
per-read color values from upstream analysis. Orographer reads these tags when
present but still displays alignments without them in an unassigned or neutral
style.

### Reference dotplot

The clickable reference dotplot thumbnail in the global top bar is a
self-similarity plot generated from the reference sequence for all displayed
regions. When multiple Paraphase coordinates are supplied, the reference
sequences are concatenated in display order for the dotplot only. The regions
are not biologically adjacent; boundary markers in the dotplot show where one
region ends and the next begins.

Each axis represents the same concatenated reference blocks. The sequence is
split into fixed bins, and each bin is represented by canonical 9-mers sampled
every 10 bp. Canonical k-mers treat a k-mer and its reverse complement as the
same signal, so both direct repeats and inverted repeats can produce
off-diagonal matches.

Dense diagonal signal is expected because each block matches itself.
Off-diagonal signal within a block can indicate repeats, duplicated sequence,
inversions, or other local structures that can make read placement ambiguous.
Dotplots are computed only when reference sequence is available and the total
concatenated span is within the configured 10 Mb dotplot size limit. Dotplot
generation can also be disabled with `--no-dotplot`.

### Gene track

The gene track appears when a bgzip-indexed GTF or GFF3 annotation file is
provided. It shows genes and exons overlapping the displayed region. Clicking gene
features opens details for the selected annotation, including gene name,
feature type, coordinates, strand, and source annotation fields when available.
Orographer selects one transcript per gene to display. Representative transcript
selection uses the first available criterion in this order:

1. MANE Select
2. RefSeq Select
3. Ensembl Canonical
4. APPRIS principal
5. GENCODE Primary
6. Most exons

If multiple transcripts satisfy the same criterion, Orographer prefers the
transcript with the most exons, then applies stable coordinate and transcript
ID tie-breakers.

The gene track is shared across samples for a region. In multi-BAM views, the
primary BAM's parsed annotation is used for the region-level gene track because
the annotation is coordinate-based rather than sample-specific.

### Multi-BAM visualization

When multiple BAMs are provided, each sample is rendered as a separate stacked
row for the same displayed region. Additional BAMs supplied with `--other-bam`
are shown above the primary `--bam`, and the primary BAM is shown at the
bottom. Sample labels appear above each sample's tracks when labels are
provided.

Each sample row has its own coverage, VCF, read, haplotype, and phase-set
evidence, but the rows share the same genomic x-axis. Zooming, panning, and
coordinate jumps keep the rows aligned so evidence can be compared across
samples at the same locus. Multiple requested coordinates are displayed as
separate region panels rather than one discovered multi-locus event.

## Read Evidence

### Which reads are included

Paraphase mode fetches alignments overlapping the requested region from each
BAM. Secondary alignments are skipped. Primary and supplementary alignments are
retained when they pass the shared mapping-quality threshold, currently MAPQ 5.
Reads with SA tags are expanded into their split or supplementary alignment
segments so all visible pieces of the read can be reviewed together.

When a read has multiple records at the same reference start, Orographer
deduplicates that read/start combination so the same alignment is not drawn
twice. Supplementary records can update segment HP and PS tags, aligned blocks,
and read-level variant annotations when they provide more complete information
than the primary record's SA-derived representation.

Unlike complex SV mode, Paraphase mode is not restricted to reads with complex
or split evidence. It displays all qualifying reads in the requested interval,
making it useful for whole-region haplotype review and comparison against
coverage and variant tracks.

### Read representation

Each read is drawn as one or more arrow-like alignment segments on the genomic
x-axis. The arrow direction indicates read orientation. A single read can have
multiple segments when supplementary or SA-tagged alignments are present. Each
segment is clipped to the visible coordinate range so reads extending beyond
the requested interval remain readable inside the panel.

If a YC color tag is present and contains a valid RGB value, that color is used
for the segment. Otherwise, Paraphase segments use the neutral Paraphase read
color rather than strand-specific red or blue. This preserves upstream
Paraphase coloring when available while keeping uncolored reads visually
consistent.

Mismatches, insertions, and deletions are drawn on top of the alignment when
reference and CIGAR evidence are available. At wider zoom levels, these
annotations are compact markers so dense regions remain readable. When the
visible range is below 1 kb, mismatch and insertion markers switch to text
labels, such as the alternate base for a mismatch or an insertion count label.

### Alignment numbering

Small numbered labels identify individual alignment segments for a read. The
numbers are assigned in read-coordinate order using forward-strand read
positions, so they provide a stable visual handle for split or supplementary
alignments. Clicking a number opens an alignment detail modal with the read
name, alignment number, coordinates, strand, haplotype, sample label, and all
displayed alignment coordinates for that read.

Alignment numbers are especially useful when a read has several placements or
when two samples show different segmentation for similar read evidence. Labels
can be hidden with the "Hide algn numbers" checkbox when the plot is dense.
When a read is selected, labels for unselected reads are faded or hidden so the
selected read is easier to follow.

### Read names and haplotypes

Reads are grouped by display name and haplotype. When a segment has a usable HP
tag, Orographer appends an `_HP` suffix to the display read key, such as
`readA_HP1`. This allows the same raw read name to be shown separately if
different segments carry different haplotype assignments. Reads without HP tags
remain under the raw read name and are treated as unassigned.

The alignment detail modal reports haplotype as `HP:n` when a haplotype tag is
available and "Unassigned" otherwise. This makes tag state explicit in the
popup rather than requiring the user to infer it only from row placement or read
name.

## Haplotype And Phase

### Haplotype blocks

Reads are grouped into haplotype blocks using HP tags when available. HP1, HP2,
and unassigned reads are displayed in separate lanes, with pink horizontal
dividers separating haplotype groups. The haplotype label appears at the left
side of each block and uses the shared Orographer haplotype color palette.

The row order is derived from the available reads in the sample. In Paraphase
mode, the coverage calculation tracks HP1 and HP2 separately, so those
haplotypes remain central to the display even when unassigned reads are also
present. Reads without usable HP tags are placed in an unassigned block rather
than being forced into a phased lane.

### Phase set delimiters

Phase set delimiters are vertical pink markers drawn inside haplotype blocks
using PS tag information from the alignments. They mark visible phase-set starts
and observed phase-block extents for the displayed evidence. The markers are
limited to the y-range of the corresponding haplotype block, so phase
boundaries are read in the context of the haplotype evidence they describe.

Phase-set starts and ends are included only when the tagged segment overlaps the
visible region and the haplotype tag is usable. Reads without HP or PS tags do
not contribute phase-set delimiter markers. The markers can be hidden with the
"Hide phaseset markers" checkbox.

### Coverage by haplotype

Paraphase coverage includes total depth and haplotype-specific HP1 and HP2
series. The total series includes all qualifying alignments, including
unassigned reads. Haplotype-specific series are drawn for reads whose HP tags
match the tracked haplotypes. This makes allele balance and haplotype-specific
dropouts visible above the read alignments.

Coverage is computed in fixed 10 bp bins. Adjacent aligned blocks are merged
when the intervening CIGAR gap is smaller than 20 bp, smoothing over small
indels while preserving larger gaps as visible dips. This keeps the coverage
track focused on regional support rather than every tiny CIGAR interruption.

## Coverage And Variants

### Coverage plots

Each sample row can include a coverage plot linked to the same genomic x-axis
as the read track. The filled total-depth series summarizes all qualifying
alignments across the requested interval. Haplotype-specific series are drawn as
colored lines, and clicking a haplotype coverage line can emphasize that series
and its label.

Coverage is based on reference-consuming aligned bases. Inserted sequence does
not add reference-position coverage at the insertion site. Large deletion or
skip gaps are preserved as coverage dips, while small CIGAR gaps below the
smoothing threshold may be bridged in the displayed coverage series.

### Read-level variant annotations

When reference sequence is available, Orographer compares aligned read bases to
the reference and extracts mismatches, insertions, and deletions for displayed
segments. Mismatches appear as colored markers using IGV-style alternate-base
colors. Insertions appear as dark markers sized by inserted sequence length,
and deletions appear as grey line segments over the deleted reference span.

One-base insertions and deletions can be hidden with the "Hide 1bp INDELs"
checkbox. When the view is zoomed in below 1 kb, mismatch and insertion markers
switch to text labels so base-level evidence can be read directly. At wider
views, compact markers are used to keep dense regions legible.

### VCF track

When a VCF is provided, variants overlapping the displayed region are shown in a
separate VCF track linked to the same genomic x-axis as the read track. Variants
are drawn as small triangle markers at their genomic position. Deletions are
placed at the center of the deleted interval.

VCF marker color reflects the parsed variant type: SNPs use the alternate-base
color, insertions use the insertion color, deletions use the deletion color,
and unknown variants use a neutral grey. Clicking a VCF marker opens a details
modal with the sample label, coordinates, variant type, alternate allele, and
haplotype values parsed from the VCF when present.

## Interaction And Navigation

### Navigation features

The coordinate controls let you zoom, pan, and jump to a specific interval.
The plus and minus buttons zoom in and out around the current view. The side
arrow buttons pan left or right. Buttons are disabled when their action would
not change the view, such as panning at the edge of the original region or
zooming out when already fully zoomed out.

The coordinate input accepts genomic coordinates in `chrom:start-end` format,
including comma-formatted positions such as `chr8:127,000,000-146,000,000`.
The current view size is shown below the navigation buttons.

The Bokeh toolbar also provides built-in plot navigation. Dragging a box across
the plot zooms the view to the selected rectangle. The reset control restores
the original x-range and y-range for the plot and clears transient view
changes.

### Read selection options

Clicking a read alignment selects all displayed segments for that read within
the plotted sample row. Selected reads keep strong visual emphasis, while
unselected reads are faded so the selected molecule is easier to trace across
its alignment blocks.

The global "Select reads" control can be used to search for read names
directly. This is useful when reviewing a known molecule or when tracing a read
from external analysis output. The "Clear selected" button clears the current
read selection and is disabled when no reads are selected.

### Display controls

The "Show cursor guide" checkbox controls the IGV-style vertical cursor guide.
When enabled, moving the mouse over a plot shows a vertical guide line across
linked tracks so positions can be compared precisely between coverage, reads,
variants, genes, and the genomic axis.

Additional display checkboxes reduce clutter:

- "Hide 1bp INDELs" hides one-base insertion and deletion annotations.
- "Hide algn numbers" hides alignment number labels.
- "Hide phaseset markers" hides vertical phase-set delimiter lines.
