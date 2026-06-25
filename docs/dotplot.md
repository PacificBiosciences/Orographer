# Reference Self-Identity Dotplot

## Table of Contents

- [Overview](#overview)
- [Algorithm](#algorithm)
  - [K-mer extraction and canonicalization](#k-mer-extraction-and-canonicalization)
  - [Bin assignment](#bin-assignment)
  - [Counting and thresholding](#counting-and-thresholding)
  - [Low-complexity filtering](#low-complexity-filtering)
  - [Adaptive threshold](#adaptive-threshold)
- [Interpreting the Display](#interpreting-the-display)
  - [The diagonal](#the-diagonal)
  - [Off-diagonal black regions](#off-diagonal-black-regions)
  - [Inverted repeats](#inverted-repeats)
  - [White regions](#white-regions)
  - [Gray regions (N-masked)](#gray-regions-n-masked)
- [Multi-Region Mode](#multi-region-mode)
  - [What the concatenated view shows](#what-the-concatenated-view-shows)
  - [Pink boundary lines](#pink-boundary-lines)
  - [Axis labels](#axis-labels)
- [Limitations](#limitations)

---

## Overview

The reference self-identity dotplot compares a genomic reference sequence against
itself to reveal internal structure: tandem repeats, segmental duplications,
inverted repeats, and assembly gaps. This dotplot is created by exact k-mer
matching at word length 15, aggregated into a 512 × 512 bin grid. Each dark
pixel in the plot represents a region where two parts of the sequence share
many identical 15-mers — specifically, regions at or above ~90% sequence
identity, where HiFi reads are at risk of misalignment.

---

## Algorithm

### K-mer extraction and canonicalization

Every 15-mer in the reference sequence is extracted and converted to its
**canonical form**: the lexicographically smaller of the 15-mer and its reverse
complement. Both a 15-mer and its reverse complement map to the same canonical
key, which makes the dotplot sensitive to inverted repeats. K-mers containing
ambiguous bases (N) are skipped.

### Bin assignment

The sequence of length N is divided into a 512 × 512 grid. Each position i maps
to bin `i * 512 // N`. Every canonical k-mer records which bin it was observed
in.

### Adaptive threshold

The minimum number of shared k-mers required to call a bin pair a match scales
with bin size:

```
min_shared = max(3, bin_size // 5)
```

where `bin_size = N / 512`. For a 2 Mb region, each bin is roughly 3,906 bp and
`min_shared ≈ 781`. For a 136 kb region, each bin is roughly 266 bp and
`min_shared ≈ 53`. For a 33 kb region, each bin is roughly 65 bp and
`min_shared ≈ 13`.


### Counting and thresholding

For each canonical k-mer, the set of bins it occupies is collected. Every pair
of bins that share the same k-mer gets one count increment. After processing all
k-mers, bin pair (i, j) has a count equal to the number of distinct canonical
k-mers shared between those two bins.

A bin pair is marked as a match (black pixel) when its count meets or exceeds the
adaptive minimum threshold (see below). The main diagonal is always marked
wherever at least one valid k-mer exists, because any bin is perfectly identical
to itself.

### Low-complexity filtering

K-mers that appear in more than 128 distinct bins are discarded before counting.
A k-mer spread across more than one quarter of the 512-bin grid is a
low-complexity sequence (e.g., a simple dinucleotide repeat).

A consequence of this filter is that highly repetitive regions — pericentromeric
satellite arrays, for example — may produce few off-diagonal matches even when
the region is internally repetitive. In those cases the diagonal may still be
marked, but the off-diagonal structure is suppressed.


---

## Interpreting the Display

### The diagonal

The main diagonal (bottom-left to top-right) represents the sequence compared
against itself. Its continuity is a quick indicator of assembly quality — 
interruptions indicate assembly gaps.

### Off-diagonal black regions

A black region at position (i, j) means that bins i and j share many identical
15-mers. Common patterns:

| Pattern | Cause |
|---------|-------|
| Parallel lines offset from the diagonal | Tandem repeat: the same sequence unit appears multiple times at regular intervals |
| Filled rectangle off the diagonal | Segmental duplication: a large block appears at two locations |
| Isolated spots | Short shared motifs |
| Grid of evenly-spaced dots | Periodic tandem repeat with period smaller than one bin |

### Inverted repeats

An inverted repeat appears as signal on the **anti-diagonal** (running from
top-right to bottom-left).

### White regions

A white pixel at position (i, j) means those two bins share fewer canonical
k-mers than the threshold. For off-diagonal positions in a
self-identity plot, white is expected unless the region contains internal
repeats.

### Gray regions (N-masked)

Light gray pixels indicate that one or both bins in a comparison contain no
valid sequence (all `N`s).

---

## Multi-Region Mode

When more than one coordinate region is provided, the dotplot concatenates the
reference sequences end-to-end and displays a single combined matrix.

### Pink boundary lines

Pink lines mark the boundaries between concatenated regions. Each line is
placed at the first bin edge that falls at or after the sequence boundary, so it
aligns with the visual pixel grid rather than falling inside a mixed bin. In
some cases a pink line may appear slightly offset from where you expect the
boundary to be; this is correct behavior when the boundary falls in the middle
of a bin.

---

## Limitations

**Satellite and pericentromeric repeats.** The low-complexity filter (k-mers in
more than 128 bins) removes most satellite k-mers. In regions dominated by
satellite DNA, nearly all k-mers are suppressed and the off-diagonal structure
disappears. The diagonal remains marked wherever the reference has sequence.

**Resolution.** The 512-bin grid means that features smaller than one bin
(roughly `sequence_length / 512` base pairs) may not be resolved. At 136 kb,
each bin covers roughly 266 bp; a 100 bp repeat unit would fall within a
single bin and appear as a solid block rather than a periodic pattern.