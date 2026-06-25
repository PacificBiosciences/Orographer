"""Reference self-identity dotplot computation using gepard-style k-mer exact matching."""

import logging

import numpy as np

logger = logging.getLogger(__name__)

MAX_DOTPLOT_REGION_BP = 10_000_000
_RESOLUTION = 512
_KMER = 15
# k-mers appearing in more than this many unique bins are suppressed as low-complexity noise.
_MAX_KMER_BINS = _RESOLUTION // 4
# Adaptive threshold: require at least this many distinct k-mers to confirm a bin pair.
# At k=15, the expected fraction of shared k-mers at identity p is p^15.
# Requiring bin_size // 5 shared k-mers (20% of the bin) targets ~90%+ identity —
# the range where HiFi aligners (minimap2 k=19 seeds) can misplace reads.
_MIN_SHARED_KMERS_FLOOR = 3  # minimum regardless of bin size
_MIN_SHARED_KMERS_PER_BP = 5  # one required k-mer per this many base pairs of bin size

_RC_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")
_MATCH_COLOR = np.array([0, 0, 0, 255], dtype=np.uint8)
_BACKGROUND_COLOR = np.array([255, 255, 255, 255], dtype=np.uint8)
_N_MASKED_COLOR = np.array([216, 216, 216, 255], dtype=np.uint8)


def _canonical(kmer: str) -> str:
    """Return the lexicographically smaller of kmer and its reverse complement."""
    rc = kmer.translate(_RC_TABLE)[::-1]
    return kmer if kmer <= rc else rc


def _reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a sequence."""
    return sequence.translate(_RC_TABLE)[::-1]


def compute_self_identity(
    sequence: str,
    resolution: int = _RESOLUTION,
    k: int = _KMER,
) -> np.ndarray:
    """Return a binary self-identity matrix via gepard-style k-mer exact matching.

    Each position i in the sequence is binned to resolution×resolution coordinates.
    Canonical k-mers (handling inverted repeats) that appear in 2–_MAX_KMER_BINS
    distinct bins contribute a count to each shared bin pair. Only pairs confirmed
    by the adaptive minimum (scaled to bin size) are marked (1.0). Short sequences
    have smaller bins with fewer k-mers, so the threshold is scaled down to avoid
    over-filtering weak but real homology. The main diagonal is always marked wherever
    at least one valid (non-N) k-mer is present. Returns a float32 matrix of shape
    (resolution, resolution).
    """
    seq = sequence.upper()
    n = len(seq)
    if n < k:
        return np.zeros((resolution, resolution), dtype=np.float32)

    bin_size = n / resolution
    min_shared = max(_MIN_SHARED_KMERS_FLOOR, int(bin_size / _MIN_SHARED_KMERS_PER_BP))

    occupied_bins: set[int] = set()
    kmer_bins: dict[str, set[int]] = {}

    for i in range(n - k + 1):
        kmer = seq[i : i + k]
        if "N" in kmer:
            continue
        canon = _canonical(kmer)
        bin_i = i * resolution // n
        occupied_bins.add(bin_i)
        kmer_bins.setdefault(canon, set()).add(bin_i)

    counts = np.zeros((resolution, resolution), dtype=np.uint16)
    for bins in kmer_bins.values():
        if len(bins) < 2 or len(bins) > _MAX_KMER_BINS:
            continue
        b = np.fromiter(bins, dtype=np.int32, count=len(bins))
        counts[np.ix_(b, b)] += 1

    matrix = (counts >= min_shared).astype(np.float32)

    for bin_i in occupied_bins:
        matrix[bin_i, bin_i] = 1.0

    # Bins with no valid k-mers get their entire row and column set to -1.0 so the
    # N-masked region is visible as a band, not just a 1-pixel diagonal dot.
    n_masked = np.array(
        [i for i in range(resolution) if i not in occupied_bins], dtype=np.intp
    )
    if n_masked.size > 0:
        matrix[n_masked, :] = -1.0
        matrix[:, n_masked] = -1.0

    return matrix


def _matrix_to_rgba(matrix: np.ndarray) -> np.ndarray:
    """Convert a similarity matrix to uint32 RGBA for Bokeh image_rgba.

    Values: 1.0 = black (match), 0.0 = white (no match), -1.0 = light gray (N-masked).
    """
    h, w = matrix.shape
    is_n_masked = matrix < 0.0
    intensity = np.clip(matrix, 0.0, 1.0)
    rgba_float = _BACKGROUND_COLOR.astype(np.float32) + intensity[:, :, None] * (
        _MATCH_COLOR.astype(np.float32) - _BACKGROUND_COLOR.astype(np.float32)
    )
    rgba = rgba_float.astype(np.uint8)
    rgba[is_n_masked] = _N_MASKED_COLOR
    rgba = np.flipud(rgba)
    return rgba.view(np.uint32).reshape(h, w)


def compute_repeat_density(matrix: np.ndarray) -> np.ndarray:
    """Return per-bin off-diagonal match count as a float32 array of shape (resolution,).

    Sums each column of the matrix, excluding the diagonal and N-masked sentinels (-1.0).
    The result is a 1-D measure of how many other bins each bin matches — useful as a
    compact repeat-density track in the main view.
    """
    clean = np.where(matrix < 0.0, 0.0, matrix)
    np.fill_diagonal(clean, 0.0)
    return clean.sum(axis=0).astype(np.float32)


def compute_dotplot_image(
    sequence: str,
    resolution: int = _RESOLUTION,
    k: int = _KMER,
) -> np.ndarray:
    """Compute and return a uint32 RGBA image array."""
    matrix = compute_self_identity(sequence, resolution=resolution, k=k)
    return _matrix_to_rgba(matrix)
