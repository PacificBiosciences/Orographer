"""Reference self-identity dotplot computation at long-read scale."""

import logging
import math
from dataclasses import dataclass

import numpy as np

logger = logging.getLogger(__name__)

MAX_DOTPLOT_REGION_BP = 10_000_000
_RESOLUTION = 512
_KMER = 15
_STRIDE = 50
_MIN_WINDOW_BP = 1_000
_MAX_WINDOW_BP = 10_000
_MIN_IDENTITY = 0.5
_MIN_SHARED_SEEDS = 3
_MAX_SEED_BIN_FRACTION = 0.05
_MAX_CANDIDATE_ALIGNMENTS = 20_000
_MAX_ALIGNMENT_SHIFT_BP = 200
_ALIGNMENT_SHIFT_STEP_BP = 25
_MIN_ALIGNED_FRACTION = 0.3

_RC_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")
_MATCH_COLOR = np.array([31, 78, 121, 255], dtype=np.uint8)
_BACKGROUND_COLOR = np.array([255, 255, 255, 255], dtype=np.uint8)


@dataclass(frozen=True)
class _Window:
    start: int
    end: int


def _canonical(kmer: str) -> str:
    """Return the lexicographically smaller of kmer and its reverse complement."""
    rc = kmer.translate(_RC_TABLE)[::-1]
    return kmer if kmer <= rc else rc


def _reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a sequence."""
    return sequence.translate(_RC_TABLE)[::-1]


def _has_called_bases(sequence: str) -> bool:
    """Return true when the sequence window contains at least one called base."""
    return any(base in "ACGTacgt" for base in sequence)


def _choose_window_size(sequence_length: int, max_windows: int) -> int:
    """Choose a read-scale window size that keeps the dotplot matrix bounded."""
    if sequence_length <= _MIN_WINDOW_BP:
        return max(1, sequence_length)
    target_size = math.ceil(2 * sequence_length / max(1, max_windows))
    return min(_MAX_WINDOW_BP, max(_MIN_WINDOW_BP, target_size))


def _build_windows(sequence_length: int, max_windows: int) -> list[_Window]:
    """Build half-overlapping windows across the sequence."""
    if sequence_length <= 0:
        return []
    window_size = _choose_window_size(sequence_length, max_windows)
    step = max(1, window_size // 2)
    windows = []
    for start in range(0, sequence_length, step):
        end = min(sequence_length, start + window_size)
        if windows and end - start <= window_size // 4:
            break
        windows.append(_Window(start=start, end=end))
        if end == sequence_length:
            break
    return windows


def _iter_window_seeds(sequence: str, window: _Window, k: int, stride: int) -> set[str]:
    """Collect canonical k-mer seeds from one window."""
    if window.end - window.start < k:
        return set()
    seeds = set()
    last_start = window.end - k
    step = max(1, stride)
    positions = list(range(window.start, last_start + 1, step))
    positions.extend(range(last_start, window.start - 1, -step))
    for pos in positions:
        kmer = sequence[pos : pos + k]
        if "N" not in kmer and "n" not in kmer:
            seeds.add(_canonical(kmer))
    return seeds


def _index_window_seeds(seed_sets: list[set[str]]) -> dict[str, list[int]]:
    """Return an inverted seed index from canonical seed to window indices."""
    seed_to_windows: dict[str, list[int]] = {}
    for window_index, seeds in enumerate(seed_sets):
        for seed in seeds:
            seed_to_windows.setdefault(seed, []).append(window_index)
    return seed_to_windows


def _candidate_pairs(seed_sets: list[set[str]]) -> list[tuple[int, int]]:
    """Return self-pairs plus seeded off-diagonal pairs that justify scoring."""
    seed_to_windows = _index_window_seeds(seed_sets)
    max_seed_bins = max(8, int(len(seed_sets) * _MAX_SEED_BIN_FRACTION))
    pair_counts: dict[tuple[int, int], int] = {}
    for window_indices in seed_to_windows.values():
        if len(window_indices) <= 1 or len(window_indices) > max_seed_bins:
            continue
        for left_offset, left_index in enumerate(window_indices[:-1]):
            for right_index in window_indices[left_offset + 1 :]:
                pair = (left_index, right_index)
                pair_counts[pair] = pair_counts.get(pair, 0) + 1

    candidates = [
        (left_index, right_index, shared_count)
        for (left_index, right_index), shared_count in pair_counts.items()
        if shared_count >= _MIN_SHARED_SEEDS
    ]
    candidates.sort(key=lambda item: item[2], reverse=True)
    self_pairs = [(index, index) for index in range(len(seed_sets))]
    off_diagonal_pairs = [
        (left_index, right_index) for left_index, right_index, _count in candidates
    ][:_MAX_CANDIDATE_ALIGNMENTS]
    return [*self_pairs, *off_diagonal_pairs]


def _shifted_identity(left_array: np.ndarray, right_array: np.ndarray, shift: int) -> float:
    """Return matches over max window length for one relative shift."""
    if shift >= 0:
        left_slice = left_array[shift:]
        right_slice = right_array[: len(left_slice)]
    else:
        right_slice = right_array[-shift:]
        left_slice = left_array[: len(right_slice)]
    aligned_length = min(len(left_slice), len(right_slice))
    denominator = max(len(left_array), len(right_array), 1)
    if aligned_length / denominator < _MIN_ALIGNED_FRACTION:
        return 0.0
    matches = np.count_nonzero(left_slice[:aligned_length] == right_slice[:aligned_length])
    return float(matches / denominator)


def _sequence_match_ratio(left_sequence: str, right_sequence: str) -> float:
    """Return a bounded collinear identity score for two long-read-scale windows."""
    left_array = np.frombuffer(left_sequence.encode("ascii"), dtype=np.uint8)
    right_array = np.frombuffer(right_sequence.encode("ascii"), dtype=np.uint8)
    max_shift = min(_MAX_ALIGNMENT_SHIFT_BP, max(len(left_array), len(right_array)) - 1)
    shifts = list(range(-max_shift, max_shift + 1, _ALIGNMENT_SHIFT_STEP_BP))
    if 0 not in shifts:
        shifts.append(0)
    return max(_shifted_identity(left_array, right_array, shift) for shift in shifts)


def _window_identity(left_sequence: str, right_sequence: str) -> float:
    """Score direct and inverted window similarity and return the better identity."""
    if not _has_called_bases(left_sequence) or not _has_called_bases(right_sequence):
        return 0.0
    if left_sequence == right_sequence:
        return 1.0
    direct_score = _sequence_match_ratio(left_sequence, right_sequence)
    inverted_score = _sequence_match_ratio(left_sequence, _reverse_complement(right_sequence))
    return max(direct_score, inverted_score)


def compute_self_identity(
    sequence: str,
    resolution: int = _RESOLUTION,
    k: int = _KMER,
    stride: int = _STRIDE,
) -> np.ndarray:
    """Return a read-scale self-identity matrix for one concatenated reference sequence.

    The matrix uses half-overlapping windows rather than fixed base-scale bins. Candidate
    off-diagonal pairs are discovered with canonical k-mer seeds, then scored with bounded
    collinear identity so that a tiny shared motif does not masquerade as long-read-scale
    similarity.
    """
    sequence = sequence.upper()
    windows = _build_windows(len(sequence), resolution)
    if not windows:
        return np.zeros((0, 0), dtype=np.float32)

    matrix = np.zeros((len(windows), len(windows)), dtype=np.float32)
    seed_sets = [_iter_window_seeds(sequence, window, k, stride) for window in windows]
    for left_index, right_index in _candidate_pairs(seed_sets):
        left_window = windows[left_index]
        right_window = windows[right_index]
        left_sequence = sequence[left_window.start : left_window.end]
        right_sequence = sequence[right_window.start : right_window.end]
        identity = _window_identity(left_sequence, right_sequence)
        if identity >= _MIN_IDENTITY:
            matrix[left_index, right_index] = identity
            matrix[right_index, left_index] = identity

    return matrix


def _matrix_to_rgba(matrix: np.ndarray) -> np.ndarray:
    """Convert a 0..1 similarity matrix to uint32 RGBA for Bokeh image_rgba."""
    h, w = matrix.shape
    intensity = np.clip(matrix.astype(np.float32), 0.0, 1.0)
    rgba_float = _BACKGROUND_COLOR.astype(np.float32) + intensity[:, :, None] * (
        _MATCH_COLOR.astype(np.float32) - _BACKGROUND_COLOR.astype(np.float32)
    )
    rgba = rgba_float.astype(np.uint8)

    rgba = np.flipud(rgba)
    return rgba.view(np.uint32).reshape(h, w)


def compute_dotplot_image(
    sequence: str,
    resolution: int = _RESOLUTION,
    k: int = _KMER,
) -> np.ndarray:
    """Compute and return a uint32 RGBA image array."""
    matrix = compute_self_identity(sequence, resolution=resolution, k=k)
    return _matrix_to_rgba(matrix)
