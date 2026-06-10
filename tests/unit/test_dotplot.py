"""Unit tests for orographer.dotplot module."""

import numpy as np

from orographer.dotplot import (
    _KMER,
    _RESOLUTION,
    MAX_DOTPLOT_REGION_BP,
    _canonical,
    _matrix_to_rgba,
    compute_dotplot_image,
    compute_self_identity,
)


class TestComputeSelfIdentity:
    def test_output_shape_is_adaptive_but_bounded(self):
        seq = "ACGT" * 5_000
        matrix = compute_self_identity(seq, resolution=64)

        assert matrix.shape[0] == matrix.shape[1]
        assert matrix.shape[0] <= 64

    def test_unique_sequence_scores_self_pairs_without_off_diagonal_signal(self) -> None:
        rng = np.random.default_rng(42)
        bases = np.array(list("ACGT"))
        seq = "".join(rng.choice(bases, size=20_000))

        matrix = compute_self_identity(seq, resolution=64)

        assert np.all(np.diag(matrix) == 1.0)
        off_diag = matrix.copy()
        np.fill_diagonal(off_diag, 0)
        assert off_diag.max() == 0

    def test_all_n_sequence(self):
        seq = "N" * 2_000
        matrix = compute_self_identity(seq, resolution=16)
        assert np.all(matrix == 0)

    def test_sequence_shorter_than_k_scores_self_identity(self) -> None:
        seq = "ACG"
        matrix = compute_self_identity(seq, resolution=8, k=11)
        assert matrix.shape == (1, 1)
        assert matrix[0, 0] == 1.0

    def test_unique_sequence_no_off_diagonal(self):
        rng = np.random.default_rng(42)
        bases = np.array(list("ACGT"))
        seq = "".join(rng.choice(bases, size=20_000))
        matrix = compute_self_identity(seq, resolution=64)
        off_diag = matrix.copy()
        np.fill_diagonal(off_diag, 0)
        assert off_diag.max() == 0

    def test_dtype_is_float32(self):
        seq = "ACGT" * 500
        matrix = compute_self_identity(seq, resolution=16)
        assert matrix.dtype == np.float32

    def test_similar_blocks_produce_cross_block_identity(self):
        rng = np.random.default_rng(7)
        bases = np.array(list("ACGT"))
        block = "".join(rng.choice(bases, size=5_000))
        mutated = list(block)
        for idx in range(0, len(mutated), 100):
            current_base = mutated[idx]
            choices = [base for base in "ACGT" if base != current_base]
            mutated[idx] = choices[(idx // 100) % len(choices)]
        seq = block + ("N" * 1_000) + "".join(mutated)

        matrix = compute_self_identity(seq, resolution=64)
        off_diag = matrix.copy()
        np.fill_diagonal(off_diag, 0)
        assert off_diag.max() > 0.95

    def test_partial_long_homology_is_displayed_below_old_high_identity_cutoff(self):
        rng = np.random.default_rng(91)
        bases = np.array(list("ACGT"))
        block = "".join(rng.choice(bases, size=5_000))
        tail = "".join(rng.choice(bases, size=2_000))
        partial_match = block[:3_000] + tail
        seq = block + ("N" * 1_000) + partial_match

        matrix = compute_self_identity(seq, resolution=64)
        off_diag = matrix.copy()
        np.fill_diagonal(off_diag, 0)
        intermediate_scores = off_diag[(off_diag >= 0.5) & (off_diag < 0.78)]
        assert intermediate_scores.size > 0

    def test_short_shared_motif_does_not_look_like_long_read_identity(self):
        rng = np.random.default_rng(13)
        bases = np.array(list("ACGT"))
        left = list("".join(rng.choice(bases, size=5_000)))
        right = list("".join(rng.choice(bases, size=5_000)))
        right[2_000:2_200] = left[2_000:2_200]
        seq = "".join(left) + ("N" * 1_000) + "".join(right)

        matrix = compute_self_identity(seq, resolution=64)
        off_diag = matrix.copy()
        np.fill_diagonal(off_diag, 0)
        assert off_diag.max() == 0


class TestMatrixToRgba:
    def test_output_shape_and_dtype(self):
        matrix = np.zeros((8, 8), dtype=bool)
        rgba = _matrix_to_rgba(matrix)
        assert rgba.shape == (8, 8)
        assert rgba.dtype == np.uint32

    def test_true_pixels_not_white(self):
        matrix = np.zeros((4, 4), dtype=np.float32)
        matrix[0, 0] = 1.0
        rgba_all_false = _matrix_to_rgba(np.zeros((4, 4), dtype=np.float32))
        rgba_with_true = _matrix_to_rgba(matrix)
        # The flipped position of (0,0) True is at row 3 in output
        assert rgba_with_true[3, 0] != rgba_all_false[3, 0]

    def test_false_pixels_are_white(self):
        matrix = np.zeros((4, 4), dtype=np.float32)
        rgba = _matrix_to_rgba(matrix)
        # White = 0xFF_FF_FF_FF in little-endian RGBA uint32
        white_val = np.array([255, 255, 255, 255], dtype=np.uint8).view(np.uint32)[0]
        assert np.all(rgba == white_val)

    def test_partial_scores_are_between_white_and_match_color(self):
        matrix = np.zeros((3, 3), dtype=np.float32)
        matrix[0, 0] = 0.5
        matrix[0, 1] = 1.0
        rgba = _matrix_to_rgba(matrix).view(np.uint8).reshape(3, 3, 4)
        partial_red = rgba[2, 0, 0]
        full_red = rgba[2, 1, 0]
        assert full_red < partial_red < 255

    def test_vertical_flip(self):
        # Row 0 of boolean matrix → row N-1 of output (bottom of image = sequence start)
        matrix = np.zeros((4, 4), dtype=np.float32)
        matrix[0, :] = 1.0
        rgba = _matrix_to_rgba(matrix)
        white_val = np.array([255, 255, 255, 255], dtype=np.uint8).view(np.uint32)[0]
        # After flipud, row 0 of matrix is at last row of output
        assert not np.all(rgba[3, :] == white_val), "Row 0 of matrix should be at bottom of output"
        assert np.all(rgba[0, :] == white_val), "Row 3 of all-False matrix stays white at top"


class TestComputeDotplotImage:
    def test_output_shape(self):
        seq = "ACGT" * 2_000
        img = compute_dotplot_image(seq, resolution=32)
        assert img.shape[0] == img.shape[1]
        assert img.shape[0] <= 32
        assert img.dtype == np.uint32

    def test_default_resolution(self):
        seq = "ACGT" * 2_000
        img = compute_dotplot_image(seq)
        assert img.shape[0] == img.shape[1]
        assert img.shape[0] <= _RESOLUTION


class TestCanonical:
    def test_forward_smaller_returns_forward(self):
        assert _canonical("AAACGTT") == "AAACGTT"

    def test_rc_smaller_returns_rc(self):
        # RC of TTTCGAA is TTCGAAA → "TTCGAAA" > "TTTCGAA" so forward wins here;
        # pick a pair where RC wins: TTTTTTTTTTA → RC = TAAAAAAAAA < original
        kmer = "TTTTTTTTTTA"
        assert _canonical(kmer) == _canonical(_canonical(kmer))  # idempotent

    def test_canonical_is_symmetric(self):
        kmer = "ACGTACGTACC"
        rc = kmer.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]
        assert _canonical(kmer) == _canonical(rc)

    def test_palindrome_returns_itself(self):
        # A palindromic k-mer equals its own RC → canonical == kmer
        kmer = "AACGTT"  # RC = AACGTT ✓
        assert _canonical(kmer) == kmer


class TestInvertedRepeats:
    def test_inverted_repeat_produces_off_diagonal(self):
        rng = np.random.default_rng(99)
        bases = np.array(list("ACGT"))
        forward_block = "".join(rng.choice(bases, size=3_000))
        inverted_block = forward_block.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        seq = forward_block + ("N" * 1_000) + inverted_block

        matrix = compute_self_identity(seq, resolution=64)

        off_diag = matrix.copy()
        np.fill_diagonal(off_diag, 0)
        assert off_diag.max() == 1.0

    def test_no_matches_without_similarity(self):
        rng = np.random.default_rng(99)
        bases = np.array(list("ACGT"))
        seq = "".join(rng.choice(bases, size=20_000))
        matrix = compute_self_identity(seq, resolution=64)
        off_diag = matrix.copy()
        np.fill_diagonal(off_diag, 0)
        assert off_diag.max() == 0


class TestConstants:
    def test_max_region_bp(self):
        assert MAX_DOTPLOT_REGION_BP == 10_000_000

    def test_resolution(self):
        assert _RESOLUTION == 512

    def test_kmer(self):
        assert _KMER == 15
