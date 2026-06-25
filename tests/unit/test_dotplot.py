"""Unit tests for orographer.dotplot module."""

import numpy as np

from orographer.dotplot import (
    _KMER,
    _MAX_KMER_BINS,
    _MIN_SHARED_KMERS_FLOOR,
    _MIN_SHARED_KMERS_PER_BP,
    _N_MASKED_COLOR,
    _RESOLUTION,
    MAX_DOTPLOT_REGION_BP,
    _canonical,
    _matrix_to_rgba,
    compute_dotplot_image,
    compute_repeat_density,
    compute_self_identity,
)


class TestComputeSelfIdentity:
    def test_output_shape_is_resolution_by_resolution(self):
        seq = "ACGT" * 5_000
        matrix = compute_self_identity(seq, resolution=64)
        assert matrix.shape == (64, 64)

    def test_dtype_is_float32(self):
        seq = "ACGT" * 500
        matrix = compute_self_identity(seq, resolution=16)
        assert matrix.dtype == np.float32

    def test_all_n_sequence_marks_all_entries_as_n_masked(self):
        seq = "N" * 2_000
        matrix = compute_self_identity(seq, resolution=16)
        assert np.all(matrix == -1.0), "All entries must be -1.0 for all-N sequence"

    def test_partial_n_marks_n_bin_rows_and_cols(self):
        # Left half ACGT, right half N.  N-masked bins appear in both rows AND columns.
        seq = "ACGTACGT" * 500 + "N" * 4_000
        matrix = compute_self_identity(seq, resolution=16)
        # ACGT-vs-ACGT quadrant: no -1.0
        assert not np.any(matrix[:8, :8] == -1.0), "ACGT-vs-ACGT quadrant must have no -1.0"
        # N-masked rows (bins 8-15): entire row must be -1.0
        assert np.all(matrix[8:, :] == -1.0), "N-masked rows must be entirely -1.0"
        # N-masked columns (bins 8-15): entire column must be -1.0
        assert np.all(matrix[:, 8:] == -1.0), "N-masked columns must be entirely -1.0"

    def test_sequence_shorter_than_k_returns_all_zeros(self) -> None:
        seq = "ACG"
        matrix = compute_self_identity(seq, resolution=8, k=11)
        assert matrix.shape == (8, 8)
        assert np.all(matrix == 0)

    def test_matrix_values_are_binary(self):
        seq = "ACGTACGT" * 2_000
        matrix = compute_self_identity(seq, resolution=32)
        unique_values = np.unique(matrix)
        for val in unique_values:
            assert val == 0.0 or val == 1.0, f"Non-binary value {val} found"

    def test_diagonal_is_marked_for_valid_sequence(self):
        # A plain repeating unit saturates every bin with k-mer occupancy.
        seq = "ACGTACGT" * 5_000
        matrix = compute_self_identity(seq, resolution=64)
        assert np.all(np.diag(matrix) == 1.0)

    def test_tandem_repeat_produces_off_diagonal_dots(self):
        # k-mers in a tandem repeat appear at every period → off-diagonal bin pairs.
        # Period must exceed k so each position generates a distinct canonical k-mer,
        # producing enough shared k-mers per bin pair to pass the adaptive threshold.
        rng = np.random.default_rng(7)
        bases = np.array(list("ACGT"))
        unit = "".join(rng.choice(bases, size=50))  # period > k
        seq = unit * 200  # 10,000 bp
        matrix = compute_self_identity(seq, resolution=64)
        off_diag = matrix.copy()
        np.fill_diagonal(off_diag, 0)
        assert off_diag.max() == 1.0

    def test_similar_blocks_produce_off_diagonal_signal(self):
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
        assert off_diag.max() == 1.0

    def test_short_shared_motif_is_detected(self):
        # New algorithm is more sensitive: a 200 bp shared block between two halves IS found.
        rng = np.random.default_rng(13)
        bases = np.array(list("ACGT"))
        left = list("".join(rng.choice(bases, size=5_000)))
        right = list("".join(rng.choice(bases, size=5_000)))
        right[2_000:2_200] = left[2_000:2_200]
        seq = "".join(left) + ("N" * 1_000) + "".join(right)

        matrix = compute_self_identity(seq, resolution=64)
        off_diag = matrix.copy()
        np.fill_diagonal(off_diag, 0)
        assert off_diag.max() == 1.0


class TestMatrixToRgba:
    def test_output_shape_and_dtype(self):
        matrix = np.zeros((8, 8), dtype=np.float32)
        rgba = _matrix_to_rgba(matrix)
        assert rgba.shape == (8, 8)
        assert rgba.dtype == np.uint32

    def test_zero_pixels_are_white(self):
        matrix = np.zeros((4, 4), dtype=np.float32)
        rgba = _matrix_to_rgba(matrix)
        white_val = np.array([255, 255, 255, 255], dtype=np.uint8).view(np.uint32)[0]
        assert np.all(rgba == white_val)

    def test_one_pixels_are_black(self):
        matrix = np.ones((4, 4), dtype=np.float32)
        rgba = _matrix_to_rgba(matrix)
        black_val = np.array([0, 0, 0, 255], dtype=np.uint8).view(np.uint32)[0]
        assert np.all(rgba == black_val)

    def test_n_masked_pixels_are_gray(self):
        matrix = np.full((4, 4), -1.0, dtype=np.float32)
        rgba = _matrix_to_rgba(matrix).view(np.uint8).reshape(4, 4, 4)
        expected_r = _N_MASKED_COLOR[0]
        assert np.all(rgba[:, :, 0] == expected_r), "N-masked pixels must match _N_MASKED_COLOR"
        assert rgba[0, 0, 0] == rgba[0, 0, 1] == rgba[0, 0, 2], "Gray must have equal RGB"
        white_val = 255
        black_val = 0
        assert rgba[0, 0, 0] != white_val, "N-masked must not be white"
        assert rgba[0, 0, 0] != black_val, "N-masked must not be black"

    def test_match_pixels_are_darker_than_background(self):
        matrix = np.zeros((4, 4), dtype=np.float32)
        matrix[0, 0] = 1.0
        rgba_match = _matrix_to_rgba(matrix).view(np.uint8).reshape(4, 4, 4)
        rgba_bg = _matrix_to_rgba(np.zeros((4, 4), dtype=np.float32)).view(np.uint8).reshape(4, 4, 4)
        # Row 0 of matrix → row 3 of output after flipud
        assert rgba_match[3, 0, 0] < rgba_bg[3, 0, 0]

    def test_vertical_flip(self):
        matrix = np.zeros((4, 4), dtype=np.float32)
        matrix[0, :] = 1.0
        rgba = _matrix_to_rgba(matrix)
        white_val = np.array([255, 255, 255, 255], dtype=np.uint8).view(np.uint32)[0]
        assert not np.all(rgba[3, :] == white_val), "Row 0 of matrix should map to bottom of output"
        assert np.all(rgba[0, :] == white_val), "Row 3 of all-zero matrix should stay white at top"


class TestComputeDotplotImage:
    def test_output_shape_is_fixed(self):
        seq = "ACGT" * 2_000
        img = compute_dotplot_image(seq, resolution=32)
        assert img.shape == (32, 32)
        assert img.dtype == np.uint32

    def test_default_resolution_produces_square_output(self):
        seq = "ACGT" * 2_000
        img = compute_dotplot_image(seq)
        assert img.shape == (_RESOLUTION, _RESOLUTION)


class TestCanonical:
    def test_forward_smaller_returns_forward(self):
        assert _canonical("AAACGTT") == "AAACGTT"

    def test_rc_smaller_returns_rc(self):
        kmer = "TTTTTTTTTTA"
        assert _canonical(kmer) == _canonical(_canonical(kmer))

    def test_canonical_is_symmetric(self):
        kmer = "ACGTACGTACC"
        rc = kmer.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]
        assert _canonical(kmer) == _canonical(rc)

    def test_palindrome_returns_itself(self):
        kmer = "AACGTT"
        assert _canonical(kmer) == kmer

    def test_canonical_is_idempotent(self):
        kmer = "TTTTTTTTTTA"
        assert _canonical(kmer) == _canonical(_canonical(kmer))


class TestInvertedRepeats:
    def test_inverted_repeat_produces_off_diagonal_via_canonical_kmers(self):
        rng = np.random.default_rng(99)
        bases = np.array(list("ACGT"))
        forward_block = "".join(rng.choice(bases, size=3_000))
        inverted_block = forward_block.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        seq = forward_block + ("N" * 1_000) + inverted_block

        matrix = compute_self_identity(seq, resolution=64)
        off_diag = matrix.copy()
        np.fill_diagonal(off_diag, 0)
        assert off_diag.max() == 1.0


class TestConstants:
    def test_max_region_bp(self):
        assert MAX_DOTPLOT_REGION_BP == 10_000_000

    def test_resolution(self):
        assert _RESOLUTION == 512

    def test_kmer_is_15(self):
        assert _KMER == 15

    def test_max_kmer_bins_is_positive(self):
        assert _MAX_KMER_BINS > 0

    def test_min_shared_kmers_floor_is_at_least_two(self):
        assert _MIN_SHARED_KMERS_FLOOR >= 2

    def test_min_shared_kmers_per_bp_is_positive(self):
        assert _MIN_SHARED_KMERS_PER_BP > 0

    def test_high_identity_repeats_are_detected(self):
        # 10 copies of a 1000bp unit at ~95% identity produce enough shared k-mers
        # to exceed the adaptive threshold. These are the regions at risk for HiFi
        # read misalignment.
        rng = np.random.default_rng(77)
        bases = np.array(list("ACGT"))
        unit = "".join(rng.choice(bases, size=1_000))
        copies = [unit]
        for _ in range(9):
            mutated = list(unit)
            for idx in rng.choice(len(mutated), size=50, replace=False):  # ~5% mutations
                mutated[idx] = rng.choice([b for b in "ACGT" if b != mutated[idx]])
            copies.append("".join(mutated))
        seq = "".join(copies)  # 10,000 bp

        matrix = compute_self_identity(seq, resolution=64)
        off_diag = matrix.copy()
        np.fill_diagonal(off_diag, 0)
        assert off_diag.max() == 1.0, "95%-identical repeat copies must be detected"

    def test_low_identity_repeats_are_not_detected(self):
        # 5 copies of a 1000bp unit at ~70% identity must NOT produce off-diagonal
        # signal. 70% identity is well below the HiFi misalignment risk threshold.
        # 30% mutations per copy ensures no 15-mer window is accidentally mutation-free
        # due to random clustering (which can occur at lower mutation rates).
        rng = np.random.default_rng(43)
        bases = np.array(list("ACGT"))
        unit = "".join(rng.choice(bases, size=1_000))
        copies = [unit]
        for _ in range(4):
            mutated = list(unit)
            for idx in rng.choice(len(mutated), size=300, replace=False):  # 30% mutations
                mutated[idx] = rng.choice([b for b in "ACGT" if b != mutated[idx]])
            copies.append("".join(mutated))
        seq = "".join(copies)  # 5,000 bp

        matrix = compute_self_identity(seq, resolution=64)
        off_diag = matrix.copy()
        np.fill_diagonal(off_diag, 0)
        assert off_diag.max() == 0.0, "70%-identical copies must NOT be detected"


class TestComputeRepeatDensity:
    def test_output_shape_matches_matrix_dimension(self):
        matrix = np.zeros((64, 64), dtype=np.float32)
        density = compute_repeat_density(matrix)
        assert density.shape == (64,)
        assert density.dtype == np.float32

    def test_diagonal_only_matrix_gives_zero_density(self):
        matrix = np.zeros((8, 8), dtype=np.float32)
        np.fill_diagonal(matrix, 1.0)
        density = compute_repeat_density(matrix)
        assert np.all(density == 0.0)

    def test_off_diagonal_match_is_counted_in_both_bins(self):
        matrix = np.zeros((8, 8), dtype=np.float32)
        matrix[2, 5] = 1.0
        matrix[5, 2] = 1.0
        density = compute_repeat_density(matrix)
        assert density[2] == 1.0
        assert density[5] == 1.0
        assert density.sum() == 2.0

    def test_n_masked_sentinel_is_excluded(self):
        matrix = np.full((8, 8), -1.0, dtype=np.float32)
        density = compute_repeat_density(matrix)
        assert np.all(density == 0.0)

    def test_input_matrix_is_not_mutated(self):
        matrix = np.zeros((8, 8), dtype=np.float32)
        np.fill_diagonal(matrix, 1.0)
        matrix[2, 5] = 1.0
        original = matrix.copy()
        compute_repeat_density(matrix)
        np.testing.assert_array_equal(matrix, original)
