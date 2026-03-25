"""Unit tests for plot external JSON schema and genomic range extraction."""

import gzip
import json
from pathlib import Path

import pytest

from tests.helpers.plot_json_schema import (
    assert_genomic_ranges_in_docs_json,
    assert_orig_region_bounds_in_docs_json,
    assert_sample_labels_in_docs_json,
    collect_genomic_like_start_end_pairs,
    collect_orig_region_bounds_from_docs_json,
    validate_plot_json_schema,
)


def test_validate_plot_json_schema_on_expected_fixture():
    p = Path(
        "tests/data/expected_outputs/paraviewer/HG002.paraphase.small/hba/"
        "chr16_171801_175500_bokeh.json.gz"
    )
    with gzip.open(p, "rt", encoding="utf-8") as f:
        data = json.load(f)
    validate_plot_json_schema(data)
    assert_orig_region_bounds_in_docs_json(data["docs_json"], [(171801, 175500)])
    assert_genomic_ranges_in_docs_json(data["docs_json"], [(171801, 175500)])


def test_validate_plot_json_schema_rejects_wrong_top_level():
    with pytest.raises(AssertionError, match="Expected keys"):
        validate_plot_json_schema({"docs_json": {}, "foo": []})


def test_validate_plot_json_schema_rejects_empty_docs_json():
    with pytest.raises(AssertionError, match="non-empty"):
        validate_plot_json_schema(
            {
                "docs_json": {},
                "render_items": [
                    {"docid": "x", "roots": [], "root_ids": []},
                ],
            }
        )


def test_validate_plot_json_schema_rejects_empty_render_items():
    with pytest.raises(AssertionError, match="render_items"):
        validate_plot_json_schema(
            {
                "docs_json": {"uuid": {"roots": [{"type": "object"}]}},
                "render_items": [],
            }
        )


def test_collect_genomic_like_start_end_pairs_filters_small_ranges():
    pairs = collect_genomic_like_start_end_pairs(
        {"a": {"start": 10, "end": 20}, "b": {"start": 171801, "end": 175500}}
    )
    assert (171801, 175500) in pairs
    assert (10, 20) not in pairs


def test_collect_genomic_always_include_expected_small_range():
    pairs = collect_genomic_like_start_end_pairs(
        {"x": {"start": 100, "end": 500}},
        always_include_ranges={(100, 500)},
    )
    assert (100, 500) in pairs


def test_collect_orig_region_bounds_from_nested_lists():
    docs = {
        "roots": [
            [
                ["orig_start", 1_000_000],
                ["orig_end", 1_000_100],
                ["other", 1],
            ],
            [["orig_start", 2_000_000], ["orig_end", 2_000_050]],
        ]
    }
    assert collect_orig_region_bounds_from_docs_json(docs) == [
        (1_000_000, 1_000_100),
        (2_000_000, 2_000_050),
    ]
    assert_orig_region_bounds_in_docs_json(docs, [(1_000_000, 1_000_100), (2_000_000, 2_000_050)])


def test_assert_sample_labels_in_docs_json_synthetic():
    docs = {
        "a": {
            "type": "object",
            "name": "Div",
            "attributes": {"text": "TrackA"},
        },
        "b": {
            "type": "object",
            "name": "ColumnDataSource",
            "attributes": {
                "data": {
                    "x": [1],
                    "sample_label": ["TrackA"],
                }
            },
        },
        "c": {
            "name": "ColumnDataSource",
            "attributes": {"data": {"sample_label": ["", "TrackB", ""]}},
        },
        "d": {
            "name": "Div",
            "attributes": {"text": "TrackB"},
        },
    }
    assert_sample_labels_in_docs_json(docs, ["TrackA", "TrackB"])


def test_assert_sample_labels_missing_div_raises():
    docs = {
        "b": {
            "name": "ColumnDataSource",
            "attributes": {"data": {"sample_label": ["OnlyCDS"]}},
        }
    }
    with pytest.raises(AssertionError, match="Div"):
        assert_sample_labels_in_docs_json(docs, ["OnlyCDS"])


def test_assert_sample_labels_bokeh_map_encoded_cds():
    """Bokeh 3 embeds ColumnDataSource.data as {type: map, entries: [[key, val], ...]}."""
    docs = {
        "d1": {"name": "Div", "attributes": {"text": "L1"}},
        "cds": {
            "name": "ColumnDataSource",
            "attributes": {
                "data": {
                    "type": "map",
                    "entries": [
                        ["x", [1]],
                        ["sample_label", ["L1"]],
                    ],
                }
            },
        },
    }
    assert_sample_labels_in_docs_json(docs, ["L1"])
