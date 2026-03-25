"""
Deeper validation of orographer plot external JSON (docs_json + render_items).

Bokeh embed JSON stores x-axis DataRange bounds as objects with numeric ``start``/``end``.
We collect pairs that look like genomic coordinates (large integers) to avoid matching
plot y-ranges or other small numeric pairs.
"""

from __future__ import annotations

from contextlib import suppress
from typing import Any

# Ignore small numeric ranges (y axes, padding, etc.)
_MIN_GENOMIC_LIKE = 50_000


def validate_plot_json_schema(data: dict[str, Any]) -> None:
    """Assert top-level shape matches orographer's gzip JSON contract."""
    if set(data.keys()) != {"docs_json", "render_items"}:
        raise AssertionError(f"Expected keys {{docs_json, render_items}}, got {set(data.keys())}")

    docs_json = data["docs_json"]
    if not isinstance(docs_json, dict) or not docs_json:
        raise AssertionError("docs_json must be a non-empty dict")

    for doc_key, doc in docs_json.items():
        if not isinstance(doc_key, str) or not doc_key:
            raise AssertionError(f"Invalid docs_json key: {doc_key!r}")
        if not isinstance(doc, dict):
            raise AssertionError(f"docs_json[{doc_key!r}] must be a dict")
        if "roots" not in doc:
            raise AssertionError(f"docs_json[{doc_key!r}] missing 'roots'")
        roots = doc["roots"]
        if not isinstance(roots, list) or not roots:
            raise AssertionError(f"docs_json[{doc_key!r}].roots must be non-empty list")

    items = data["render_items"]
    if not isinstance(items, list) or not items:
        raise AssertionError("render_items must be a non-empty list")

    for i, item in enumerate(items):
        if not isinstance(item, dict):
            raise AssertionError(f"render_items[{i}] must be a dict")
        for k in ("docid", "roots", "root_ids"):
            if k not in item:
                raise AssertionError(f"render_items[{i}] missing {k!r}")


def collect_genomic_like_start_end_pairs(
    obj: Any,
    *,
    always_include_ranges: set[tuple[int, int]] | None = None,
) -> set[tuple[int, int]]:
    """
    Recursively collect (start, end) pairs that look like genomic x-ranges.

    Pairs listed in ``always_include_ranges`` are included even when below the
    usual genomic magnitude floor (for testing small loci without weakening
    filtering for unrelated small numeric ranges).
    """
    always = always_include_ranges or set()
    out: set[tuple[int, int]] = set()

    def walk(o: Any) -> None:
        if isinstance(o, dict):
            if "start" in o and "end" in o:
                s, e = o["start"], o["end"]
                if (
                    isinstance(s, (int, float))
                    and isinstance(e, (int, float))
                    and (
                        not (isinstance(s, float) or isinstance(e, float))
                        or (s == int(s) and e == int(e))
                    )
                ):
                    si, ei = round(s), round(e)
                    if (
                        ei >= si
                        and (ei - si) <= 50_000_000
                        and (
                            (si, ei) in always
                            or (si >= _MIN_GENOMIC_LIKE and ei >= _MIN_GENOMIC_LIKE)
                        )
                    ):
                        out.add((si, ei))
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)

    walk(obj)
    return out


def collect_orig_region_bounds_from_docs_json(docs_json: dict[str, Any]) -> list[tuple[int, int]]:
    """
    Collect (orig_start, orig_end) from Bokeh-embedded CustomJS args.

    Orographer passes these for coord / zoom callbacks; they match the plotted
    1-based inclusive region bounds regardless of axis padding heuristics.
    """
    found: list[tuple[int, int]] = []

    def walk(o: Any) -> None:
        if isinstance(o, list):
            for i in range(len(o) - 1):
                a, b = o[i], o[i + 1]
                if (
                    isinstance(a, list)
                    and len(a) == 2
                    and a[0] == "orig_start"
                    and isinstance(b, list)
                    and len(b) == 2
                    and b[0] == "orig_end"
                ):
                    with suppress(TypeError, ValueError):
                        found.append((int(a[1]), int(b[1])))
            for v in o:
                walk(v)
        elif isinstance(o, dict):
            for v in o.values():
                walk(v)

    walk(docs_json)
    return found


def assert_orig_region_bounds_in_docs_json(
    docs_json: dict[str, Any],
    expected_ranges: list[tuple[int, int]],
    *,
    min_occurrences_per_range: int = 1,
) -> None:
    """Assert each expected (start, end) appears in CustomJS orig_start/orig_end args."""
    from collections import Counter

    collected = collect_orig_region_bounds_from_docs_json(docs_json)
    counts = Counter(collected)
    for r in expected_ranges:
        if counts[r] < min_occurrences_per_range:
            raise AssertionError(
                f"Expected orig region bounds {r} at least "
                f"{min_occurrences_per_range}x in docs_json (CustomJS args), "
                f"found {counts[r]}. Distinct bounds seen: {sorted(counts.items())}"
            )


def assert_genomic_ranges_in_docs_json(
    docs_json: dict[str, Any],
    expected_ranges: list[tuple[int, int]],
) -> None:
    """Assert each (start, end) appears as a plausible x-range or matches expected exactly."""
    always = set(expected_ranges)
    pairs = collect_genomic_like_start_end_pairs(docs_json, always_include_ranges=always)
    missing = [r for r in expected_ranges if r not in pairs]
    if missing:
        raise AssertionError(
            f"Expected genomic x-ranges {missing} not found in docs_json. "
            f"Sample collected pairs: {sorted(pairs)[:20]}..."
        )


def assert_sample_labels_in_docs_json(docs_json: dict[str, Any], labels: list[str]) -> None:
    """
    Assert each label is wired into the plot JSON in two stable places:

    1. Bokeh ``Div`` with ``attributes.text`` exactly equal to the label (track title).
    2. At least one ``ColumnDataSource.attributes.data['sample_label']`` list entry.
    """
    if not labels:
        raise AssertionError("labels must be non-empty")
    div_ok = dict.fromkeys(labels, False)
    cds_hits: set[str] = set()

    def walk(o: Any) -> None:
        if isinstance(o, dict):
            if o.get("name") == "Div":
                t = (o.get("attributes") or {}).get("text")
                if isinstance(t, str) and t in div_ok:
                    div_ok[t] = True
            if o.get("name") == "ColumnDataSource":
                data = (o.get("attributes") or {}).get("data")
                if isinstance(data, dict):
                    if data.get("type") == "map" and isinstance(data.get("entries"), list):
                        for entry in data["entries"]:
                            if (
                                isinstance(entry, list)
                                and len(entry) >= 2
                                and entry[0] == "sample_label"
                                and isinstance(entry[1], list)
                            ):
                                for x in entry[1]:
                                    if isinstance(x, str) and x:
                                        cds_hits.add(x)
                    else:
                        sl = data.get("sample_label")
                        if isinstance(sl, list):
                            for x in sl:
                                if isinstance(x, str) and x:
                                    cds_hits.add(x)
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)

    walk(docs_json)
    missing_div = [lab for lab, ok in div_ok.items() if not ok]
    if missing_div:
        raise AssertionError(f"Missing Bokeh Div title(s) for sample label(s): {missing_div}")
    missing_cds = [lab for lab in labels if lab not in cds_hits]
    if missing_cds:
        raise AssertionError(
            f"Sample label(s) not present in ColumnDataSource.sample_label data: "
            f"{missing_cds} (found CDS labels: {sorted(cds_hits)})"
        )
