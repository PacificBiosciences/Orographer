import argparse
import gzip
import json
import os
import re
import sys
from collections import defaultdict

import pysam
from skimage import io
from skimage.metrics import structural_similarity as ssim
from skimage.transform import resize

# IDs we inject (modal, etc.); do not normalize when comparing Bokeh HTML/JSON
BOKEH_FIXED_IDS = frozenset(
    {
        "alignmentModal",
        "alignmentModalWrapper",
        "alignmentModalDialog",
        "closeModal",
        "modalContent",
    }
)
BOKEH_ROOT_PLACEHOLDER = "bokeh-root"
# Bokeh model IDs are like p1234; we normalize to p0, p1, ... for comparison
BOKEH_MODEL_ID_RE = re.compile(r"^p\d+$")
# Document UUIDs (docs_json keys and render_items[].roots values)
BOKEH_UUID_RE = re.compile(
    r"^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$"
)

MATCHING_SSI = 0.75


class Colors:
    RED = "\033[91m"
    GREEN = "\033[92m"
    RESET = "\033[0m"


def setup():
    parser = argparse.ArgumentParser(
        description="Compare two images to determine if essentially identical",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-r",
        "--reference-file-path",
        help="the path to a reference file: json, image, or test",
    )
    parser.add_argument(
        "-t",
        "--test-file-path",
        help="the path to a test file: json, image, or text",
    )
    args = parser.parse_args()
    return args


def compare_plots(reference_image_path, test_image_path):
    """
    Compare two images and report an error if they're dissimilar
    """
    test_image = io.imread(test_image_path)
    reference_image = io.imread(reference_image_path)
    if test_image.shape != reference_image.shape:
        test_image = resize(
            test_image, reference_image.shape, preserve_range=True, anti_aliasing=True
        ).astype(reference_image.dtype)

    # window size has to be equal to or less than the smallest side of the input images
    win_size = int(min(test_image.shape[:2] + reference_image.shape[:2]))
    # must be odd
    if win_size % 2 == 0:
        win_size -= 1
    similarity_index, _ = ssim(
        test_image, reference_image, full=True, win_size=win_size, channel_axis=-1
    )
    print("SSI:", similarity_index)
    if similarity_index < MATCHING_SSI:
        msg = (
            f"{Colors.RED}Test image {test_image_path} does not match reference "
            f"image {reference_image_path}: SSI={similarity_index}{Colors.RESET}"
        )
        print(msg, file=sys.stderr)
        sys.exit(-1)
    else:
        msg = (
            f"{Colors.GREEN}Test image {test_image_path} matches reference "
            f"image {reference_image_path}: SSI={similarity_index}{Colors.RESET}"
        )
        print(msg, file=sys.stdout)


def get_bam_alignment_counts(bamfile):
    """
    Given a bam file, return a dict of the read names mapped to
    the number of alignments for each from the bam
    """
    counts = defaultdict(int)
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            counts[read.query_name] += 1
    return counts


def compare_bams(reference_bam_path, test_bam_path):
    """
    Compare two BAM files and report an error if the reads present are not the same
    """
    test_read_counts = get_bam_alignment_counts(test_bam_path)
    reference_read_counts = get_bam_alignment_counts(reference_bam_path)
    error_message = None
    for read_name in test_read_counts:
        if read_name not in reference_read_counts:
            error_message = (
                f"{Colors.RED}Test BAM {test_bam_path} includes read {read_name} "
                f"not in reference BAM {reference_bam_path}{Colors.RESET}"
            )
            break
        elif test_read_counts[read_name] != reference_read_counts[read_name]:
            error_message = (
                f"{Colors.RED}Test BAM {test_bam_path} read {read_name} count "
                f"{test_read_counts[read_name]} does not match reference BAM "
                f"{reference_bam_path} count {reference_read_counts[read_name]} "
                f"{Colors.RESET}"
            )
            break
    if not error_message:
        for read_name in reference_read_counts:
            if read_name not in test_read_counts:
                error_message = (
                    f"{Colors.RED}Test BAM {test_bam_path} does not include "
                    f"read {read_name} in reference BAM "
                    f"{reference_bam_path}{Colors.RESET}"
                )
                break
    if error_message:
        print(
            error_message,
            file=sys.stderr,
        )
        sys.exit(-1)
    else:
        msg = (
            f"{Colors.GREEN}Test BAM {test_bam_path} matches reference BAM "
            f"{reference_bam_path}{Colors.RESET}"
        )
        print(msg, file=sys.stdout)


def normalize_bokeh_html(html):
    """
    Replace Bokeh-generated div/script IDs and cache-bust query params so HTML
    from different runs compares equal.
    """

    # Normalize random Bokeh embed root/script IDs
    def repl(m):
        val = m.group(1)
        if val in BOKEH_FIXED_IDS:
            return m.group(0)
        return f'id="{BOKEH_ROOT_PLACEHOLDER}"'

    html = re.sub(r'id="([a-zA-Z0-9_-]{12,})"', repl, html)
    # Normalize cache-bust ?v=<mtime> on JSON URL (differs every run)
    html = re.sub(r"\?v=\d+", "?v=0", html)
    return html


def _collect_bokeh_model_ids(obj, ids=None):
    """Collect all Bokeh model ID strings (p1234) in a structure."""
    if ids is None:
        ids = set()
    if isinstance(obj, dict):
        for k, v in obj.items():
            if isinstance(k, str) and BOKEH_MODEL_ID_RE.match(k):
                ids.add(k)
            _collect_bokeh_model_ids(k, ids)
            _collect_bokeh_model_ids(v, ids)
    elif isinstance(obj, list):
        for x in obj:
            _collect_bokeh_model_ids(x, ids)
    elif isinstance(obj, str) and BOKEH_MODEL_ID_RE.match(obj):
        ids.add(obj)
    return ids


def _replace_bokeh_model_ids(obj, id_map):
    """Replace Bokeh model IDs with canonical p0, p1, ... using id_map."""
    if isinstance(obj, dict):
        return {
            (
                id_map[k]
                if isinstance(k, str) and BOKEH_MODEL_ID_RE.match(k)
                else _replace_bokeh_model_ids(k, id_map)
            ): _replace_bokeh_model_ids(v, id_map)
            for k, v in obj.items()
        }
    if isinstance(obj, list):
        return [_replace_bokeh_model_ids(x, id_map) for x in obj]
    if isinstance(obj, str) and BOKEH_MODEL_ID_RE.match(obj):
        return id_map.get(obj, obj)
    return obj


def _normalize_bokeh_docs_json(docs_json):
    """Normalize docs_json so Bokeh model IDs (p1234, roots keys, etc.) are canonical."""
    ids = _collect_bokeh_model_ids(docs_json)
    id_map = {old: f"p{i}" for i, old in enumerate(sorted(ids))}
    return _replace_bokeh_model_ids(docs_json, id_map)


def _collect_uuid_strings(obj, out=None):
    """Collect all strings that match UUID format (dict keys and values) in obj."""
    if out is None:
        out = set()
    if isinstance(obj, dict):
        for k, v in obj.items():
            if isinstance(k, str) and BOKEH_UUID_RE.match(k):
                out.add(k)
            _collect_uuid_strings(k, out)
            _collect_uuid_strings(v, out)
    elif isinstance(obj, list):
        for x in obj:
            _collect_uuid_strings(x, out)
    elif isinstance(obj, str) and BOKEH_UUID_RE.match(obj):
        out.add(obj)
    return out


def _replace_doc_uuids(obj, uuid_map):
    """Replace document UUID strings (keys and values) with canonical doc0, doc1, ..."""
    if isinstance(obj, dict):
        return {
            (uuid_map[k] if isinstance(k, str) and k in uuid_map else k): _replace_doc_uuids(
                v, uuid_map
            )
            for k, v in obj.items()
        }
    if isinstance(obj, list):
        return [_replace_doc_uuids(x, uuid_map) for x in obj]
    if isinstance(obj, str) and obj in uuid_map:
        return uuid_map[obj]
    return obj


def normalize_bokeh_json(obj):
    """
    Replace Bokeh embed keys and auto-generated IDs so JSON from different
    runs compares equal: root_id, elementid, docid, docs_json keys (UUIDs),
    docs_json model IDs (p1234). Caller should then replace remaining UUID
    strings (e.g. in render_items[].roots) via _collect_uuid_strings +
    _replace_doc_uuids.
    """
    if isinstance(obj, dict):
        result = {}
        for k, v in obj.items():
            if k in ("root_id", "elementid", "docid"):
                result[k] = BOKEH_ROOT_PLACEHOLDER
            elif k == "docs_json":
                # Top-level keys are document UUIDs; normalize to doc0, doc1, ...
                if not isinstance(v, dict):
                    result[k] = normalize_bokeh_json(v)
                else:
                    sorted_keys = sorted(v.keys(), key=str.lower)
                    new_docs = {}
                    doc_i = 0
                    for key in sorted_keys:
                        if BOKEH_UUID_RE.match(key):
                            new_docs[f"doc{doc_i}"] = _normalize_bokeh_docs_json(v[key])
                            doc_i += 1
                        else:
                            new_docs[key] = normalize_bokeh_json(v[key])
                    result[k] = new_docs
            else:
                result[k] = normalize_bokeh_json(v)
        return result
    if isinstance(obj, list):
        return [normalize_bokeh_json(x) for x in obj]
    return obj


def _json_diff_paths(ref, test, path=()):
    """
    Recursively compare ref and test; yield (path_str, ref_val, test_val) for
    each leaf where they differ. path is a tuple of keys/indices.
    """

    def path_str(p):
        parts = []
        for seg in p:
            if isinstance(seg, int):
                parts.append(f"[{seg}]")
            else:
                parts.append(str(seg) if not parts else f".{seg}")
        return "".join(parts) if parts else "(root)"

    if type(ref) is not type(test):
        yield (path_str(path), ref, test)
        return
    if isinstance(ref, dict):
        all_keys = set(ref) | set(test)
        for k in sorted(all_keys):
            sub = (*path, k)
            if k not in ref:
                yield (path_str(sub), "(missing in reference)", test[k])
            elif k not in test:
                yield (path_str(sub), ref[k], "(missing in test)")
            else:
                yield from _json_diff_paths(ref[k], test[k], sub)
    elif isinstance(ref, list):
        if len(ref) != len(test):
            yield (path_str(path), f"len={len(ref)}", f"len={len(test)}")
        for i in range(min(len(ref), len(test))):
            yield from _json_diff_paths(ref[i], test[i], (*path, i))
    else:
        if ref != test:
            yield (path_str(path), ref, test)


def compare_html(reference_html_path, test_html_path):
    """
    Compare two HTML files (e.g. Bokeh output), normalizing random Bokeh IDs.
    """
    with open(reference_html_path, encoding="utf-8") as ref_fh:
        ref_text = normalize_bokeh_html(ref_fh.read().strip())
    with open(test_html_path, encoding="utf-8") as test_fh:
        test_text = normalize_bokeh_html(test_fh.read().strip())

    if ref_text != test_text:
        msg = (
            f"{Colors.RED}Test HTML {test_html_path} does not match reference "
            f"HTML {reference_html_path}{Colors.RESET}"
        )
        print(msg, file=sys.stderr)
        sys.exit(-1)
    msg = (
        f"{Colors.GREEN}Test HTML {test_html_path} matches reference "
        f"HTML {reference_html_path}{Colors.RESET}"
    )
    print(msg, file=sys.stdout)


def compare_jsons(reference_json_path, test_json_path, compressed=False):
    """
    Compare two JSON files (optionally gzip-compressed); error if they don't match.
    """
    with (
        (
            gzip.open(reference_json_path, "rt", encoding="utf-8")
            if compressed
            else open(reference_json_path, encoding="utf-8")
        ) as ref_fh,
        (
            gzip.open(test_json_path, "rt", encoding="utf-8")
            if compressed
            else open(test_json_path, encoding="utf-8")
        ) as test_fh,
    ):
        ref_json = json.load(ref_fh)
        try:
            test_json = json.load(test_fh)
        except json.decoder.JSONDecodeError:
            print(
                f"{Colors.RED}Test json {test_json_path} is empty or malformed {Colors.RESET}",
                file=sys.stderr,
            )
            sys.exit(-1)
    ref_norm = normalize_bokeh_json(ref_json)
    ref_uuids = _collect_uuid_strings(ref_norm)
    ref_uuid_map = {u: f"doc{i}" for i, u in enumerate(sorted(ref_uuids, key=str.lower))}
    if ref_uuid_map:
        ref_norm = _replace_doc_uuids(ref_norm, ref_uuid_map)

    test_norm = normalize_bokeh_json(test_json)
    test_uuids = _collect_uuid_strings(test_norm)
    test_uuid_map = {u: f"doc{i}" for i, u in enumerate(sorted(test_uuids, key=str.lower))}
    if test_uuid_map:
        test_norm = _replace_doc_uuids(test_norm, test_uuid_map)
    ref_canonical = json.dumps(ref_norm, sort_keys=True)
    test_canonical = json.dumps(test_norm, sort_keys=True)
    if ref_canonical != test_canonical:
        with open("tmp1", "w") as fh:
            fh.write(json.dumps(ref_norm, sort_keys=True, indent=2))
        with open("tmp2", "w") as fh:
            fh.write(json.dumps(test_norm, sort_keys=True, indent=2))
        msg = (
            f"{Colors.RED}Test json {test_json_path} does not match reference "
            f"json {reference_json_path}{Colors.RESET}"
        )
        print(msg, file=sys.stderr)

        def _val_summary(v):
            if isinstance(v, dict):
                return f"<dict len={len(v)}>"
            if isinstance(v, list):
                return f"<list len={len(v)}>"
            s = repr(v)
            return s if len(s) <= 60 else s[:57] + "..."

        diffs = list(_json_diff_paths(ref_norm, test_norm))
        print("Mismatch paths (reference vs test):", file=sys.stderr)
        for path, rv, tv in diffs:
            print(f"  {path}: {_val_summary(rv)} vs {_val_summary(tv)}", file=sys.stderr)
        sys.exit(-1)
    else:
        msg = (
            f"{Colors.GREEN}Test json {test_json_path} matches reference "
            f"json {reference_json_path}{Colors.RESET}"
        )
        print(msg, file=sys.stdout)


def compare_text_files(reference_text_path, test_text_path):
    """
    Compare two text files from svtopo and error if they don't match
    """
    with open(reference_text_path, encoding="utf-8") as ref_fh:
        ref_text = ref_fh.read().strip()
    with open(test_text_path, encoding="utf-8") as test_fh:
        test_text = test_fh.read().strip()

    if ref_text != test_text:
        msg = (
            f"{Colors.RED}Test text file {test_text_path} does not match "
            f"reference text file {reference_text_path}{Colors.RESET}"
        )
        print(msg, file=sys.stderr)
        sys.exit(-1)
    else:
        msg = (
            f"{Colors.GREEN}Test text file {test_text_path} matches reference "
            f"text file {reference_text_path}{Colors.RESET}"
        )
        print(msg, file=sys.stdout)


def check_exists(ref_file_path, test_file_path):
    ref_exists = os.path.isfile(ref_file_path)
    test_exists = os.path.isfile(test_file_path)
    error_msg = None
    if ref_exists and not test_exists:
        error_msg = (
            f"{Colors.RED}Test file {test_file_path} does not exist, "
            f"reference file {ref_file_path} does{Colors.RESET}"
        )
    elif not ref_exists and test_exists:
        error_msg = (
            f"{Colors.RED}Test file {test_file_path} exists but has no "
            f"matching file at {ref_file_path}{Colors.RESET}"
        )
    if error_msg:
        print(
            error_msg,
            file=sys.stdout,
        )
        sys.exit(-1)


def main():
    args = setup()
    ref_file_path = args.reference_file_path
    test_file_path = args.test_file_path

    check_exists(ref_file_path, test_file_path)

    if ref_file_path.endswith(".json.gz") or ref_file_path.endswith("json"):
        compare_jsons(
            reference_json_path=ref_file_path,
            test_json_path=test_file_path,
            compressed=ref_file_path.endswith(".json.gz"),
        )
    elif ref_file_path.endswith("png"):
        compare_plots(reference_image_path=ref_file_path, test_image_path=test_file_path)
    elif ref_file_path.endswith("bam"):
        compare_bams(reference_bam_path=ref_file_path, test_bam_path=test_file_path)
    elif ref_file_path.endswith(".html"):
        compare_html(reference_html_path=ref_file_path, test_html_path=test_file_path)
    else:
        compare_text_files(reference_text_path=ref_file_path, test_text_path=test_file_path)


main()
