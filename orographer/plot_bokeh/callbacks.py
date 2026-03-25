"""JavaScript callbacks and HTML/export logic for Bokeh plots."""

import gzip
import json
import logging
import os
import re
import sys

from bokeh.io import output_file as bokeh_output_file
from bokeh.models import CustomJS
from bokeh.plotting import save

from .utils import COORD_INPUT_CENTER_OFFSET_PX, load_javascript

logger = logging.getLogger(__name__)


def get_exon_click_callback(source):
    """Get JavaScript callback for exon clicks (shows modal)."""
    return CustomJS(
        args={"source": source},
        code=load_javascript("exon_click_callback.js"),
    )


def get_vcf_variant_click_callback(source):
    """Get JavaScript callback for VCF variant clicks (shows modal)."""
    return CustomJS(
        args={"source": source},
        code=load_javascript("vcf_variant_click_callback.js"),
    )


def get_arrow_tap_callback(arrow_source):
    """Get JavaScript callback for arrow tap selection (uses manual hit detection)."""
    return CustomJS(
        args={"source": arrow_source},
        code=load_javascript("arrow_tap_callback.js"),
    )


def get_arrow_tap_callback_multi_region(arrow_source, all_arrow_sources):
    """Get JS callback for arrow tap with cross-region highlighting."""
    return CustomJS(
        args={"source": arrow_source, "all_sources": all_arrow_sources},
        code=load_javascript("arrow_tap_multi_region_callback.js"),
    )


def get_number_click_callback(source):
    """Get JavaScript callback for number label clicks (shows modal)."""
    return CustomJS(
        args={"source": source},
        code=load_javascript("number_click_callback.js"),
    )


def get_modal_html():
    """Get HTML string for the details modal (used for alignments and variants)."""
    modal_css = (
        "display:none;position:fixed;z-index:1000;left:0;top:0;width:100%;"
        "height:100%;overflow:auto;background-color:rgba(0,0,0,0.4);"
        "align-items:center;justify-content:center;"
    )
    wrapper_css = "display:flex;align-items:center;justify-content:center;max-height:90vh;"
    dialog_css = "background:#fefefe;padding:20px;border:1px solid #888;width:400px;"
    close_css = "color:#aaa;float:right;font-size:28px;font-weight:bold;cursor:pointer;"
    parts = [
        "\n    <!-- Modal for details (alignments and variants) -->\n",
        '    <div id="alignmentModal" style="',
        modal_css,
        '">\n        <div id="alignmentModalWrapper" style="',
        wrapper_css,
        '">\n            <div id="alignmentModalDialog" style="',
        dialog_css,
        '">\n                <span id="closeModal" style="',
        close_css,
        '">&times;</span>\n                <h2 style="margin-top: 0; color: #333;">'
        'Details</h2>\n                <div id="modalContent" '
        'style="font-size: 14px; line-height: 1.8;"></div>\n            '
        "</div>\n        ",
        "</div>\n    </div>\n    \n    <script>\n    ",
        load_javascript("modal.js"),
        "\n    </script>\n    ",
    ]
    return "".join(parts)


def extract_bokeh_json(html_string):
    """
    Extract docs_json and render_items from Bokeh HTML.
    Bokeh 3.x stores docs_json in a script tag with an ID, and render_items inline.
    Returns (docs_json, render_items) as Python dicts, or (None, None) if not found.
    """
    script_id_pattern = r"document\.getElementById\('([^']+)'\)\.textContent"
    script_id_match = re.search(script_id_pattern, html_string)
    if not script_id_match:
        return None, None
    script_id = script_id_match.group(1)

    script_tag_pattern = rf'<script[^>]*id="{re.escape(script_id)}"[^>]*>(.*?)</script>'
    script_match = re.search(script_tag_pattern, html_string, re.DOTALL)
    if not script_match:
        return None, None
    try:
        docs_json = json.loads(script_match.group(1))
    except json.JSONDecodeError:
        return None, None

    render_items_pattern = r"render_items\s*=\s*(\[.*?\]);"
    render_items_match = re.search(render_items_pattern, html_string, re.DOTALL)
    if not render_items_match:
        return None, None
    try:
        render_items = json.loads(render_items_match.group(1))
        return docs_json, render_items
    except json.JSONDecodeError:
        return None, None


def save_plot_with_modal(layout, output_file, prefix):
    """Save the plot layout to HTML with external JSON and inject modal HTML."""
    file_dir = os.path.dirname(output_file)
    if file_dir:
        os.makedirs(file_dir, exist_ok=True)

    title = f"Orographer Plot: {prefix}" if prefix else "Orographer Plot"

    bokeh_output_file(output_file, title=title)
    save(layout, filename=output_file)

    with open(output_file) as infile:
        html_string = infile.read()

    docs_json, render_items = extract_bokeh_json(html_string)
    if docs_json is None or render_items is None:
        modal_html = get_modal_html()
        html_string = html_string.replace("</body>", modal_html + "</body>")
        with open(output_file, "w") as outfile:
            outfile.write(html_string)
        return

    json_file = os.path.splitext(output_file)[0] + ".json.gz"
    json_data = {"docs_json": docs_json, "render_items": render_items}
    with gzip.open(json_file, "wt", encoding="utf-8") as f:
        json.dump(json_data, f)

    json_filename = os.path.basename(json_file)
    # Prevent stale JSON cache from mismatching newly generated HTML root IDs.
    json_url = f"{json_filename}?v={int(os.path.getmtime(json_file))}"
    script_id_pattern = r"const docs_json\s*=\s*document\.getElementById\('([^']+)'\)\.textContent;"
    script_id_match = re.search(script_id_pattern, html_string)
    if not script_id_match:
        logger.error("Error: Could not find script tag pattern in HTML. Expected Bokeh 3.x.")
        sys.exit(1)
    script_id = script_id_match.group(1)

    script_tag_pattern = rf'<script[^>]*id="{re.escape(script_id)}"[^>]*>.*?</script>'
    html_string = re.sub(script_tag_pattern, "", html_string, flags=re.DOTALL)

    replacement_fetch = load_javascript(
        "embed_replace.js",
        replacements={
            "__JSON_FILENAME__": json_url,
            "__OUTPUT_BASENAME__": os.path.basename(output_file),
            "__CENTER_OFFSET_WIDE__": str(COORD_INPUT_CENTER_OFFSET_PX),
        },
    )

    full_pattern = (
        r"const docs_json\s*=\s*document\.getElementById\("
        + re.escape(f"'{script_id}'")
        + r"\)\.textContent;\s*"
        + r"const render_items\s*=\s*\[.*?\];\s*"
        + r"root\.Bokeh\.embed\.embed_items\(docs_json,\s*render_items\);"
    )
    full_match = re.search(full_pattern, html_string, re.DOTALL)
    if not full_match:
        logger.error("Error: Could not find Bokeh code pattern in HTML after removing tag.")
        sys.exit(1)

    new_html = (
        html_string[: full_match.start()]
        + replacement_fetch
        + "\n"
        + html_string[full_match.end() :]
    )
    modal_html = get_modal_html()
    new_html = new_html.replace("</body>", modal_html + "</body>")
    with open(output_file, "w") as outfile:
        outfile.write(new_html)
