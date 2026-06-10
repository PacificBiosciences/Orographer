#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

cd "${REPO_ROOT}"

PYTHON_BIN="${PYTHON_BIN:-.venv_e2e/bin/python}"
if [[ ! -x "${PYTHON_BIN}" ]]; then
    PYTHON_BIN="${PYTHON:-python}"
fi

# Python unit tests cover pure Python logic without browser or Node dependencies.
echo "==> Running Python unit tests"
"${PYTHON_BIN}" -m pytest tests/unit "$@"

# JavaScript unit tests exercise standalone callback helpers with mocked Bokeh objects.
echo "==> Running JavaScript unit tests"
NODE_BIN="${NODE_BIN:-node}"
if command -v "${NODE_BIN}" >/dev/null 2>&1; then
    "${NODE_BIN}" tests/js/run_js_unit_tests.mjs
elif command -v pixi >/dev/null 2>&1; then
    pixi run node tests/js/run_js_unit_tests.mjs
else
    echo "Skipping JS unit tests: node is not available on PATH." >&2
fi

# Browser tests load generated Bokeh HTML in Playwright-backed browsers.
echo "==> Running Playwright browser tests"
if "${PYTHON_BIN}" -c "import playwright.sync_api" >/dev/null 2>&1; then
    "${PYTHON_BIN}" -m pytest tests/browser "$@"
elif command -v pixi >/dev/null 2>&1; then
    pixi run python -m pytest tests/browser "$@"
else
    "${PYTHON_BIN}" -m pytest tests/browser "$@"
fi
