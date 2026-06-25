#!/usr/bin/env bash
set -euo pipefail

pixi run python -m pytest tests/browser "$@"
