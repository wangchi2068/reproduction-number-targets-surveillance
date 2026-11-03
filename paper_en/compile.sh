#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
PYTHON_BIN="${PYTHON:-python}"
"$PYTHON_BIN" "$SCRIPT_DIR/../experiment/code/build_papers.py" --lang en
