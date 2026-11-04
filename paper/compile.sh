#!/bin/bash
set -e

# This script is meant to be run from within the paper output directory
# where main.tex, bibliography.bib, and sections/ already exist.
#
# Engine selection (automatic):
#   fontspec or ctexart -> xelatex (required for exact manuscript fonts)
#   legacy templates   -> pdflatex
#
# Guard: CJK prose under a non-ctexart class is a setup error — fail fast
# with a clear message instead of producing a broken PDF.

if grep -Eq '^[^%]*\\documentclass[^{]*\{ctexart\}|^[^%]*\\usepackage([[:space:]]*\[[^]]*\])?[[:space:]]*\{fontspec\}' main.tex; then
    ENGINE=xelatex
else
    ENGINE=pdflatex
    if python3 - <<'EOF'
import sys
from pathlib import Path

for f in [Path("main.tex"), *sorted(Path("sections").rglob("*.tex"))]:
    try:
        text = f.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        continue
    if any("一" <= c <= "鿿" for c in text):
        print(f"CJK text found in {f}", file=sys.stderr)
        sys.exit(0)
sys.exit(1)
EOF
    then
        echo "ERROR: Chinese text found but \\documentclass is not ctexart." >&2
        echo "Chinese papers must switch to ctexart (see the CHINESE PAPER" >&2
        echo "SWITCH comment at the top of main.tex); pdflatex cannot" >&2
        echo "typeset CJK. Aborting before a broken compile." >&2
        exit 1
    fi
fi

echo "Engine: $ENGINE"

if ! python3 -c 'import pypdf, PIL' 2>/dev/null; then
    echo "ERROR: manuscript QA requires pypdf and Pillow." >&2
    echo "Install qa/requirements.txt (source bundle) or the equivalent runtime packages." >&2
    exit 4
fi

# --- Manuscript source QA ----------------------------------------------------
MANUSCRIPT_QA=""
if [ -f "$PWD/qa/manuscript_qa.py" ]; then
    MANUSCRIPT_QA="$PWD/qa/manuscript_qa.py"
fi
PROBE=$PWD
for _ in 1 2 3 4 5 6 7 8; do
    [ -n "$MANUSCRIPT_QA" ] && break
    if [ -f "$PROBE/.claude/skills/_shared/scripts/manuscript_qa.py" ]; then
        MANUSCRIPT_QA="$PROBE/.claude/skills/_shared/scripts/manuscript_qa.py"
        break
    fi
    PROBE=$(dirname "$PROBE")
done

if [ -z "$MANUSCRIPT_QA" ]; then
    echo "MANUSCRIPT-QA: required manuscript_qa.py not found above $PWD" >&2
    exit 4
elif ! python3 "$MANUSCRIPT_QA" --source-root "$PWD" --source-only --json manuscript-source-qa.json; then
    echo "MANUSCRIPT-QA: source defects found; fix them before compiling." >&2
    exit 4
fi

echo "First $ENGINE run..."
$ENGINE -interaction=nonstopmode -halt-on-error main.tex

echo "Running bibtex..."
if [ -f bibliography.bib ]; then
    bibtex main
fi

echo "Second $ENGINE run..."
$ENGINE -interaction=nonstopmode -halt-on-error main.tex

echo "Final $ENGINE run..."
$ENGINE -interaction=nonstopmode -halt-on-error main.tex

echo "Compilation complete. PDF: main.pdf"

# --- Render-QA (mechanical) --------------------------------------------------
# Text-on-text overlap and per-character CJK wrapping are computable from the
# PDF's own geometry, and an agent skimming an A4 page raster has already
# missed both in production. Reported on every compile so the numbers are
# unavoidable; the hard stop itself lives in the quality gate (G5.5a / G5.6a).
FIGQA=""
if [ -f "$PWD/qa/figqa.py" ]; then
    FIGQA="$PWD/qa/figqa.py"
fi
PROBE=$PWD
for _ in 1 2 3 4 5 6 7 8; do
    [ -n "$FIGQA" ] && break
    if [ -f "$PROBE/.claude/skills/paper-writer/scripts/figqa.py" ]; then
        FIGQA="$PROBE/.claude/skills/paper-writer/scripts/figqa.py"
        break
    fi
    PROBE=$(dirname "$PROBE")
done

if [ -z "$FIGQA" ]; then
    echo "RENDER-QA: required figqa.py not found above $PWD" >&2
    exit 3
elif [ ! -s main.pdf ]; then
    echo "RENDER-QA: main.pdf is missing or empty" >&2
    exit 3
elif python3 "$FIGQA" main.pdf --json figqa.json; then
    echo "RENDER-QA: clean (0 defects)"
else
    echo "RENDER-QA: DEFECTS FOUND — listed above, also in figqa.json." >&2
    echo "  Fix the figure/table SOURCE (see 02-*figures*.md), recompile," >&2
    echo "  re-run. This is hard gate G5.5a/G5.6a; do not ship it." >&2
    echo "  main.pdf was written but MUST NOT be delivered as-is." >&2
    exit 3   # distinct from LaTeX failure: the PDF exists, it is just defective
fi

if [ -n "$MANUSCRIPT_QA" ]; then
    if python3 "$MANUSCRIPT_QA" --source-root "$PWD" --pdf main.pdf --json manuscript-qa.json --allow-font-fallback; then
        echo "MANUSCRIPT-QA: clean"
    else
        echo "MANUSCRIPT-QA: PDF typography or layout defects found." >&2
        echo "  Inspect manuscript-qa.json, fix the source, and recompile." >&2
        exit 4
    fi
fi
