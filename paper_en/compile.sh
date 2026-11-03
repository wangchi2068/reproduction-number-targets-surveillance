#!/bin/bash
set -e

# Build paper_en. Both manuscripts require XeLaTeX (fontspec + unicode-math);
# pdflatex fails by design inside main.tex. Prefer a self-contained Tectonic
# binary when one is available, and fall back to a local XeLaTeX + BibTeX toolchain.
#
# Tectonic bundles the TeX engine but not the system text faces, so the font probe
# below runs before either engine: main.tex requires Times New Roman, with a
# metric-compatible Liberation Serif fallback.

# --- font probe --------------------------------------------------------------
have_font() {
    local pat
    for pat in "$@"; do
        if command -v fc-list >/dev/null 2>&1 && fc-list 2>/dev/null | grep -qi "${pat%\*}"; then
            return 0
        fi
        local d
        for d in /c/Windows/Fonts /mnt/c/Windows/Fonts \
                 "$HOME/AppData/Local/Microsoft/Windows/Fonts" \
                 "$HOME/Library/Fonts" /usr/share/fonts /usr/local/share/fonts; do
            [ -d "$d" ] || continue
            if find "$d" -maxdepth 4 -iname "$pat" -print -quit 2>/dev/null | grep -q .; then
                return 0
            fi
        done
    done
    return 1
}

if have_font "times*.ttf" "LiberationSerif-Regular.*"; then
    echo "Font: Times New Roman or Liberation Serif found."
else
    echo "ERROR: neither Times New Roman nor Liberation Serif was found." >&2
    echo "main.tex needs one of them; Tectonic does not supply system text faces." >&2
    exit 5
fi

# Tectonic: explicit override via $TECTONIC, then PATH, then a common local install.
TECTONIC_BIN=""
if [ -n "${TECTONIC:-}" ] && [ -x "${TECTONIC}" ]; then
    TECTONIC_BIN="${TECTONIC}"
elif command -v tectonic >/dev/null 2>&1; then
    TECTONIC_BIN="$(command -v tectonic)"
elif [ -x "/d/Tools/Tectonic/tectonic" ]; then
    TECTONIC_BIN="/d/Tools/Tectonic/tectonic"
elif [ -x "/d/Tools/Tectonic/tectonic.exe" ]; then
    TECTONIC_BIN="/d/Tools/Tectonic/tectonic.exe"
fi

if [ -n "$TECTONIC_BIN" ]; then
    echo "Engine: Tectonic ($TECTONIC_BIN)"
    "$TECTONIC_BIN" main.tex
    cp -f main.pdf paper.pdf
    echo "Compilation complete. Output: paper.pdf"
    exit 0
fi

if ! command -v xelatex >/dev/null 2>&1; then
    echo "ERROR: neither Tectonic nor XeLaTeX was found." >&2
    echo "Install Tectonic (https://tectonic-typesetting.github.io) or a TeX distribution" >&2
    echo "that provides xelatex and bibtex, then rerun this script." >&2
    exit 4
fi

echo "Engine: XeLaTeX"

echo "First xelatex run..."
xelatex -interaction=nonstopmode -halt-on-error main.tex

echo "Running bibtex..."
if [ -f bibliography.bib ]; then
    bibtex main
fi

echo "Second xelatex run..."
xelatex -interaction=nonstopmode -halt-on-error main.tex

echo "Final xelatex run..."
xelatex -interaction=nonstopmode -halt-on-error main.tex

cp -f main.pdf paper.pdf
echo "Compilation complete. Output: paper.pdf"
