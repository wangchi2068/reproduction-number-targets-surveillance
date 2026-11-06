"""Static cross-reference, citation, figure and table checks for the manuscript.

Usage:  python code/check_manuscript.py <repo-root or paper-folder>
The paper lives in <root>/paper (main.tex, sections/*.tex, bibliography.bib, figures/, tables/).
Exits 1 if any reference, citation, figure, table or label problem is found (unused bib keys are reported only).
"""
import json
import re
import sys
from pathlib import Path

ARG = Path(sys.argv[1] if len(sys.argv) > 1 else '.')
FOLDER = ARG if (ARG / 'main.tex').exists() else ARG / 'paper'
CITE = re.compile(r'\\(?:cite|citep|citet|citealp|citealt|citeauthor|citeyear|citeyearpar|parencite|textcite)\*?'
                  r'(?:\[[^\]]*\])*\{([^}]+)\}')


def tex_sources(folder):
    files = [folder / 'main.tex'] + sorted((folder / 'sections').glob('*.tex')) + sorted((folder / 'tables').glob('*.tex'))
    return '\n'.join(p.read_text(encoding='utf-8') for p in files if p.exists())


def check(folder):
    text = tex_sources(folder)
    text_nc = re.sub(r'(?<!\\)%.*', '', text)  # drop comments
    labels = re.findall(r'\\label\{([^}]+)\}', text_nc)
    refs = set(re.findall(r'\\(?:ref|eqref|autoref|cref|Cref|pageref)\{([^}]+)\}', text_nc))
    cites = set()
    for group in CITE.findall(text_nc):
        cites.update(x.strip() for x in group.split(',') if x.strip())
    bib = folder / 'bibliography.bib'
    keys = set(re.findall(r'@\w+\{([^,\s]+),', bib.read_text(encoding='utf-8'))) if bib.exists() else set()
    figures = set(re.findall(r'includegraphics(?:\[[^\]]*\])?\{([^}]+)\}', text_nc))
    tables = set(re.findall(r'\\input\{(tables/[^}]+)\}', text_nc))

    def fig_exists(f):
        p = folder / f
        return p.exists() or any(p.with_suffix(s).exists() for s in ('.pdf', '.png', '.jpg'))

    return {
        'paper_folder': str(folder),
        'missing_refs': sorted(refs - set(labels)),
        'missing_citations': sorted(cites - keys),
        'unused_bib_keys': sorted(keys - cites),
        'missing_figures': sorted(f for f in figures if not fig_exists(f)),
        'missing_tables': sorted(t for t in tables if not (folder / (t if t.endswith('.tex') else t + '.tex')).exists()),
        'duplicate_labels': sorted({x for x in labels if labels.count(x) > 1}),
        'n_citations': len(cites), 'n_bib_keys': len(keys), 'n_labels': len(labels), 'n_refs': len(refs),
    }


report = check(FOLDER)
print(json.dumps(report, indent=2, ensure_ascii=False))
bad = any(report[k] for k in ('missing_refs', 'missing_citations', 'missing_figures', 'missing_tables', 'duplicate_labels'))
sys.exit(1 if bad else 0)
