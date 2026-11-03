"""Static cross-reference, citation, figure and table checks for both manuscripts."""
import json
import re
import sys
from pathlib import Path

ROOT = Path(sys.argv[1] if len(sys.argv) > 1 else '.')
PAPERS = {'en': ROOT / 'paper_en', 'zh': ROOT / 'paper_zh'}
report = {}
for lang, folder in PAPERS.items():
    text = '\n'.join(p.read_text(encoding='utf-8') for p in sorted((folder / 'sections').glob('*.tex')))
    main = (folder / 'main.tex').read_text(encoding='utf-8')
    labels = set(re.findall(r'\\label\{([^}]+)\}', text + main))
    refs = set(re.findall(r'\\(?:ref|eqref)\{([^}]+)\}', text + main))
    cites = set()
    for group in re.findall(r'\\citep?\{([^}]+)\}', text):
        cites.update(x.strip() for x in group.split(','))
    keys = set(re.findall(r'@\w+\{([^,]+),', (folder / 'bibliography.bib').read_text(encoding='utf-8')))
    figures = set(re.findall(r'includegraphics(?:\[[^\]]*\])?\{([^}]+)\}', text))
    tables = set(re.findall(r'\\input\{tables/([^}]+)\}', text))
    macros = set(re.findall(r'\\([A-Za-z]+Point|[A-Za-z]+Interval|Nra|Ngs|Nweeks|Nstates|FEvalid|FEBootstrap|ThinBootstrap|ThinPointReps|ThinInnerReps|ThinMCMax|ThinMCStatus|FEBlockTwoValid|FEBlockThreeValid)\b', text + main))
    defined = set(re.findall(r'\\newcommand\{\\([A-Za-z]+)\}', (folder / 'generated' / 'numbers.tex').read_text(encoding='utf-8'))) if (folder / 'generated' / 'numbers.tex').exists() else set()
    report[lang] = {
        'missing_refs': sorted(refs - labels),
        'missing_citations': sorted(cites - keys),
        'unused_bib_keys': sorted(keys - cites),
        'missing_figures': sorted(f for f in figures if not (folder / f).exists()),
        'missing_tables': sorted(t for t in tables if not (folder / 'tables' / (t + '.tex')).exists()),
        'undefined_number_macros': sorted(m for m in macros if m not in defined),
        'duplicate_labels': sorted({x for x in labels if (text + main).count('\\label{' + x + '}') > 1}),
    }
print(json.dumps(report, indent=2, ensure_ascii=False))
