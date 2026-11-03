"""Prune bibliography to keys actually cited in the English manuscript."""
import re
import sys
from pathlib import Path

root = Path(sys.argv[1])
text = '\n'.join(p.read_text(encoding='utf-8') for p in sorted((root / 'paper_en/sections').glob('*.tex')))
cites = set()
for group in re.findall(r'\\citep?\{([^}]+)\}', text):
    cites.update(x.strip() for x in group.split(','))
print('cited keys', len(cites))
for lang in ('en', 'zh'):
    path = root / f'paper_{lang}/bibliography.bib'
    raw = path.read_text(encoding='utf-8')
    blocks = re.split(r'(?m)(?=^@)', raw)
    head = ''.join(b for b in blocks if not b.strip().startswith('@'))
    entries = [b for b in blocks if b.strip().startswith('@')]
    kept, removed = [], []
    for block in entries:
        match = re.match(r'@\w+\{([^,]+),', block)
        key = match.group(1).strip() if match else ''
        if key in cites:
            kept.append(block if block.endswith('\n') else block + '\n')
        else:
            removed.append(key)
    path.write_text(head + ''.join(kept), encoding='utf-8')
    print(lang, 'kept', len(kept), 'removed', len(removed), removed)
print('identical', (root / 'paper_en/bibliography.bib').read_bytes() == (root / 'paper_zh/bibliography.bib').read_bytes())
