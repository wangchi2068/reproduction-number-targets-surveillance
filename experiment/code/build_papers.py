"""Compile the two manuscripts and retain final-pass diagnostics outside source folders."""
from __future__ import annotations
import argparse
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

REPO=Path(__file__).resolve().parents[2]


def main():
    p=argparse.ArgumentParser();p.add_argument('--lang',choices=('en','zh','both'),default='both');p.add_argument('--output-dir',type=Path,default=REPO/'validation/build')
    args=p.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    engine=os.environ.get('TECTONIC') or shutil.which('tectonic')
    if not engine:raise SystemExit('Set TECTONIC to a Tectonic executable or install it on PATH. System Times New Roman/Liberation Serif and Chinese fonts are also required.')
    previous = args.output_dir / 'build_report.json'
    if previous.is_file():
        report = json.loads(previous.read_text(encoding='utf-8'))
        report['engine'] = str(engine)
        report['resumed_at'] = datetime.now().astimezone().isoformat()
    else:
        report = {'started_at': datetime.now().astimezone().isoformat(), 'engine': str(engine), 'languages': {}}
    failed = False
    for lang in (('en','zh') if args.lang=='both' else (args.lang,)):
        source=REPO/f'paper_{lang}'
        out=args.output_dir/lang;out.mkdir(exist_ok=True)
        command=[engine,'--keep-logs','--keep-intermediates','--print','--outdir',str(out),'main.tex']
        proc=subprocess.run(command,cwd=source,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        text=proc.stdout.decode('utf-8',errors='replace')
        (out/'console.log').write_text(text,encoding='utf-8')
        log=(out/'main.log').read_text(encoding='utf-8',errors='replace') if (out/'main.log').exists() else ''
        blg=(out/'main.blg').read_text(encoding='utf-8',errors='replace') if (out/'main.blg').exists() else ''
        patterns={'overfull_hbox':r'Overfull \\hbox','undefined_references':r'(?:Reference|Citation).*undefined|There were undefined references',
                  'missing_characters':r'Missing character:','latex_errors':r'^!|^.*?Emergency stop',
                  'duplicate_destinations':r'Object @.*already defined'}
        findings={k:re.findall(v,log,re.M|re.I) for k,v in patterns.items()}
        biberrors=bool(re.search(r'There (?:was|were) \d+ error|I was expecting|Illegal, another|Repeated entry',blg))
        success=proc.returncode==0 and (out/'main.pdf').exists() and not any(findings.values()) and not biberrors
        record={'command':command,'returncode':proc.returncode,'passed':success,'final_pass_findings':findings,'bibtex_errors':biberrors,
                'system_font_warnings':[line for line in text.splitlines() if 'absolute path' in line or 'Fontconfig error' in line],
                'console_log':str((out/'console.log').relative_to(REPO)) if out.is_relative_to(REPO) else str(out/'console.log')}
        if (out/'main.pdf').exists():
            from pypdf import PdfReader
            pdf=PdfReader(out/'main.pdf');record['pages']=len(pdf.pages)
            record['unresolved_question_marks']=sum((page.extract_text() or '').count('??') for page in pdf.pages)
            if record['unresolved_question_marks']:success=False;record['passed']=False
            if success:
                shutil.copy2(out/'main.pdf',source/'paper.pdf')
                record['pdf_sha256']=hashlib.sha256((source/'paper.pdf').read_bytes()).hexdigest()
        report['languages'][lang]=record;failed|=not success
        print(lang,'PASS' if success else 'FAIL',record.get('pages'),{k:len(v) for k,v in findings.items()},flush=True)
    report['completed_at']=datetime.now().astimezone().isoformat();report['passed']=not failed
    (args.output_dir/'build_report.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
    if failed:raise SystemExit(1)

if __name__=='__main__':main()
