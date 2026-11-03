"""Run every retained analysis with content-addressed stage validation.

Fresh runs use a new output directory and never consult legacy outputs. A resumed
stage must match its complete source/input/configuration fingerprint and all
recorded output hashes. This runner does not claim successful stages it did not run.
"""
from __future__ import annotations
import argparse
from datetime import datetime
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

CODE=Path(__file__).resolve().parent
REPO=CODE.parents[1]
EXP=CODE.parent


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()


def fingerprint(config):
    sources={p.relative_to(REPO).as_posix():digest(p) for p in sorted(CODE.glob('*.py')) if p.name.startswith('revision_') or p.name in ('s03b_thinning.py','s03c_state_fe.py')}
    inputs={p.relative_to(REPO).as_posix():digest(p) for p in sorted((EXP/'data/external').glob('*')) if p.is_file()}
    environment={name:importlib.metadata.version(name) for name in ('numpy','scipy','pandas','matplotlib','threadpoolctl')}
    payload={'sources':sources,'inputs':inputs,'configuration':config,'environment':environment,'python':sys.version}
    payload['fingerprint']=hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest()
    return payload


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',type=Path,default=EXP/'analysis/revised')
    p.add_argument('--fresh',action='store_true');p.add_argument('--resume',action='store_true')
    p.add_argument('--quick',action='store_true',help='Development run, never a publication result')
    p.add_argument('--render',action='store_true');p.add_argument('--compile',action='store_true')
    args=p.parse_args()
    if args.fresh and args.resume:raise SystemExit('Choose fresh or resume')
    out=args.output_dir.resolve()
    if args.fresh and out.exists():
        archive=out.parent/(out.name+'_archived_'+datetime.now().strftime('%Y%m%dT%H%M%S'))
        if archive.exists():raise SystemExit('Archive path already exists')
        shutil.move(str(out),str(archive))
    out.mkdir(parents=True,exist_ok=True);logs=out/'logs';logs.mkdir(exist_ok=True)
    q=args.quick
    cfg={'pooled_bootstrap':20 if q else 2000,'secondary_bootstrap':10 if q else 1000,
         'fe_bootstrap':10 if q else 1000,'thinning_bootstrap':8 if q else 1000,
         'point_reps':20 if q else 400,'inner_reps':8 if q else 50,'max_inner_reps':8 if q else 400,
         'calibration_panels':2 if q else 200,'calibration_bootstrap':20 if q else 2000,
         'design_replicates':5 if q else 500,'theory_replicates':2000 if q else 200000,'profile':'quick' if q else 'publication'}
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',PYTHONIOENCODING='utf-8',MPLBACKEND='Agg')
    fp=fingerprint(cfg);statefile=out/'run_manifest.json'
    state=json.loads(statefile.read_text(encoding='utf-8')) if statefile.exists() else {}
    if args.resume and state and state.get('fingerprint')!=fp['fingerprint']:
        raise SystemExit('Resume fingerprint differs: use --fresh or a new output directory')
    if state and not args.resume:raise SystemExit('Existing manifest: use --resume, --fresh or a new output directory')
    state={**fp,'started_at':datetime.now().astimezone().isoformat(),'stages':state.get('stages',{}),'status':'running'}
    def save():statefile.write_text(json.dumps(state,indent=2,ensure_ascii=False),encoding='utf-8')
    save()
    from revision_data import audit_inputs
    audit_inputs(out)
    common=['--output-dir',str(out)]
    stages=[
      ('empirical','revision_empirical.py',common+['--bootstrap',str(cfg['pooled_bootstrap']),'--secondary-bootstrap',str(cfg['secondary_bootstrap'])],['empirical.json','pooled_draws.npz','centering_draws_RA.npz','centering_draws_G.npz','pooled_block2_draws.npz','pooled_block3_draws.npz']),
      ('fixed_effects','s03c_state_fe.py',common+['--bootstrap',str(cfg['fe_bootstrap'])],['fixed_effects.json','fixed_effects_draws.npz']),
      ('fixed_effects_block2','s03c_state_fe.py',['--output-dir',str(out/'block2'),'--bootstrap',str(cfg['fe_bootstrap']),'--block-length','2'],['block2/fixed_effects.json','block2/fixed_effects_draws.npz']),
      ('fixed_effects_block3','s03c_state_fe.py',['--output-dir',str(out/'block3'),'--bootstrap',str(cfg['fe_bootstrap']),'--block-length','3'],['block3/fixed_effects.json','block3/fixed_effects_draws.npz']),
      ('thinning','s03b_thinning.py',common+['--bootstrap',str(cfg['thinning_bootstrap']),'--point-reps',str(cfg['point_reps']),'--inner-reps',str(cfg['inner_reps']),'--max-inner-reps',str(cfg['max_inner_reps'])],['thinning.json','thinning_draws.npz']),
      ('theory','revision_theory.py',common+['--replicates',str(cfg['theory_replicates'])],['theory.json']),
      ('secondary','revision_secondary.py',common+['--bootstrap',str(cfg['secondary_bootstrap'])],['secondary.json']),
      ('design','revision_design.py',common+['--replicates',str(cfg['design_replicates'])],['design.json','design_draws.npz']),
      ('calibration','revision_calibration_parallel.py',common+['--replicates',str(cfg['calibration_panels']),'--bootstrap',str(cfg['calibration_bootstrap']),'--empirical-draws',str(out/'pooled_draws.npz'),'--empirical-points',str(out/'empirical.json')],['calibration.json','calibration_panels.jsonl','calibration_inference.json','calibration_draws.csv'])]
    for name,script,cmdargs,outputs in stages:
        old=state['stages'].get(name)
        if args.resume and old and old.get('status')=='passed' and old.get('outputs') and all((out/f).is_file() and digest(out/f)==h for f,h in old['outputs'].items()):
            print('verified resume',name,flush=True);continue
        cmd=[sys.executable,'-B',str(CODE/script),*cmdargs]
        started=datetime.now().astimezone().isoformat();t=time.perf_counter()
        print('running',name,flush=True)
        with (logs/f'{name}.log').open('w',encoding='utf-8') as log:
            proc=subprocess.run(cmd,cwd=CODE,env=env,stdout=log,stderr=subprocess.STDOUT)
        missing=[f for f in outputs if not (out/f).is_file()]
        record={'command':cmd,'started_at':started,'ended_at':datetime.now().astimezone().isoformat(),'seconds':time.perf_counter()-t,
                'returncode':proc.returncode,'status':'passed' if proc.returncode==0 and not missing else 'failed','missing_outputs':missing,
                'outputs':{f:digest(out/f) for f in outputs if (out/f).is_file()},'log':f'logs/{name}.log'}
        state['stages'][name]=record;save()
        if record['status']!='passed':state['status']='failed';save();raise SystemExit(f'Stage {name} failed; see {logs/name}.log')
    state['status']='analyses_passed';state['finished_at']=datetime.now().astimezone().isoformat();save()
    if args.render:
        subprocess.run([sys.executable,'-B',str(CODE/'render_revision.py'),'--results-dir',str(out),'--repo-root',str(REPO)],env=env,check=True)
    if args.compile:
        for lang in ('en','zh'):
            subprocess.run(['bash',str(REPO/f'paper_{lang}/compile.sh')],cwd=REPO/f'paper_{lang}',env=env,check=True)
    print('completed',out,flush=True)


if __name__=='__main__':main()
