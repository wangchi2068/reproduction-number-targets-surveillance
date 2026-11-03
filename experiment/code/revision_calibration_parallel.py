"""Run deterministic disjoint calibration batches and validate their merge."""
from __future__ import annotations
import argparse,json,os
from datetime import datetime
from pathlib import Path
import subprocess,sys,time


def main():
    p=argparse.ArgumentParser();p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--replicates',type=int,default=200)
    p.add_argument('--bootstrap',type=int,default=2000);p.add_argument('--empirical-draws',type=Path,required=True);p.add_argument('--empirical-points',type=Path,required=True)
    p.add_argument('--workers',type=int,default=4)
    a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
    script=Path(__file__).with_name('revision_calibration.py')
    workers=min(a.workers,a.replicates);base,remainder=divmod(a.replicates,workers)
    batches=[];running=[];handles=[];start=0;begin=time.perf_counter()
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',PYTHONIOENCODING='utf-8')
    for i in range(workers):
        n=base+(i<remainder);folder=a.output_dir/'calibration_batches'/f'b{i}';folder.mkdir(parents=True,exist_ok=True)
        log=(folder/'execution.log').open('w',encoding='utf-8');handles.append(log)
        cmd=[sys.executable,'-B',str(script),'--output-dir',str(folder),'--replicates',str(n),'--replicate-start',str(start),'--simulate-only']
        running.append(subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,env=env));batches.append({'path':str(folder),'command':cmd,'replicate_start':start,'replicates':n});start+=n
    for i,process in enumerate(running):
        code=process.wait();handles[i].close();batches[i]['returncode']=code
    manifest={'started_profile':'independent deterministic replicate batches','batches':batches,'seconds':time.perf_counter()-begin,'completed_at':datetime.now().astimezone().isoformat()}
    (a.output_dir/'calibration_execution.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    if any(r['returncode'] for r in batches):raise SystemExit('Calibration batch failed; see per-batch execution.log')
    cmd=[sys.executable,'-B',str(script),'--output-dir',str(a.output_dir),'--bootstrap',str(a.bootstrap),'--empirical-draws',str(a.empirical_draws),'--empirical-points',str(a.empirical_points),'--merge-batches',*[r['path'] for r in batches]]
    subprocess.run(cmd,env=env,check=True)

if __name__=='__main__':main()
