"""Prespecified supplementary season and participation sensitivities."""
from __future__ import annotations
import argparse
from datetime import datetime
import json
from pathlib import Path
import numpy as np
import pandas as pd
from revision_data import EXP,RAW,load_sample,load_raw,build_panel,wl_R,sha256
from revision_stats import fit_loss
from revision_empirical import resample_counts,finite_summary


def run(output_dir,bootstrap=1000,seed=2026110306):
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    result={'meta':{'bootstrap':bootstrap,'seed':seed,'raw_sha256':sha256(RAW),'code_sha256':sha256(__file__)},'seasons':{},'participation':{},'summer':{}}
    for stat in ('RA','G'):
        frame=load_sample(stat)
        season_results={}
        for s,e in frame.groupby('season',sort=True):
            e=e.reset_index(drop=True);f=fit_loss(e.x,e.D);draws=[]
            rng=np.random.default_rng(np.random.SeedSequence([seed,0 if stat=='RA' else 1,int(s[:4])]))
            for b in range(bootstrap):
                w,_=resample_counts({stat:e},rng)[stat];g=fit_loss(e.x,e.D,w)
                draws.append(g['a'] if g['converged'] else np.nan)
            season_results[s]={'n':len(e),'a':finite_summary(f['a'],draws),'point':f}
        result['seasons'][stat]=season_results
        e=frame[(frame.coverage>0)&(frame.coverage_next>0)].reset_index(drop=True)
        adjusted_ct=e.C_t/(e.coverage/100);adjusted_cn=e.C_next/(e.coverage_next/100)
        if stat=='RA':
            xa=wl_R((adjusted_cn+.5)/(adjusted_ct+.5))
        else:
            raw=load_raw();lookup=raw.set_index(['state','week_end']).coverage
            c1=np.array([lookup.get((s,t-pd.Timedelta(days=7)),np.nan) for s,t in zip(e.state,e.week_end)])/100
            c2=np.array([lookup.get((s,t-pd.Timedelta(days=14)),np.nan) for s,t in zip(e.state,e.week_end)])/100
            good=(c1>0)&(c2>0)&np.isfinite(c1)&np.isfinite(c2)
            e=e[good].reset_index(drop=True);adjusted_ct=adjusted_ct.to_numpy()[good];adjusted_cn=adjusted_cn.to_numpy()[good]
            xa=adjusted_cn/(.65*adjusted_ct+.25*e.C_prev1.to_numpy()/c1[good]+.10*e.C_prev2.to_numpy()/c2[good])
        f=fit_loss(e.x,e.D);g=fit_loss(xa,e.D)
        draws=[];rng=np.random.default_rng(np.random.SeedSequence([seed,2,0 if stat=='RA' else 1]))
        for b in range(bootstrap):
            w,_=resample_counts({stat:e},rng)[stat];fo=fit_loss(e.x,e.D,w);fa=fit_loss(xa,e.D,w)
            draws.append(fo['a']-fa['a'] if fo['converged'] and fa['converged'] else np.nan)
        result['participation'][stat]={'n':len(e),'coverage_median':float(e.coverage.median()),'original_a':f['a'],'adjusted_a':g['a'],
                                     'delta_original_minus_adjusted':finite_summary(f['a']-g['a'],draws),
                                     'scale':'original scale held fixed; adjusts statistic counts using facility reporting participation',
                                     'scope':'descriptive paired sensitivity; does not measure all ascertainment or identify mandate effects'}
    p=build_panel();p=p[(~p.any_missing)&(p.D_t>=5)&(p.C_next>0)].copy();p['logG']=np.log(p.C_next/p.D_t)
    for year in (2023,2024,2025):
        e=p[p.week_end.between(f'{year}-05-01',f'{year}-10-31')]
        result['summer'][str(year)]={'n':len(e),'median_count_scale':float(e.D_t.median()),'coverage_median':float(e.coverage.median()),
                                    'iqr_logG':float(np.percentile(e.logG,75)-np.percentile(e.logG,25))}
    result['summer']['note']='One voluntary summer, 2024; unadjusted descriptive summaries only, not causal or independent replications.'
    result['meta']['completed_at']=datetime.now().astimezone().isoformat()
    (out/'secondary.json').write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output-dir',type=Path,default=EXP/'analysis/revised');p.add_argument('--bootstrap',type=int,default=1000)
    a=p.parse_args();run(a.output_dir,a.bootstrap)
