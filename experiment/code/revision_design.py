"""Illustrative curve-discrimination design simulation, not a biological floor test."""
from __future__ import annotations
import argparse
from datetime import datetime
import json
from pathlib import Path
import numpy as np
from scipy.stats import beta
from revision_data import EXP, load_sample, RAW, sha256
from revision_empirical import fit_power
from revision_stats import fit_affine, relative_loss


def binomial_interval(k,n):
    return [float(beta.ppf(.025,k,n-k+1)) if k else 0.,float(beta.ppf(.975,k+1,n-k)) if k<n else 1.]


def make_design(frame,seasons,multiplier,rng):
    labels=sorted(frame.season.unique()); states={s:i for i,s in enumerate(sorted(frame.state.unique()))}
    d=[];si=[];wi=[];lengths=[];offset=0
    for j in range(seasons):
        source=labels[j%len(labels)]
        f=frame[frame.season==source];weeks={w:i for i,w in enumerate(sorted(f.week_end.unique()))}
        d.extend((f.D.to_numpy()*multiplier).tolist());si.extend(f.state.map(states));wi.extend(f.week_end.map(weeks)+offset)
        lengths.append(len(weeks));offset+=len(weeks)
    return np.asarray(d),np.asarray(si,int),np.asarray(wi,int),lengths


def multiplier_noise(sc,wc,lengths,rng,state_sd=.35,week_sd=.35,individual_sd=.9,phi=.25):
    state=rng.normal(0,state_sd,int(sc.max()+1));week=np.zeros(sum(lengths));offset=0
    for size in lengths:
        z=rng.standard_normal(size);z[0]*=week_sd
        for j in range(1,size):z[j]=phi*z[j-1]+np.sqrt(1-phi*phi)*week_sd*z[j]
        week[offset:offset+size]=z;offset+=size
    eps=rng.normal(0,individual_sd,len(sc))
    variance=state_sd**2+week_sd**2+individual_sd**2
    return np.exp(state[sc]+week[wc]+eps-.5*variance)


def contrast(loss,d):
    a=fit_affine(loss,d);p=fit_power(loss,d)
    ma=a['a']+a['b']/d;mp=p['coefficient']*d**(-p['exponent'])
    if not a['converged'] or not p['converged'] or np.any(ma<=0):return np.nan
    return float(np.mean(loss/mp+np.log(mp))-np.mean(loss/ma+np.log(ma)))


def run(output_dir,replicates=500,seed=2026110304):
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    frame=load_sample();loss,_=relative_loss(frame.x.to_numpy())
    af=fit_affine(loss,frame.D.to_numpy());pw=fit_power(loss,frame.D.to_numpy())
    records=[];all_draws={}
    for ns in (4,8,12,20):
        for mult in (1,2,4):
            d,sc,wc,lengths=make_design(frame,ns,mult,np.random.default_rng(seed))
            mu0=pw['coefficient']*d**(-pw['exponent']);mu1=af['a']+af['b']/d
            draws=[]
            for stage in range(3):
                z=np.empty(replicates)
                for r in range(replicates):
                    rng=np.random.default_rng(np.random.SeedSequence([seed,ns,mult,stage,r]))
                    y=(mu1 if stage==2 else mu0)*multiplier_noise(sc,wc,lengths,rng)
                    z[r]=contrast(y,d)
                draws.append(z)
            cal,null,alternative=draws
            crit=float(np.nanquantile(cal,.95));good0=np.isfinite(null);good1=np.isfinite(alternative)
            k0=int(np.sum(null[good0]>crit));k1=int(np.sum(alternative[good1]>crit))
            rec={'seasons':ns,'count_multiplier':mult,'rows':len(d),'critical_value':crit,'calibration_valid':int(np.isfinite(cal).sum()),
                 'null_valid':int(good0.sum()),'alternative_valid':int(good1.sum()),
                 'size':k0/int(good0.sum()),'size_ci95':binomial_interval(k0,int(good0.sum())),
                 'power':k1/int(good1.sum()),'power_ci95':binomial_interval(k1,int(good1.sum()))}
            rec['size_compatible_with_nominal']=bool(rec['size_ci95'][0]<=.05<=rec['size_ci95'][1])
            records.append(rec)
            for name,z in zip(('calibration','null','alternative'),draws):all_draws[f'S{ns}_m{mult}_{name}']=z
            print('design',ns,mult,'size',rec['size'],'power',rec['power'],flush=True)
    result={'meta':{'replicates_per_stage':replicates,'seed':seed,'raw_sha256':sha256(RAW),'code_sha256':sha256(__file__),
                     'completed_at':datetime.now().astimezone().isoformat()},
            'generator':{'description':'Conditional relative-loss curves on recycled observed season/count designs, crossed mean-one lognormal multipliers. No inferred transmission truth.',
                         'state_sd':.35,'week_sd':.35,'individual_sd':.9,'weekly_ar':.25,
                         'affine':af,'power':pw,'nominal_alpha':.05,'calibration_and_evaluation_independent':True,
                         'scope':'Illustrative model-conditional discrimination; not a universal required season count, not retrospective observed power.'},'designs':records}
    np.savez(out/'design_draws.npz',**all_draws)
    (out/'design.json').write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output-dir',type=Path,default=EXP/'analysis/revised');p.add_argument('--replicates',type=int,default=500)
    a=p.parse_args();run(a.output_dir,a.replicates)
