"""Recomputed pooled inference, honest validation and loss-definition diagnostics."""
from __future__ import annotations

import argparse
from datetime import datetime
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

from revision_data import EXP, RAW, load_sample, sha256
from revision_stats import fit_affine, fit_loss, relative_loss


def finite_summary(point, values, *, ratio=False):
    a = np.asarray(values, float)
    good = np.isfinite(a)
    out = {"value": float(point) if point is not None and np.isfinite(point) else None,
           "ci95": np.percentile(a[good], [2.5,97.5]).tolist() if good.any() else None,
           "valid_draws": int(good.sum()), "total_draws": len(a),
           "nonpositive_fraction": float(np.mean(a[good] <= 0)) if good.any() else None}
    if ratio and (not good.all() or out["nonpositive_fraction"] > .025):
        out["ci95"] = None
        out["interval_status"] = "unstable denominator; no finite ratio interval reported"
    return out


def fit_power(loss, d, weights=None):
    loss, d = np.asarray(loss,float),np.asarray(d,float)
    w=np.ones(len(loss)) if weights is None else np.asarray(weights,float)
    keep=(w>0)&np.isfinite(loss)&np.isfinite(d)&(d>0)
    loss,d,w=loss[keep],d[keep],w[keep]; w=w/w.sum(); ld=np.log(d)
    def objective(gamma):
        scaled=loss*np.exp(gamma*ld)
        c=max(float(w@scaled),1e-12)
        return np.log(c)-gamma*(w@ld)+1
    opt=minimize_scalar(objective,bounds=(0.,3.),method='bounded',options={'xatol':1e-10})
    gamma=float(opt.x); c=float(w@(loss*np.exp(gamma*ld)))
    return {'coefficient':c,'exponent':gamma,'converged':bool(opt.success),'boundary':bool(gamma<1e-5 or gamma>2.99999)}


def curve_validation(frame, group):
    x,d=frame.x.to_numpy(),frame.D.to_numpy(); rows=[]
    for label in sorted(frame[group].unique()):
        test=(frame[group].to_numpy()==label); train=~test
        loss,center=relative_loss(x[train]); af=fit_affine(loss,d[train]); pw=fit_power(loss,d[train])
        held=((x[test]-center)/center)**2
        ma=af['a']+af['b']/d[test]; mp=pw['coefficient']*d[test]**(-pw['exponent'])
        valid=bool(af['converged'] and pw['converged'] and np.all(ma>0) and np.all(mp>0))
        if valid:
            sa=float(np.mean(held/ma+np.log(ma)));sp=float(np.mean(held/mp+np.log(mp)))
            mse_a=float(np.mean((held-ma)**2));mse_p=float(np.mean((held-mp)**2))
        else:sa=sp=mse_a=mse_p=None
        rows.append({'fold':str(label),'test_n':int(test.sum()),'train_n':int(train.sum()),'training_center':float(center),
                     'converged_and_positive':valid,'score_affine':sa,'score_power':sp,
                     'power_minus_affine':sp-sa if valid else None,'mse_affine':mse_a,'mse_power':mse_p,
                     'affine_fit':af,'power_fit':pw})
    valid=[r for r in rows if r['converged_and_positive']]
    total=sum(r['test_n'] for r in valid)
    return {'group':group,'folds':rows,'valid_folds':len(valid),'total_folds':len(rows),
            'mean_score_difference':sum(r['test_n']*r['power_minus_affine'] for r in valid)/total if total else None,
            'affine_better_folds':sum(r['power_minus_affine']>0 for r in valid),
            'score':'mean L/mu+log(mu); lower is better; positive difference favors affine',
            'uncertainty_note':'Fold values share training data; no independent-fold standard error or four-season significance test.'}


def group_center(values, groups, weights):
    codes,_=pd.factorize(groups,sort=True)
    den=np.bincount(codes,weights=weights)
    num=np.bincount(codes,weights=weights*values)
    means=np.divide(num,den,out=np.zeros_like(num),where=den>0)
    return means[codes]


def centered_loss(frame, mode, weights=None):
    x=frame.x.to_numpy();w=np.ones(len(x)) if weights is None else np.asarray(weights,float)
    xb=float(np.average(x,weights=w))
    if mode=='pooled':center=np.full(len(x),xb)
    elif mode=='season':center=group_center(x,frame.season,w)
    elif mode=='season_linear':
        center=np.zeros(len(x))
        for season in sorted(frame.season.unique()):
            mask=(frame.season.to_numpy()==season); positive=mask&(w>0)
            if not positive.any():continue
            t=(frame.week_end-frame.loc[mask,'week_end'].min()).dt.days.to_numpy()/7
            X=np.column_stack([np.ones(len(x)),t])
            beta=np.linalg.lstsq(X[positive]*np.sqrt(w[positive,None]),x[positive]*np.sqrt(w[positive]),rcond=None)[0]
            center[mask]=X[mask]@beta
    elif mode in ('week','week_leave_state_out'):
        center=group_center(x,frame.week_end,w)
        if mode=='week_leave_state_out':
            codes,_=pd.factorize(frame.week_end,sort=True)
            den=np.bincount(codes,weights=w)[codes]-w
            num=np.bincount(codes,weights=w*x)[codes]-w*x
            center=np.divide(num,den,out=np.full(len(x),np.nan),where=den>0)
    else:raise ValueError(mode)
    loss=((x-center)/xb)**2
    valid=np.isfinite(loss)&(w>0)
    if not valid.any():return None
    return loss,valid,xb


def vr_week(frame, weights=None, week_multiplicity=None):
    x=frame.x.to_numpy();w=np.ones(len(x)) if weights is None else weights
    xb=float(np.average(x,weights=w));d=x/xb-1
    codes,unique=pd.factorize(frame.week_end,sort=True)
    den=np.bincount(codes,weights=w)
    means=np.divide(np.bincount(codes,weights=w*d),den,out=np.zeros(len(unique)),where=den>0)
    valid=den>0
    ww=np.ones(len(unique)) if week_multiplicity is None else week_multiplicity
    ww=ww[valid];m=means[valid]
    between=float(np.average((m-np.average(m,weights=ww))**2,weights=ww))
    total=float(np.average((d-np.average(d,weights=w))**2,weights=w))
    return between/total if total>0 else np.nan


def resample_counts(frames, rng, block_length=1):
    union=pd.concat([f[['state','week_end','season']] for f in frames.values()]).drop_duplicates()
    states=sorted(union.state.unique());weeks=sorted(union.week_end.unique())
    sm=rng.multinomial(len(states),np.ones(len(states))/len(states))
    if block_length==1:
        wm=rng.multinomial(len(weeks),np.ones(len(weeks))/len(weeks))
    else:
        wm=np.zeros(len(weeks),int);wi={t:i for i,t in enumerate(weeks)}
        for season,g in union[['week_end','season']].drop_duplicates().groupby('season',sort=True):
            indices=np.array([wi[t] for t in sorted(g.week_end)]);n=len(indices)
            selected=[]
            while len(selected)<n:
                start=int(rng.integers(n))
                selected.extend(indices[(start+np.arange(block_length))%n].tolist())
            wm+=np.bincount(selected[:n],minlength=len(weeks))
    sdict=dict(zip(states,sm));wdict=dict(zip(weeks,wm))
    result={}
    for key,f in frames.items():
        w=f.state.map(sdict).to_numpy()*f.week_end.map(wdict).to_numpy()
        result[key]=(w,np.array([wdict[t] for t in sorted(f.week_end.unique())]))
    return result


def run(output_dir, bootstrap=2000, secondary_bootstrap=1000, seed=2026110301):
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    frames={'RA':load_sample('RA'),'G':load_sample('G'),'RA32':load_sample('RA',3.2),'RA32source':load_sample('RA',3.2,gi_sd=2.1)}
    result={'meta':{'analysis':'pooled, foldwise validation and centering','seed':seed,'bootstrap':bootstrap,
                    'secondary_bootstrap':secondary_bootstrap,'raw_sha256':sha256(RAW),'code_sha256':sha256(__file__)},
            'pooled':{},'validation':{},'centering':{},'block_sensitivity':{}}
    names=('a','b','ratio','xbar')
    draws={f'{key}_{v}':np.full(bootstrap,np.nan) for key in frames for v in names}
    for key,e in frames.items():
        f=fit_loss(e.x.to_numpy(),e.D.to_numpy())
        result['pooled'][key]={'point':f,'n':len(e),'states':int(e.state.nunique()),'weeks':int(e.week_end.nunique())}
        result['validation'][key]={g:curve_validation(e,g) for g in ('state','season')}
        loss,_=relative_loss(e.x.to_numpy())
        pw=fit_power(loss,e.D.to_numpy())
        bins=pd.DataFrame({'D':e.D,'L':loss,'bin':pd.qcut(e.D,12,labels=False,duplicates='drop')}).groupby('bin').agg(D=('D','mean'),L=('L','mean'),n=('L','size')).reset_index()
        result['pooled'][key]['power_fit']=pw
        result['pooled'][key]['bins']=bins.to_dict('records')
    rng=np.random.default_rng(seed)
    for b in range(bootstrap):
        weights=resample_counts(frames,rng)
        for key,e in frames.items():
            w,_=weights[key]
            f=fit_loss(e.x.to_numpy(),e.D.to_numpy(),w)
            if f['converged']:
                for v in names:draws[f'{key}_{v}'][b]=f[v] if f[v] is not None else np.nan
        if (b+1)%250==0:print('pooled bootstrap',b+1,flush=True)
    for key in frames:
        for name in names:
            result['pooled'][key][name]=finite_summary(result['pooled'][key]['point'][name],draws[f'{key}_{name}'],ratio=name=='ratio')
    np.savez(out/'pooled_draws.npz',**draws)
    modes=('pooled','season','season_linear','week','week_leave_state_out')
    for key in ('RA','G'):
        e=frames[key];point={};dd={mode:np.full((secondary_bootstrap,2),np.nan) for mode in modes}
        for mode in modes:
            L,valid,xb=centered_loss(e,mode)
            point[mode]=fit_affine(L[valid],e.D.to_numpy()[valid])
        for b in range(secondary_bootstrap):
            w,wm=resample_counts({key:e},np.random.default_rng(np.random.SeedSequence([seed,2,0 if key=='RA' else 1,b])))[key]
            for mode in modes:
                z=centered_loss(e,mode,w)
                if z is None:continue
                L,valid,_=z;f=fit_affine(L[valid],e.D.to_numpy()[valid],w[valid])
                if f['converged']:dd[mode][b]=[f['a'],f['b']]
        result['centering'][key]={mode:{'a':finite_summary(point[mode]['a'],dd[mode][:,0]),
                                               'b':finite_summary(point[mode]['b'],dd[mode][:,1]),
                                               'point_fit':point[mode]} for mode in modes}
        result['centering'][key]['week_change']=finite_summary(point['week']['a']-point['pooled']['a'],dd['week'][:,0]-dd['pooled'][:,0])
        result['centering'][key]['week_reduction_fraction']=finite_summary(1-point['week']['a']/point['pooled']['a'],1-dd['week'][:,0]/dd['pooled'][:,0])
        result['centering'][key]['VR_week']=float(vr_week(e))
        np.savez(out/f'centering_draws_{key}.npz',**dd)
    for length in (2,3):
        bd={k:np.full((secondary_bootstrap,2),np.nan) for k in ('RA','G')}
        rng=np.random.default_rng(np.random.SeedSequence([seed,3,length]))
        for b in range(secondary_bootstrap):
            weights=resample_counts(frames,rng,length)
            for key in bd:
                e=frames[key];f=fit_loss(e.x.to_numpy(),e.D.to_numpy(),weights[key][0])
                if f['converged']:bd[key][b]=[f['a'],f['b']]
        result['block_sensitivity'][str(length)]={key:{name:finite_summary(result['pooled'][key]['point'][name],bd[key][:,i]) for i,name in enumerate(('a','b'))} for key in bd}
        np.savez(out/f'pooled_block{length}_draws.npz',**bd)
    result['common_sample_n']=len(pd.merge(frames['RA'][['state','week_end']],frames['G'][['state','week_end']],on=['state','week_end']))
    result['meta']['completed_at']=datetime.now().astimezone().isoformat()
    (out/'empirical.json').write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False),encoding='utf-8')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output-dir',type=Path,default=EXP/'analysis/revised')
    parser.add_argument('--bootstrap',type=int,default=2000)
    parser.add_argument('--secondary-bootstrap',type=int,default=1000)
    args=parser.parse_args()
    run(args.output_dir,args.bootstrap,args.secondary_bootstrap)
