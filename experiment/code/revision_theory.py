"""Known-truth checks of the stated model, distinct from empirical validation."""
from __future__ import annotations
import argparse,json
from datetime import datetime
from pathlib import Path
import numpy as np
from revision_data import EXP,sha256


def run(output_dir,replicates=200000,seed=2026110307):
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    rb,vr,k,rho=1.15,.04,.35,.25
    nvalues=np.unique(np.rint(np.geomspace(8,20000,22)).astype(int));rows=[]
    for j,n in enumerate(nvalues):
        m=rho*n
        for tau in (0.,.15):
            rng=np.random.default_rng(np.random.SeedSequence([seed,j,int(tau*100)]))
            R=rng.gamma(rb*rb/vr,vr/rb,replicates)
            if tau:
                concentration=(1-rho)/(rho*tau*tau)-1
                p=rng.beta(rho*concentration,(1-rho)*concentration,replicates)
            else:p=np.full(replicates,rho)
            counts=rng.negative_binomial(n*k,k/(k+p*R))
            estimate=counts/m;vh=tau*tau*(rb*rb+vr);A=rb+rho*(1+tau*tau)*(rb*rb+vr)/k
            for target,errors,exact in [('current',(estimate-R)**2,vh+A/m),('process',(estimate-rb)**2,vr+vh+A/m)]:
                rows.append({'m':float(m),'n':int(n),'tau':tau,'target':target,'analytic':float(exact),
                             'monte_carlo':float(errors.mean()),'mc_se':float(errors.std(ddof=1)/np.sqrt(replicates))})
    finite=[]
    T=10;repT=replicates//10
    for j,m in enumerate((10.,50.,500.)):
        for tau in (0.,.15):
            rng=np.random.default_rng(np.random.SeedSequence([seed,10,j,int(tau*100)]))
            R=rng.gamma(rb*rb/vr,vr/rb,(repT,T));n=int(m/rho)
            if tau:
                con=(1-rho)/(rho*tau*tau)-1;p=rng.beta(rho*con,(1-rho)*con,(repT,T))
            else:p=np.full_like(R,rho)
            obs=rng.negative_binomial(n*k,k/(k+p*R))/m
            err=(obs.mean(axis=1)-R.mean(axis=1))**2
            vh=tau*tau*(rb*rb+vr);A=rb+rho*(1+tau*tau)*(rb*rb+vr)/k
            finite.append({'m':m,'T':T,'tau':tau,'analytic':(vh+A/m)/T,'monte_carlo':float(err.mean()),'mc_se':float(err.std(ddof=1)/np.sqrt(repT))})
    result={'meta':{'seed':seed,'replicates':replicates,'finite_replicates':repT,'code_sha256':sha256(__file__),'completed_at':datetime.now().astimezone().isoformat()},
            'parameters':{'Rbar':rb,'v_R':vr,'k':k,'rho':rho},'checks':rows,'finite_mean_checks':finite,
            'scope':'Checks formulas under their own fixed-exposure, independent reporting and conditional-unbiasedness assumptions; not empirical mechanism validation.'}
    (out/'theory.json').write_text(json.dumps(result,indent=2),encoding='utf-8');return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output-dir',type=Path,default=EXP/'analysis/revised');p.add_argument('--replicates',type=int,default=200000)
    a=p.parse_args();run(a.output_dir,a.replicates)
