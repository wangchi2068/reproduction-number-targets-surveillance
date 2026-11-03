"""Nested controlled thinning of the frozen RA sample; import performs no IO.

Each unique needed state-week count is drawn exactly once per realization and
reused wherever it occurs. Each realization is centred and FIT independently;
point estimates average 400 fitted coefficients and each bootstrap averages at
least 50 fits. Covariance, failures, denominator signs and inner Monte Carlo
error are retained. The Poisson reference is explicitly a delta approximation.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from revision_stats import (VERSION, bootstrap_weights, ci_summary, file_hash,
                            fit_loss, fit_loss_batch, sample_hash, save_draws, write_json)

QS = np.array([.8,.6,.5,.4,.3,.2,.1])
SUBSETS = ("full", "highCt50", "both50")
SEED = 20261103
METHOD_CODES = {"":0,"batch_damped_fgls":1,"damped_fgls":2,"gamma_objective_endpoint_fallback":3,"fgls":4,"failed":5}
FAILURE_CODES = {"":0,"no_positive_weight_rows":1,"invalid_relative_loss":2,"rank_deficient_slope":3,
                 "zero_loss_has_no_positive_gamma_mean_solution":4,"gamma_objective_fallback_failed":5}


def count_vector(frame):
    """Minimal unique count vector and index maps; reject conflicting overlaps."""
    records = {}
    for state, week, ct, cn in frame[["state","week_end","C_t","C_next"]].itertuples(index=False, name=None):
        for date, count in ((pd.Timestamp(week),ct),(pd.Timestamp(week)+pd.Timedelta(days=7),cn)):
            if not np.isfinite(count) or count < 0 or count != int(count):
                raise ValueError("Thinning requires finite nonnegative integer counts")
            key = (state,date)
            if key in records and records[key] != int(count):
                raise ValueError(f"Conflicting shared count at {key}")
            records[key] = int(count)
    keys = sorted(records)
    lookup = {key:i for i,key in enumerate(keys)}
    current = np.array([lookup[(s,pd.Timestamp(t))] for s,t in zip(frame.state,frame.week_end)])
    following = np.array([lookup[(s,pd.Timestamp(t)+pd.Timedelta(days=7))] for s,t in zip(frame.state,frame.week_end)])
    return keys, np.array([records[k] for k in keys],dtype=np.int64), current, following


def thinning_realizations(counts, current, following, q, reps, rng, active=None):
    """Draw only keys required by active rows; omitted rows receive NaN."""
    active = np.ones(len(current),bool) if active is None else np.asarray(active,bool)
    needed = np.unique(np.r_[current[active],following[active]])
    draws = rng.binomial(counts[needed],q,size=(reps,len(needed)))
    reverse = np.full(len(counts),-1,dtype=int)
    reverse[needed] = np.arange(len(needed))
    ct, cn = np.full((reps,len(current)),np.nan), np.full((reps,len(current)),np.nan)
    ct[:,active] = draws[:,reverse[current[active]]]
    cn[:,active] = draws[:,reverse[following[active]]]
    return ct,cn


def poisson_reference(frame, w=None, mean=3.6, sd=None):
    """Approximate bP = weighted mean[kappa(g)^2*(1+(Ct+.5)/(Cn+.5))].

    kappa=d log(RA)/d log(g)=(mean/7)/(1+(sd^2/mean)*log(g)/7).
    This plug-in elasticity approximation replaces row-specific RA/xbar by one
    and uses half-count stabilized inverse ratios. It is not the exact delta
    variance of the normalized nonlinear statistic, nor an exact Poisson
    invariant. It is recomputed with each sample's row multiplicities.
    """
    sd = mean*1.6/3.6 if sd is None else sd
    ct,cn = frame.C_t.to_numpy(float),frame.C_next.to_numpy(float)
    weights = np.ones(len(frame)) if w is None else np.asarray(w,float)
    use = weights > 0
    if not use.any():
        return np.nan
    g = (cn[use]+.5)/(ct[use]+.5)
    kappa = (mean/7)/(1+(sd**2/mean)*np.log(g)/7)
    return float(np.average(kappa**2*(1+1/g),weights=weights[use]))


def subset_masks(frame):
    return [np.ones(len(frame),bool), frame.C_t.to_numpy()>=50,
            (frame.C_t.to_numpy()>=50)&(frame.C_next.to_numpy()>=50)]


def _baseline(frame, w, masks):
    fits, references = [],[]
    for mask in masks:
        fits.append(fit_loss(frame.x.to_numpy()[mask],frame.C_t.to_numpy()[mask],w[mask]))
        references.append(poisson_reference(frame.loc[mask],w[mask]))
    return fits,np.array(references)


def _inner(frame, counts, current, following, masks, w, reps, seed, phase, draw_id, chunk):
    from revision_data import wl_R
    shape = (len(masks),len(QS),reps)
    coefficients = np.full(shape+(2,),np.nan)
    converged, iterations = np.zeros(shape,bool),np.zeros(shape,int)
    min_mu, methods, failures = np.full(shape,np.nan),np.zeros(shape,dtype=np.uint8),np.zeros(shape,dtype=np.uint8)
    for qi,q in enumerate(QS):
        rng = np.random.default_rng(np.random.SeedSequence([seed,phase,draw_id,qi,chunk]))
        ct,cn = thinning_realizations(counts,current,following,q,reps,rng,w>0)
        x = wl_R((cn+.5)/(ct+.5))
        for si,mask in enumerate(masks):
            if not np.any(w[mask]>0):
                failures[si,qi] = FAILURE_CODES["no_positive_weight_rows"]
                continue
            f = fit_loss_batch(x[:,mask],frame.C_t.to_numpy()[mask],w[mask])
            coefficients[si,qi,:,0],coefficients[si,qi,:,1] = f["a"],f["b"]
            converged[si,qi],iterations[si,qi],min_mu[si,qi] = f["converged"],f["iterations"],f["min_mu"]
            methods[si,qi] = [METHOD_CODES.get(m,255) for m in f["method"]]
            failures[si,qi] = [FAILURE_CODES.get(m,255) for m in f["failure"]]
            if np.any(methods[si,qi]==255) or np.any(failures[si,qi]==255):
                raise ValueError(f"Unmapped fit diagnostic: {set(f['method'])}, {set(f['failure'])}")
    return dict(coefficients=coefficients,converged=converged,iterations=iterations,min_mu=min_mu,methods=methods,failures=failures)


def _average(coefficients):
    # An incomplete realization set is not silently averaged over successes.
    return np.mean(coefficients,axis=-2)


def _rate_specs():
    specs = {"full":np.ones(len(QS),bool),"noq01":QS>=.2}
    specs.update({f"leave_q{q:g}":QS!=q for q in QS})
    return specs


def _metrics(coefficients, baseline, bp):
    """Map mean per-q coefficients to paired changes, rates and ratios."""
    result = {}
    inv = 1/QS-1
    for si,subset in enumerate(SUBSETS):
        a1,b1 = baseline[si]["a"],baseline[si]["b"]
        a1 = np.nan if a1 is None else a1
        b1 = np.nan if b1 is None else b1
        aa,bb = coefficients[si,:,0],coefficients[si,:,1]
        db,da = bb-b1,aa-a1
        result[f"{subset}/unthinned_a"] = a1
        result[f"{subset}/unthinned_b"] = b1
        result[f"{subset}/bP"] = bp[si]
        for qi,q in enumerate(QS):
            prefix = f"{subset}/q{q:g}"
            for name,value in (("a",aa[qi]),("b",bb[qi]),("c",db[qi]/inv[qi]),("delta_a",da[qi])):
                result[f"{prefix}/{name}"] = value
            result[f"{prefix}/c_over_bP"] = db[qi]/inv[qi]/bp[si] if bp[si]>0 else np.nan
        for label,mask in _rate_specs().items():
            c = float(inv[mask]@db[mask]/(inv[mask]@inv[mask]))
            da_rate = float(inv[mask]@da[mask]/(inv[mask]@inv[mask]))
            result[f"{subset}/{label}/c"] = c
            result[f"{subset}/{label}/delta_a_rate"] = da_rate
            result[f"{subset}/{label}/c_over_bP"] = c/bp[si] if bp[si]>0 else np.nan
        # Paired uncertainty for q exclusions themselves, not only separate CIs.
        for label in _rate_specs():
            if label != "full":
                result[f"{subset}/{label}/c_change_vs_full"] = result[f"{subset}/{label}/c"]-result[f"{subset}/full/c"]
    return result


def _metric_realizations(coefficients, baseline, bp, names):
    return np.array([[v.get(k,np.nan) for k in names]
                     for v in (_metrics(coefficients[:,:,r,:],baseline,bp) for r in range(coefficients.shape[2]))])


def _mc_variance(coefficients, baseline, bp, names):
    realizations = _metric_realizations(coefficients,baseline,bp,names)
    return np.var(realizations,axis=0,ddof=1)/len(realizations) if len(realizations)>1 else np.full(len(names),np.nan)


def _simultaneous(values, point):
    """Centred bootstrap max-|z| bands over all seven q values jointly."""
    values,point = np.asarray(values,float),np.asarray(point,float)
    use = np.isfinite(values).all(axis=1)
    if use.sum()<2 or not np.isfinite(point).all():
        return dict(ci95=None,complete_draws=int(use.sum()),finite_fraction=float(use.mean()),failure="insufficient_complete_draws")
    v = values[use]
    sd = v.std(axis=0,ddof=1)
    centred = v-v.mean(axis=0)
    z = np.divide(centred,sd,out=np.zeros_like(centred),where=sd>0)
    critical = float(np.quantile(np.max(np.abs(z),axis=1),.95))
    return dict(ci95=np.column_stack([point-critical*sd,point+critical*sd]),critical=critical,
                standard_errors=sd,complete_draws=int(use.sum()),finite_fraction=float(use.mean()),
                interval_scope="all_draws" if use.all() else "conditional_on_complete_draws",
                method="centred_bootstrap_max_absolute_standardized_deviation")


def known_count_check(seed=SEED,reps=100):
    """Known Poisson count generator, followed by the identical thinning fitter.

    This small conditional mechanism check is not a coverage or rejection-rate
    calibration. Latent intensities are deterministic state-specific exponential
    trajectories; observations are sampled once from Poisson(latent intensity).
    """
    from revision_data import wl_R
    rng = np.random.default_rng(np.random.SeedSequence([seed,3]))
    dates = pd.date_range("2020-10-03",periods=9,freq="7D")
    rows,latent,observed = [],[],[]
    for s in range(12):
        lam = (70+25*s)*1.12**np.arange(9)
        counts = rng.poisson(lam)
        latent.append(lam);observed.append(counts)
        for t in range(8):
            rows.append((f"K{s:02d}",dates[t],"known",counts[t],counts[t+1],counts[t],float(wl_R((counts[t+1]+.5)/(counts[t]+.5))),counts[t]))
    frame = pd.DataFrame(rows,columns=["state","week_end","season","C_t","C_next","D_t","x","D"])
    _,counts,current,following = count_vector(frame)
    masks = subset_masks(frame)
    base,bp = _baseline(frame,np.ones(len(frame)),masks)
    inner = _inner(frame,counts,current,following,masks,np.ones(len(frame)),reps,seed,4,0,0)
    metrics = _metrics(_average(inner["coefficients"]),base,bp)
    return dict(generator="independent_Poisson_known_latent_exponential_counts",
                seed_components=[seed,3],thinning_seed_phase=4,repetitions=reps,
                interpretation="conditional mechanism check only; not interval-coverage validation",
                full_sample={k:v for k,v in metrics.items() if k.startswith("full/")},
                fit_failures=int((~inner["converged"]).sum())),dict(known_latent=np.array(latent),known_observed=np.array(observed),
                    known_inner_coefficients=inner["coefficients"],known_inner_converged=inner["converged"])


def run_analysis(frame, bootstrap=1000, point_reps=400, inner_reps=50, seed=SEED,
                 block_length=1, mc_tolerance=.10, max_inner_reps=200, progress=False):
    if min(bootstrap,point_reps,inner_reps)<1 or mc_tolerance<=0:
        raise ValueError("Positive bootstrap/repetition counts and MC tolerance required")
    max_inner_reps = max(inner_reps,max_inner_reps)
    masks = subset_masks(frame)
    if any(not m.any() for m in masks):
        raise ValueError("A prespecified count subset is empty")
    keys,counts,current,following = count_vector(frame)
    base,bp = _baseline(frame,np.ones(len(frame)),masks)
    point_inner = _inner(frame,counts,current,following,masks,np.ones(len(frame)),point_reps,seed,0,0,0)
    point = _metrics(_average(point_inner["coefficients"]),base,bp)
    names = list(point)
    nstates,nweeks = frame.state.nunique(),frame.week_end.nunique()
    sm_all,wm_all = np.zeros((bootstrap,nstates),int),np.zeros((bootstrap,nweeks),int)
    row_weights,baseline_fits,references,all_inner = [],[],[],[]
    for b in range(bootstrap):
        rng = np.random.default_rng(np.random.SeedSequence([seed,2,b]))
        w,sm,wm = bootstrap_weights(frame,rng,block_length)
        sm_all[b],wm_all[b] = sm,wm
        row_weights.append(w)
        baseline,reference = _baseline(frame,w,masks)
        baseline_fits.append(baseline);references.append(reference)
        all_inner.append(_inner(frame,counts,current,following,masks,w,inner_reps,seed,1,b,0))
        if progress and (b+1)%max(1,bootstrap//10)==0:
            print(f"thinning bootstrap {b+1}/{bootstrap}, R={inner_reps}",flush=True)
    actual_reps,chunk_sizes,adaptation = inner_reps,[inner_reps],[]
    # Monitor the primary response summaries only; every subset-by-q rate is
    # reported, but demanding a fixed relative precision for all of them would
    # let one noisy secondary cell dictate the whole budget.
    monitor = [i for i,k in enumerate(names) if k in ("full/full/c","full/noq01/c","full/q0.1/c","highCt50/full/c","both50/full/c")]
    while True:
        metrics = [_metrics(_average(r["coefficients"]),f,p) for r,f,p in zip(all_inner,baseline_fits,references)]
        values = np.array([[v[k] for k in names] for v in metrics])
        mc_var = np.array([_mc_variance(r["coefficients"],f,p,names) for r,f,p in zip(all_inner,baseline_fits,references)])
        total_var = np.full(len(names),np.nan)
        inner_var = np.full(len(names),np.nan)
        for j in range(len(names)):
            use = np.isfinite(values[:,j]) & np.isfinite(mc_var[:,j])
            if use.sum()>1:
                total_var[j] = np.var(values[use,j],ddof=1)
                inner_var[j] = np.mean(mc_var[use,j])
        sampling_var = np.maximum(total_var-inner_var,0)
        # Precision criterion (protocol amendment 2026-11-03): the original
        # MC-SE/sampling-SE ratio grows with the outer budget and is unattainable
        # at any feasible inner budget. Use a scale-free relative precision:
        # MC standard error no more than 5% of the point estimate.
        scale = np.maximum(np.abs(np.array([point[k] for k in names])),1e-6)
        ratios = np.sqrt(np.divide(inner_var,scale**2,out=np.full(len(names),np.inf),where=scale>0))
        monitored = ratios[monitor]
        worst = float(np.max(monitored)) if len(monitored) else np.inf
        adaptation.append(dict(inner_reps=actual_reps,max_mc_se_over_sampling_se=worst,
                               tolerance=mc_tolerance,complete_metrics=int(np.isfinite(monitored).sum())))
        if worst<=mc_tolerance or actual_reps>=max_inner_reps or bootstrap<10:
            break
        extra = min(actual_reps,max_inner_reps-actual_reps)
        for b in range(bootstrap):
            addition = _inner(frame,counts,current,following,masks,row_weights[b],extra,seed,1,b,len(chunk_sizes))
            for key in all_inner[b]:
                all_inner[b][key] = np.concatenate([all_inner[b][key],addition[key]],axis=2)
        actual_reps += extra;chunk_sizes.append(extra)
        if progress:
            print(f"uniform MC refinement to R={actual_reps} (tolerance={mc_tolerance})",flush=True)
    point_mc = _mc_variance(point_inner["coefficients"],base,bp,names)
    summaries = {k:ci_summary(values[:,j],point[k]) for j,k in enumerate(names)}
    for name,summary in summaries.items():
        if name.endswith("/c_over_bP"):
            summary["denominator"] = dict(name=name.split('/')[0]+"/bP",**summaries[name.split('/')[0]+"/bP"])
    bands = {}
    for subset in SUBSETS:
        for metric in ("c","delta_a"):
            family = [f"{subset}/q{q:g}/{metric}" for q in QS]
            ids = [names.index(k) for k in family]
            bands[f"{subset}/{metric}"] = dict(metrics=family,**_simultaneous(values[:,ids],[point[k] for k in family]))
    # Complete-case joint covariance, never an undocumented pairwise covariance.
    joint_names = [k for k in names if k.endswith("/c") or k.endswith("/delta_a")]
    joint_ids = [names.index(k) for k in joint_names]
    complete = np.isfinite(values[:,joint_ids]).all(axis=1)
    covariance = np.cov(values[complete][:,joint_ids],rowvar=False) if complete.sum()>1 else np.full((len(joint_ids),len(joint_ids)),np.nan)
    known,known_arrays = known_count_check(seed,min(point_reps,100))
    full_budget = bootstrap>=1000 and point_reps>=400 and inner_reps>=50
    adequate = bool(worst<=mc_tolerance)
    result = dict(metadata=dict(version=VERSION,analysis="thinning",seed=seed,bootstrap=bootstrap,
                    point_reps=point_reps,inner_reps=actual_reps,requested_inner_reps=inner_reps,
                    status="full_budget" if full_budget else "pilot_small_budget",
                    sample_rows=len(frame),sample_sha256=sample_hash(frame),source_code_sha256=file_hash(__file__),
                    stats_code_sha256=file_hash(Path(__file__).with_name("revision_stats.py")),
                    config=dict(qs=QS,block_length=block_length,mc_tolerance=mc_tolerance,max_inner_reps=max_inner_reps,
                                chunk_sizes=chunk_sizes,subsets={s:int(m.sum()) for s,m in zip(SUBSETS,masks)},
                                estimate="mean_of_individually_centered_and_fitted_realizations",
                                eligibility="frozen_before_thinning",regressor="original_C_t"),
                    random_streams=dict(point="SeedSequence([seed,0,0,q_index,0])",
                                        bootstrap_thinning="SeedSequence([seed,1,bootstrap_index,q_index,chunk_index])",
                                        bootstrap_weights="SeedSequence([seed,2,bootstrap_index])",
                                        count_keys="Only unique current/following keys required by positive-weight rows")),
                estimates=summaries,simultaneous_bands=bands,
                poisson_reference=dict(formula="weighted_mean(kappa(g)^2*(1+(Ct+0.5)/(Cnext+0.5)))",g="(Cnext+0.5)/(Ct+0.5)",
                    kappa="(mean/7)/(1+(sd^2/mean)*log(g)/7)",mean=3.6,sd=1.6,
                    note="Approximate plug-in elasticity reference; RA/xbar is approximated by one and half-count ratios stabilize inverses. Not exact nonlinear normalized Poisson variance. Recomputed per draw."),
                point_baseline_fits=base,known_count_model=known,
                covariance=dict(file="thinning_draws.npz",array="joint_covariance",metric_names=joint_names,
                                complete_draws=int(complete.sum()),scope="complete_draws_only"),
                monte_carlo=dict(tolerance=mc_tolerance,adequate=adequate,adaptation=adaptation,
                    rule="Uniform doubling of all bootstrap inner counts, up to cap, based on MC SE / estimated sampling SE; no significance rule. B<10 only diagnoses.",
                    status="tolerance_met" if adequate else "tolerance_not_met_intervals_include_material_MC_error",
                    metrics={k:dict(point_mc_se=np.sqrt(point_mc[j]),bootstrap_mean_mc_se=np.sqrt(inner_var[j]),
                                    estimated_sampling_se=np.sqrt(sampling_var[j]),mc_to_sampling_se=ratios[j]) for j,k in enumerate(names)}),
                failure_policy="Any failed inner fit makes that q average undefined; no silent averaging over successful fits.")
    converged = np.array([r["converged"] for r in all_inner])
    result["fit_diagnostics"] = dict(method_codes=METHOD_CODES,failure_codes=FAILURE_CODES,
                                    point_attempted=int(point_inner["converged"].size),point_converged=int(point_inner["converged"].sum()),
                                    bootstrap_attempted=int(converged.size),bootstrap_converged=int(converged.sum()),
                                    failed_bootstrap_averages=int(np.sum(~converged.all(axis=-1))))
    arrays = dict(values=values,metric_names=np.array(names),point_values=np.array([point[k] for k in names]),
                  point_inner_coefficients=point_inner["coefficients"],point_inner_converged=point_inner["converged"],
                  point_inner_iterations=point_inner["iterations"],point_inner_min_mu=point_inner["min_mu"],
                  point_inner_failures=point_inner["failures"],point_inner_methods=point_inner["methods"],
                  state_multiplicities=sm_all,week_multiplicities=wm_all,
                  state_ids=np.unique(frame.state).astype(str),week_ids=np.unique(frame.week_end).astype(str),
                  count_states=np.array([str(k[0]) for k in keys]),count_weeks=np.array([str(k[1]) for k in keys]),
                  counts=counts,current_indices=current,following_indices=following,
                  row_state_ids=frame.state.to_numpy().astype(str),row_week_ids=frame.week_end.to_numpy().astype(str),
                  baseline_coefficients=np.array([[[np.nan if f[k] is None else f[k] for k in ("a","b")] for f in ff] for ff in baseline_fits]),
                  bootstrap_bP=np.array(references),inner_mc_variance=mc_var,
                  joint_covariance=covariance,covariance_metric_names=np.array(joint_names),**known_arrays)
    for key in all_inner[0]:
        arrays["bootstrap_inner_"+key] = np.array([r[key] for r in all_inner])
    return result,arrays


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path,default=Path(__file__).resolve().parents[1]/"analysis/revised")
    parser.add_argument("--bootstrap",type=int,default=1000)
    parser.add_argument("--point-reps",type=int,default=400)
    parser.add_argument("--inner-reps",type=int,default=50)
    parser.add_argument("--seed",type=int,default=SEED)
    parser.add_argument("--block-length",type=int,default=1)
    parser.add_argument("--mc-tolerance",type=float,default=.10)
    parser.add_argument("--max-inner-reps",type=int,default=200)
    args = parser.parse_args(argv)
    from revision_data import load_sample
    frame = load_sample(stat="RA",gi_mean=3.6)
    result,arrays = run_analysis(frame,args.bootstrap,args.point_reps,args.inner_reps,args.seed,args.block_length,
                                args.mc_tolerance,args.max_inner_reps,progress=True)
    draw_path = args.output_dir/"thinning_draws.npz"
    save_draws(draw_path,**arrays)
    result["metadata"]["draws"] = dict(file=draw_path.name,sha256=file_hash(draw_path))
    write_json(args.output_dir/"thinning.json",result)
    print(f"thinning: B={args.bootstrap}, point R={args.point_reps}, inner R={result['metadata']['inner_reps']}, {result['metadata']['status']}; {args.output_dir}")
    return result


if __name__ == "__main__":
    main()
