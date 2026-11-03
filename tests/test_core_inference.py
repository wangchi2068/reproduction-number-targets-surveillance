"""Numerical regression tests for corrective inference; no production outputs."""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

CODE = Path(__file__).resolve().parents[1]/"experiment/code"
sys.path.insert(0,str(CODE))

from revision_stats import (bootstrap_weights,ci_summary,fit_affine,fit_loss,
                            fit_loss_batch,relative_loss,save_draws)
from s03c_state_fe import (fit_design_diagnostic,fit_group_affine,group_wls,
                           multiplicity_median,state_summaries)
from s03b_thinning import (QS,_average,_metrics,count_vector,poisson_reference,
                           thinning_realizations)


def toy_frame():
    rows = []
    for si,state in enumerate(("AL","NY","WY")):
        for t,week in enumerate(pd.date_range("2020-10-03",periods=7,freq="7D")):
            c = 55+12*si+8*t
            rows.append((state,week,"2020",c,c+8,c-3,1+.02*si+.03*t))
    return pd.DataFrame(rows,columns=["state","week_end","season","C_t","C_next","D_t","x"])


def test_affine_signed_intercept_and_degenerate_failure():
    D = np.arange(1.,15.)
    f = fit_affine(-.02+1/D,D)
    assert f["converged"] and f["min_mu"]>0
    assert f["a"] == pytest.approx(-.02,abs=1e-10)
    assert f["b"] == pytest.approx(1,abs=1e-10)
    assert f["ratio"] is None
    for y,d in ((np.ones(4),np.ones(4)),(np.zeros(4),np.arange(1.,5.))):
        f = fit_affine(y,d)
        assert not f["converged"] and f["a"] is None and f["b"] is None
        assert f["failure"]


def test_relative_loss_and_fit_match_explicit_duplication():
    rng = np.random.default_rng(52)
    x = rng.lognormal(.1,.15,80)
    D = rng.uniform(10,500,80)
    weights = rng.integers(0,5,80)
    duplicate = np.repeat(np.arange(80),weights)
    loss,xbar = relative_loss(x,weights)
    repeated,repeat_mean = relative_loss(x[duplicate])
    np.testing.assert_allclose(loss[duplicate],repeated,rtol=1e-12,atol=1e-15)
    assert xbar == pytest.approx(repeat_mean)
    f,g = fit_loss(x,D,weights),fit_loss(x[duplicate],D[duplicate])
    assert f["converged"] and g["converged"]
    np.testing.assert_allclose([f["a"],f["b"]],[g["a"],g["b"]],rtol=1e-7,atol=1e-9)


def test_group_closed_form_matches_direct_design_weighted_ls():
    rng = np.random.default_rng(63)
    codes = np.repeat(np.arange(4),20)
    z,y,w = rng.uniform(size=80),rng.gamma(2,.5,80),rng.uniform(.1,3,80)
    X = np.column_stack([np.eye(4)[codes],z])
    root = np.sqrt(w)
    direct = np.linalg.lstsq(X*root[:,None],y*root,rcond=None)[0]
    np.testing.assert_allclose(group_wls(y,z,codes,w,4),direct,rtol=1e-11,atol=1e-11)


def test_group_fgls_matches_direct_design_and_row_duplication():
    rng = np.random.default_rng(64)
    states = np.repeat(["AK","DC","NY","WY"],25)
    _,codes = np.unique(states,return_inverse=True)
    D = rng.uniform(30,500,len(states))
    y = (.02+codes*.004+1.3/D)*rng.gamma(2,.5,len(states))
    w = rng.integers(1,4,len(states))
    group = fit_group_affine(y,D,states,w)
    X = np.column_stack([np.eye(4)[codes],1/D])
    direct = fit_design_diagnostic(y,X,w)
    duplicate = np.repeat(np.arange(len(y)),w)
    repeated = fit_group_affine(y[duplicate],D[duplicate],states[duplicate])
    assert group["converged"] and direct["converged"] and repeated["converged"]
    np.testing.assert_allclose(np.r_[group["intercepts"],group["b"]],direct["coefficients"],rtol=3e-6,atol=1e-8)
    np.testing.assert_allclose(group["intercepts"],repeated["intercepts"],rtol=3e-7,atol=1e-9)
    assert group["b"] == pytest.approx(repeated["b"],rel=1e-7)


def test_missing_states_not_zero_and_summaries_count_occurrences():
    e = toy_frame()
    names,codes = np.unique(e.state,return_inverse=True)
    sm = np.array([2,0,1])
    wm = np.array([1,2,0,1,2,1,0])
    _,wc = np.unique(e.week_end,return_inverse=True)
    w = sm[codes]*wm[wc]
    loss,_ = relative_loss(e.x.to_numpy(),w)
    fit = fit_group_affine(loss,e.C_t,e.state,w)
    assert fit["converged"]
    assert np.isnan(fit["intercepts"][1]) and not fit["active"][1]
    result = state_summaries(fit,e,w,sm)
    occurrences = np.repeat(fit["intercepts"],sm)
    assert result["mean"] == pytest.approx(occurrences.mean())
    assert result["median"] == pytest.approx(np.median(occurrences))
    duplicate = np.repeat(np.arange(len(e)),w.astype(int))
    assert result["row"] == pytest.approx(np.mean(fit["intercepts"][codes[duplicate]]))
    assert result["case"] == pytest.approx(np.average(fit["intercepts"][codes[duplicate]],weights=e.C_t.to_numpy()[duplicate]))
    assert multiplicity_median([1,9,12],[1,1,0]) == 5


def test_sorted_bootstrap_multiplicity_mapping_and_season_blocks():
    e = toy_frame().iloc[::-1].reset_index(drop=True)
    states,sc = np.unique(e.state,return_inverse=True)
    weeks,wc = np.unique(e.week_end,return_inverse=True)
    w,sm,wm = bootstrap_weights(e,np.random.default_rng(99))
    np.testing.assert_array_equal(w,sm[sc]*wm[wc])
    assert len(sm)==len(states) and len(wm)==len(weeks)
    assert sm.sum()==len(states) and wm.sum()==len(weeks)
    later = e.copy()
    later.week_end += pd.Timedelta(days=364)
    later.season = "2021"
    e = pd.concat([e,later],ignore_index=True)
    weeks = np.sort(e.week_end.unique())
    for length in (2,3):
        _,sm,wm = bootstrap_weights(e,np.random.default_rng(99),length)
        assert sm.sum()==3
        assert wm[:7].sum()==7 and wm[7:].sum()==7


def test_signed_ci_does_not_positive_filter():
    s = ci_summary([-2,-1,1,2,None],1)
    assert s["ci95"][0]<0<s["ci95"][1]
    assert s["finite_fraction"] == .8
    assert s["nonpositive_fraction"] == .4
    assert s["interval_scope"] == "conditional_on_finite_draws"


def test_thinning_reuses_overlap_and_only_needed_counts():
    e = toy_frame()
    keys,counts,current,following = count_vector(e)
    assert len(keys)==3*8
    ct,cn = thinning_realizations(counts,current,following,.4,20,np.random.default_rng(10))
    for row in range(len(e)-1):
        if e.state.iloc[row]==e.state.iloc[row+1]:
            np.testing.assert_array_equal(cn[:,row],ct[:,row+1])
    active = np.zeros(len(e),bool);active[0]=True
    class CountingRNG:
        def binomial(self,n,p,size):
            assert len(n)==2 and size==(3,2)
            return np.broadcast_to(n,size).copy()
    ct,cn = thinning_realizations(counts,current,following,.4,3,CountingRNG(),active)
    assert np.isnan(ct[:,1:]).all() and np.isfinite(cn[:,0]).all()
    bad = e.copy();bad.loc[1,"C_t"] += 1
    with pytest.raises(ValueError,match="Conflicting"):
        count_vector(bad)


def test_poisson_reference_derivative_and_weighting():
    e = toy_frame()
    w = np.arange(1,len(e)+1)
    g = (e.C_next.to_numpy()+.5)/(e.C_t.to_numpy()+.5)
    kappa = (3.6/7)/(1+(1.6**2/3.6)*np.log(g)/7)
    from revision_data import wl_R
    eps = 1e-5
    derivative = (np.log(wl_R(g*np.exp(eps)))-np.log(wl_R(g*np.exp(-eps))))/(2*eps)
    np.testing.assert_allclose(kappa,derivative,rtol=1e-8)
    assert poisson_reference(e,w) == pytest.approx(np.average(kappa**2*(1+1/g),weights=w))
    assert poisson_reference(e,w) != pytest.approx(poisson_reference(e),rel=1e-6)


def test_batch_individual_fits_equal_scalar_not_average_loss():
    rng = np.random.default_rng(77)
    D = rng.uniform(5,400,90)
    x = rng.lognormal(.2,np.linspace(.15,.5,8)[:,None],size=(8,90))
    w = rng.integers(0,4,len(D))
    batch = fit_loss_batch(x,D,w)
    single = [fit_loss(row,D,w) for row in x]
    assert batch["converged"].all() and all(f["converged"] for f in single)
    expected = np.array([[f["a"],f["b"]] for f in single])
    np.testing.assert_allclose(np.column_stack([batch["a"],batch["b"]]),expected,rtol=2e-6,atol=1e-8)
    averaged_losses = np.mean([relative_loss(row,w)[0] for row in x],axis=0)
    wrong = fit_affine(averaged_losses,D,w)
    assert not np.allclose(expected.mean(axis=0),[wrong["a"],wrong["b"]],rtol=1e-3)
    coefficients = np.zeros((3,7,2,2));coefficients[0,0,0,0]=np.nan
    assert np.isnan(_average(coefficients)[0,0,0])


def test_q_exclusions_and_intercept_changes_are_paired():
    inv = 1/QS-1
    coefficients = np.zeros((3,7,2))
    coefficients[:,:,0] = .02+.001*inv
    coefficients[:,:,1] = .4+.5*inv
    baseline = [dict(a=.02,b=.4)]*3
    metrics = _metrics(coefficients,baseline,np.ones(3)*.45)
    for subset in ("full","highCt50","both50"):
        assert metrics[f"{subset}/full/c"] == pytest.approx(.5)
        assert metrics[f"{subset}/noq01/c"] == pytest.approx(.5)
        assert metrics[f"{subset}/noq01/c_change_vs_full"] == pytest.approx(0)
        assert metrics[f"{subset}/full/delta_a_rate"] == pytest.approx(.001)


def test_deterministic_draw_files(tmp_path):
    a,b = tmp_path/"one.npz",tmp_path/"two.npz"
    arrays = dict(x=np.arange(7),labels=np.array(["A","B"]))
    save_draws(a,**arrays);save_draws(b,**arrays)
    assert a.read_bytes()==b.read_bytes()
    with np.load(a,allow_pickle=False) as data:
        np.testing.assert_array_equal(data["x"],arrays["x"])


def test_imports_do_not_open_outputs_or_run_analysis(tmp_path):
    script = r'''
import builtins, pathlib, sys
sys.dont_write_bytecode=True
sys.path.insert(0,sys.argv[1])
original=builtins.open
path_open=pathlib.Path.open
def checked(file, mode='r', *args, **kwargs):
    if any(c in mode for c in 'wax+'):
        raise AssertionError('import attempted write: '+str(file))
    return original(file,mode,*args,**kwargs)
def checked_path(self,mode='r',*args,**kwargs):
    if any(c in mode for c in 'wax+'):
        raise AssertionError('import attempted write: '+str(self))
    return path_open(self,mode,*args,**kwargs)
builtins.open=checked
pathlib.Path.open=checked_path
pathlib.Path.mkdir=lambda *a,**k: (_ for _ in ()).throw(AssertionError('import mkdir'))
import revision_stats, s03b_thinning, s03c_state_fe
'''
    subprocess.run([sys.executable,"-B","-c",script,str(CODE)],cwd=tmp_path,check=True,
                   env=dict(os.environ,PYTHONDONTWRITEBYTECODE="1",OPENBLAS_NUM_THREADS="1"),capture_output=True,text=True)
    assert not list(tmp_path.iterdir())


def test_actual_D_t_sensitivity_uses_distinct_denominator():
    from revision_data import load_sample
    e = load_sample()
    loss,_ = relative_loss(e.x.to_numpy())
    a,b = fit_group_affine(loss,e.C_t,e.state),fit_group_affine(loss,e.D_t,e.state)
    assert a["converged"] and b["converged"]
    assert np.sum(e.D_t.to_numpy()!=e.C_t.to_numpy())==1742
    assert a["b"] == pytest.approx(.3793764702,rel=1e-6)
    assert b["b"] == pytest.approx(.369082,rel=5e-6)
    assert abs(a["b"]-b["b"])>.005


def test_gamma_objective_fallback_converges_without_positive_intercept_constraint(monkeypatch):
    import revision_stats as stats
    monkeypatch.setattr(stats,"_iterate_positive",lambda *args,**kwargs:
                        (np.array([1.,0.]),dict(converged=False,iterations=200,failure="forced_test")))
    D = np.linspace(1.,15.,100)
    f = stats.fit_affine(-.02+1/D,D)
    assert f["converged"] and f["method"]=="gamma_objective_endpoint_fallback"
    assert f["min_mu"]>0 and f["ratio"] is None
    np.testing.assert_allclose([f["a"],f["b"]],[-.02,1],rtol=1e-6,atol=1e-8)


def test_failed_diagnostics_have_no_point_estimate():
    X = np.ones((8,3))
    f = fit_design_diagnostic(np.arange(1.,9.),X)
    assert not f["converged"] and f["b"] is None
    assert f["failure"]=="rank_deficient_design"
