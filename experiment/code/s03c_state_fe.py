"""Corrected paired state-intercept inference; import performs no IO.

State-intercept FGLS uses weighted group means and a within-group slope at each
iteration. Missing states have no coefficient. Summaries target sampled state
occurrences, or resampled row/case totals. Calendar-week models are point-only
numerical diagnostics, never primary bootstrap inference.
"""
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import LinearConstraint, minimize

from revision_stats import (VERSION, _failure, _iterate_positive, _objective, _weights,
                            bootstrap_weights, ci_summary, file_hash, fit_loss,
                            relative_loss, sample_hash, save_draws, write_json)

SEED = 20261104


def group_wls(y, z, codes, w, n_groups):
    """Exact WLS for one intercept per active group and one common slope."""
    totals = np.bincount(codes, weights=w, minlength=n_groups)
    if np.any(totals <= 0):
        raise ValueError("empty_active_group")
    ybar = np.bincount(codes, weights=w*y, minlength=n_groups) / totals
    zbar = np.bincount(codes, weights=w*z, minlength=n_groups) / totals
    dz, dy = z-zbar[codes], y-ybar[codes]
    denominator = float(w @ (dz*dz))
    if denominator <= 1e-14 * float(w.sum()):
        raise ValueError("rank_deficient_within_slope")
    slope = float(w @ (dz*dy) / denominator)
    return np.r_[ybar-slope*zbar, slope]


def fit_group_affine(loss, D, groups, w=None):
    """Stable positive-mean FGLS with active original group IDs preserved.

    All input group IDs appear in the output; omitted groups have NaN intercepts
    and active=False, never zero coefficients. Signed intercepts are allowed.
    """
    names, codes = np.unique(np.asarray(groups), return_inverse=True)
    empty = dict(groups=names, intercepts=np.full(len(names), np.nan),
                 active=np.zeros(len(names), dtype=bool))
    try:
        y, D = np.asarray(loss, float), np.asarray(D, float)
        if y.ndim != 1 or y.shape != D.shape or len(y) != len(codes):
            raise ValueError("incompatible_data_shapes")
        weights = _weights(len(y), w)
        keep = weights > 0
        y, z, weights, codes = y[keep], 1/D[keep], weights[keep], codes[keep]
        active_ids, local = np.unique(codes, return_inverse=True)
        ng = len(active_ids)
        if not np.all(np.isfinite(y)) or np.any(y < 0) or not np.all(np.isfinite(z)) or np.any(z <= 0):
            raise ValueError("invalid_active_loss_or_scale")
        scale = float(np.average(y, weights=weights))
        span, lo = float(np.ptp(z)), float(z.min())
        if scale <= 0 or span <= 0:
            raise ValueError("zero_loss_or_rank_deficient_slope")
        y, weights, z = y/scale, weights/weights.sum(), (z-lo)/span
        solve = lambda p: group_wls(y, z, local, p, ng)
        solve(weights)  # rank check before iteration
        predict = lambda b: b[local]+b[-1]*z
        beta, diagnostics = _iterate_positive(y, weights, predict, solve, np.r_[np.ones(ng), 0.], maxit=300)
        method = "group_mean_within_slope_damped_fgls"
        if not diagnostics["converged"]:
            # Same gamma objective; linear constraints enforce all active fitted means.
            prior = diagnostics.copy()
            design = np.column_stack([np.eye(ng)[local], z])
            def objective_gradient(b):
                mu = design @ b
                if np.any(mu <= 0):
                    return 1e100, np.zeros(ng+1)
                return _objective(y, mu, weights), design.T @ (weights*(mu-y)/mu**2)
            opt = minimize(objective_gradient, beta, jac=True, method="SLSQP",
                           constraints=[LinearConstraint(design, 1e-10, np.inf)],
                           options={"maxiter": 600, "ftol": 1e-12})
            mu = predict(opt.x)
            score = objective_gradient(opt.x)[1]
            good = bool(opt.success and np.min(mu) > 1e-9 and np.max(np.abs(score)) < 2e-5)
            beta = opt.x
            diagnostics = dict(converged=good, iterations=prior["iterations"]+int(opt.nit),
                               failure=None if good else "gamma_objective_fallback_failed",
                               optimizer_message=str(opt.message), fgls_failure=prior,
                               scoring_residual=float(np.max(np.abs(score))))
            method = "group_gamma_objective_fallback"
        if not diagnostics["converged"]:
            return dict(_failure(diagnostics.pop("failure"), method=method,
                                 **{k:v for k,v in diagnostics.items() if k != "converged"}), **empty)
        slope = float(scale*beta[-1]/span)
        intercepts = np.full(len(names), np.nan)
        intercepts[active_ids] = scale*beta[:-1]-slope*lo
        active = np.zeros(len(names), dtype=bool)
        active[active_ids] = True
        return dict(b=slope, groups=names, intercepts=intercepts, active=active,
                    min_mu=float(np.min(predict(beta))*scale), method=method,
                    rank=ng+1, n_active=int(keep.sum()), **diagnostics)
    except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
        return dict(_failure(exc), **empty)


def fit_design_diagnostic(loss, design, w=None):
    """Direct-design FGLS for diagnostics; rank/nonconvergence are explicit."""
    try:
        y, X = np.asarray(loss, float), np.asarray(design, float)
        weights = _weights(len(y), w)
        keep = weights > 0
        y, X, weights = y[keep], X[keep], weights[keep]
        rank = int(np.linalg.matrix_rank(X))
        if rank != X.shape[1]:
            return dict(_failure("rank_deficient_design"), rank=rank, columns=int(X.shape[1]))
        scale = float(np.average(y, weights=weights))
        if scale <= 0:
            return _failure("zero_loss")
        y, weights = y/scale, weights/weights.sum()
        def solve(p):
            root = np.sqrt(p)
            beta, _, rank_now, _ = np.linalg.lstsq(X*root[:,None], y*root, rcond=None)
            if rank_now != X.shape[1]:
                raise ValueError("rank_loss_during_iteration")
            return beta
        start = np.linalg.lstsq(X, np.ones(len(y)), rcond=None)[0]
        beta, diagnostics = _iterate_positive(y, weights, lambda b:X@b, solve, start, maxit=250)
        if not diagnostics["converged"]:
            return dict(_failure(diagnostics.pop("failure"), method="direct_design_damped_fgls",
                                 **{k:v for k,v in diagnostics.items() if k != "converged"}), rank=rank)
        return dict(coefficients=beta*scale, b=float(beta[-1]*scale),
                    min_mu=float(np.min(X@beta)*scale), rank=rank,
                    method="direct_design_damped_fgls", **diagnostics)
    except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
        return _failure(exc, method="direct_design_damped_fgls")


def multiplicity_median(values, counts):
    """Median of explicitly duplicated state occurrences, including even-n mean."""
    values, counts = np.asarray(values, float), np.asarray(counts, int)
    keep = counts > 0
    values, counts = values[keep], counts[keep]
    if not len(values) or not np.all(np.isfinite(values)):
        return np.nan
    order = np.argsort(values)
    values, cumulative = values[order], np.cumsum(counts[order])
    n = int(cumulative[-1])
    return float(np.mean(values[np.searchsorted(cumulative, [(n-1)//2+1, n//2+1])]))


def state_summaries(fit, frame, row_weights, state_multiplicities):
    result = {k:np.nan for k in ("mean", "median", "row", "case")}
    if not fit["converged"]:
        return result
    names, codes = np.unique(frame.state.to_numpy(), return_inverse=True)
    if not np.array_equal(names, fit["groups"]):
        raise ValueError("state ID ordering mismatch")
    sm = np.asarray(state_multiplicities, int)
    intercepts = fit["intercepts"]
    if np.any((sm > 0) & ~np.isfinite(intercepts)):
        return result  # selected state with no observed row: the target is not estimable
    present = sm > 0
    result["mean"] = float(np.average(intercepts[present], weights=sm[present]))
    result["median"] = multiplicity_median(intercepts, sm)
    rows = np.bincount(codes, weights=row_weights, minlength=len(names))
    cases = np.bincount(codes, weights=row_weights*frame.C_t.to_numpy(), minlength=len(names))
    for name, totals in (("row", rows), ("case", cases)):
        use = totals > 0
        if use.any():
            result[name] = float(np.average(intercepts[use], weights=totals[use]))
    return result


def _ratio(numerator, denominator):
    if numerator is None or denominator is None or not np.isfinite(numerator) or not np.isfinite(denominator) or denominator <= 0:
        return np.nan
    return float(numerator/denominator)


def _one_draw(frame, row_weights, sm):
    pool = fit_loss(frame.x.to_numpy(), frame.C_t.to_numpy(), row_weights)
    try:
        loss, _ = relative_loss(frame.x.to_numpy(), row_weights)
    except ValueError as exc:
        failed = _failure(exc)
        return {}, {"pooled":pool, "state":failed, "state_D":failed}
    state = fit_group_affine(loss, frame.C_t, frame.state, row_weights)
    state_D = fit_group_affine(loss, frame.D_t, frame.state, row_weights)
    summaries = state_summaries(state, frame, row_weights, sm)
    sums_D = state_summaries(state_D, frame, row_weights, sm)
    vals = dict(pool_a=pool["a"], pool_b=pool["b"], pool_ratio=pool["ratio"],
                state_b=state["b"], state_D_b=state_D["b"], state_D_mean=sums_D["mean"])
    for key, a in summaries.items():
        vals["a_"+key] = a
        vals["r_"+key] = _ratio(state["b"], a)
    vals["state_D_ratio"] = _ratio(state_D["b"], sums_D["mean"])
    vals["paired_ratio"] = _ratio(pool["ratio"], vals["r_mean"])
    return vals, {"pooled":pool, "state":state, "state_D":state_D}


def run_analysis(frame, bootstrap=1000, seed=SEED, block_length=1):
    if bootstrap < 1:
        raise ValueError("bootstrap must be positive")
    states = np.unique(frame.state.to_numpy())
    weeks = np.unique(pd.to_datetime(frame.week_end).to_numpy())
    point, point_fits = _one_draw(frame, np.ones(len(frame)), np.ones(len(states), int))
    metrics = list(point) or ["pool_a", "pool_b", "pool_ratio", "state_b", "state_D_b", "state_D_mean",
                              "a_mean", "r_mean", "a_median", "r_median", "a_row", "r_row", "a_case", "r_case", "state_D_ratio", "paired_ratio"]
    draws = np.full((bootstrap, len(metrics)), np.nan)
    sm_all, wm_all = np.zeros((bootstrap, len(states)), int), np.zeros((bootstrap, len(weeks)), int)
    intercepts = np.full((bootstrap, 2, len(states)), np.nan)
    convergence = np.zeros((bootstrap, 3), bool)
    iterations = np.zeros((bootstrap, 3), int)
    min_mu = np.full((bootstrap, 3), np.nan)
    methods, failures = np.full((bootstrap, 3), "", dtype="U80"), np.full((bootstrap, 3), "", dtype="U160")
    rng = np.random.default_rng(seed)
    for i in range(bootstrap):
        w, sm, wm = bootstrap_weights(frame, rng, block_length)
        sm_all[i], wm_all[i] = sm, wm
        vals, fits = _one_draw(frame, w, sm)
        draws[i] = [np.nan if vals.get(k) is None else vals.get(k, np.nan) for k in metrics]
        for j, key in enumerate(("pooled", "state", "state_D")):
            f = fits[key]
            convergence[i,j], iterations[i,j] = f["converged"], f["iterations"]
            methods[i,j], failures[i,j] = f["method"], f.get("failure") or ""
            min_mu[i,j] = np.nan if f["min_mu"] is None else f["min_mu"]
            if j and "intercepts" in f:
                intercepts[i,j-1] = f["intercepts"]
    summaries = {k:ci_summary(draws[:,j], point.get(k)) for j,k in enumerate(metrics)}
    denominators = {"pool_ratio":"pool_a", "r_mean":"a_mean", "r_median":"a_median",
                    "r_row":"a_row", "r_case":"a_case", "state_D_ratio":"state_D_mean", "paired_ratio":"r_mean"}
    for ratio, denom in denominators.items():
        summaries[ratio]["denominator"] = dict(name=denom, **summaries[denom])
    loss, _ = relative_loss(frame.x.to_numpy())
    week_fit = fit_group_affine(loss, frame.C_t, frame.week_end.astype(str))
    _, si = np.unique(frame.state, return_inverse=True)
    _, wi = np.unique(frame.week_end, return_inverse=True)
    design = np.column_stack([np.eye(len(states))[si], np.eye(len(weeks))[wi,1:], 1/frame.C_t.to_numpy()])
    both_fit = fit_design_diagnostic(loss, design)
    result = dict(metadata=dict(version=VERSION, analysis="fixed_effects", seed=seed, bootstrap=bootstrap,
                                status="full_budget" if bootstrap >= 1000 else "pilot_small_B",
                                config=dict(block_length=block_length, centering="refit_weighted_mean",
                                            state_target="sampled_state_occurrences", diagnostic_scope="point_only"),
                                sample_rows=len(frame), sample_sha256=sample_hash(frame),
                                source_code_sha256=file_hash(__file__), stats_code_sha256=file_hash(Path(__file__).with_name("revision_stats.py"))),
                  estimates=summaries, point_fits=point_fits,
                  diagnostics=dict(week_only=week_fit, state_and_week=both_fit,
                                   note="Point-only numerical diagnostics; failed fits have null slopes, no confidence intervals."),
                  bootstrap_diagnostics={k:dict(converged=int(convergence[:,j].sum()), failed=int((~convergence[:,j]).sum()),
                                               failure_reasons=dict(Counter(failures[:,j][failures[:,j] != ""])),
                                               methods=dict(Counter(methods[:,j]))) for j,k in enumerate(("pooled","state","state_D"))},
                  ratio_note="Ratios require positive denominator. Missing fractions and denominator signs are retained; finite-only intervals are explicitly conditional.",
                  actual_D_t_differences=int(np.sum(frame.D_t.to_numpy() != frame.C_t.to_numpy())))
    arrays = dict(values=draws, metric_names=np.array(metrics), state_ids=states.astype(str), week_ids=weeks.astype(str),
                  state_multiplicities=sm_all, week_multiplicities=wm_all, state_intercepts=intercepts,
                  converged=convergence, iterations=iterations, min_mu=min_mu, methods=methods, failures=failures)
    return result, arrays


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parents[1]/"analysis/revised")
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--point-reps", type=int, default=400, help="Accepted for runner compatibility; state FE is deterministic")
    parser.add_argument("--inner-reps", type=int, default=50, help="Accepted for runner compatibility; state FE has no inner simulation")
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--block-length", type=int, default=1)
    args = parser.parse_args(argv)
    from revision_data import load_sample
    frame = load_sample(stat="RA", gi_mean=3.6)
    result, arrays = run_analysis(frame, args.bootstrap, args.seed, args.block_length)
    result["metadata"]["point_reps"] = 1
    result["metadata"]["inner_reps"] = 0
    draw_path = args.output_dir/"fixed_effects_draws.npz"
    save_draws(draw_path, **arrays)
    result["metadata"]["draws"] = dict(file=draw_path.name, sha256=file_hash(draw_path))
    write_json(args.output_dir/"fixed_effects.json", result)
    print(f"fixed effects: {len(frame)} rows, B={args.bootstrap}, {result['metadata']['status']}; {args.output_dir}")
    return result


if __name__ == "__main__":
    main()
