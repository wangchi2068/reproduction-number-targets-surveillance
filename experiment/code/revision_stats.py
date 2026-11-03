"""Import-safe inference for the corrective revision (no historical result reads).

FGLS solves the gamma quasi-score with multiplicity weights and positive fitted
means. Intercepts and slopes themselves are unconstrained. Nonconvergence is a
failed fit, never an estimate. Percentile intervals retain signed finite draws;
undefined ratios must be accompanied by their denominator sign summaries.
"""
from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np
import pandas as pd
from scipy.optimize import minimize

VERSION = "revision-stats-1.0"


def _weights(n, w=None):
    w = np.ones(n, dtype=float) if w is None else np.asarray(w, dtype=float)
    if w.shape != (n,) or not np.all(np.isfinite(w)) or np.any(w < 0):
        raise ValueError("Weights must be a finite, nonnegative vector matching observations")
    if not np.any(w > 0):
        raise ValueError("No positive-weight observations")
    return w


def relative_loss(x, w=None):
    """Relative squared loss, with the weighted centre refitted in every draw."""
    x = np.asarray(x, dtype=float)
    if x.ndim != 1:
        raise ValueError("x must be a vector")
    w = _weights(len(x), w)
    keep = w > 0
    if not np.all(np.isfinite(x[keep])):
        raise ValueError("Nonfinite positive-weight x")
    xbar = float(np.average(x[keep], weights=w[keep]))
    if not np.isfinite(xbar) or xbar <= 0:
        raise ValueError("Relative loss requires a positive finite weighted mean")
    with np.errstate(invalid="ignore", over="ignore"):
        loss = ((x - xbar) / xbar) ** 2
    return loss, xbar


def _failure(reason, iterations=0, method="fgls", **extra):
    return dict(a=None, b=None, ratio=None, converged=False,
                iterations=int(iterations), method=method, min_mu=None,
                failure=str(reason), **extra)


def _objective(y, mu, w):
    if np.any(mu <= 0) or not np.all(np.isfinite(mu)):
        return np.inf
    return float(np.sum(w * (y / mu + np.log(mu))))


def _iterate_positive(y, w, predict, solve, start, maxit=200, tol=1e-9):
    """Damped Fisher scoring; solve is exact weighted LS in the model space.

    y and w are normalized by their weighted mean and sum by the caller. Each
    accepted step decreases the *same* gamma objective and preserves mu>0.
    A small undamped scoring step establishes convergence, not a small damped
    step at the positive boundary. Returned unsuccessful candidates are private.
    """
    beta = np.asarray(start, dtype=float)
    mu = predict(beta)
    floor = 1e-12
    objective = _objective(y, mu, w)
    if not np.isfinite(objective):
        return beta, dict(converged=False, iterations=0, failure="invalid_initial_mean")
    for iteration in range(1, maxit + 1):
        precision = w * (np.min(mu) / mu) ** 2
        try:
            proposal = solve(precision)
        except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exc:
            return beta, dict(converged=False, iterations=iteration, failure=str(exc))
        delta = proposal - beta
        change = predict(proposal) - mu
        residual = float(np.max(np.abs(change) / np.maximum(1.0, np.abs(mu))))
        if residual <= tol and np.min(predict(proposal)) > floor:
            return proposal, dict(converged=True, iterations=iteration,
                                  objective=_objective(y, predict(proposal), w),
                                  scoring_residual=residual, failure=None)
        step = 1.0
        negative = change < 0
        if np.any(negative):
            step = min(step, .99 * float(np.min((mu[negative] - floor) / -change[negative])))
        accepted = False
        directional = float(np.sum(w * (mu - y) / mu ** 2 * change))
        for _ in range(45):
            if step <= 0:
                break
            candidate = beta + step * delta
            next_mu = predict(candidate)
            next_obj = _objective(y, next_mu, w)
            if np.min(next_mu) > floor and next_obj <= objective + 1e-4 * step * min(directional, 0) + 2e-15:
                beta, mu, objective = candidate, next_mu, next_obj
                accepted = True
                break
            step *= .5
        if not accepted:
            return beta, dict(converged=False, iterations=iteration, failure="line_search_failed",
                              scoring_residual=residual)
    return beta, dict(converged=False, iterations=maxit, failure="iteration_limit",
                      scoring_residual=residual)


def _wls_line(y, z, w):
    """Centred closed-form weighted LS; avoids the ill-scaled 1/D normal matrix."""
    total = float(np.sum(w))
    if total <= 0:
        raise ValueError("zero_precision")
    zbar = float(w @ z / total)
    ybar = float(w @ y / total)
    dz = z - zbar
    variance = float(w @ (dz * dz))
    if variance <= np.finfo(float).eps * total * 1e-6:
        raise ValueError("rank_deficient_slope")
    slope = float(w @ (dz * (y - ybar)) / variance)
    return np.array([ybar - slope * zbar, slope])


def fit_affine(loss, D, w=None):
    """Fit E(loss)=a+b/D; signed a is allowed, ratio is None unless a>0.

    Invalid data/rank/optimization produce a diagnostic dict with null a,b.
    The gamma-objective fallback uses log endpoint means, so it constrains
    fitted means at the observed 1/D extremes without constraining a or b.
    """
    try:
        y = np.asarray(loss, dtype=float)
        D = np.asarray(D, dtype=float)
        if y.ndim != 1 or D.shape != y.shape:
            raise ValueError("loss and D must be equally sized vectors")
        weights = _weights(len(y), w)
        keep = weights > 0
        y, D, weights = y[keep], D[keep], weights[keep]
        if not np.all(np.isfinite(y)) or np.any(y < 0) or not np.all(np.isfinite(D)) or np.any(D <= 0):
            raise ValueError("Require finite nonnegative loss and positive D on active rows")
        if len(y) < 2:
            raise ValueError("insufficient_observations")
        z = 1.0 / D
        lo, span = float(z.min()), float(np.ptp(z))
        if span <= np.finfo(float).eps * max(1.0, float(z.max())):
            raise ValueError("rank_deficient_slope")
        scale = float(np.average(y, weights=weights))
        if scale <= 0:
            raise ValueError("zero_loss_has_no_positive_gamma_mean_solution")
        y, weights = y / scale, weights / weights.sum()
        t = (z - lo) / span
        predict = lambda beta: beta[0] + beta[1] * t
        solve = lambda precision: _wls_line(y, t, precision)
        beta, diagnostics = _iterate_positive(y, weights, predict, solve, [1.0, 0.0])
        method = "damped_fgls"
        if not diagnostics["converged"]:
            previous = diagnostics.copy()
            endpoints = np.maximum([beta[0], beta.sum()], 1e-10)
            def objective_gradient(eta):
                means = np.exp(eta)
                mu = (1 - t) * means[0] + t * means[1]
                score = weights * (mu - y) / mu ** 2
                return _objective(y, mu, weights), np.array([
                    score @ (1 - t) * means[0], score @ t * means[1]])
            opt = minimize(objective_gradient, np.log(endpoints), jac=True,
                           method="L-BFGS-B", bounds=[(-27, 27)] * 2,
                           options={"maxiter": 1000, "ftol": 1e-14, "gtol": 1e-9, "maxls": 50})
            means = np.exp(opt.x)
            beta = np.array([means[0], means[1] - means[0]])
            objective, gradient = objective_gradient(opt.x)
            good = bool(opt.success and np.max(np.abs(gradient)) < 1e-6 and np.min(predict(beta)) > 1e-11)
            diagnostics = dict(converged=good, iterations=previous["iterations"] + int(opt.nit),
                               objective=objective, scoring_residual=float(np.max(np.abs(gradient))),
                               fgls_failure=previous, optimizer_message=str(opt.message),
                               failure=None if good else "gamma_objective_fallback_failed")
            method = "gamma_objective_endpoint_fallback"
        if not diagnostics["converged"]:
            return _failure(diagnostics.pop("failure"), method=method, **{k:v for k,v in diagnostics.items() if k != "converged"})
        b = float(scale * beta[1] / span)
        a = float(scale * beta[0] - b * lo)
        mu = predict(beta) * scale
        return dict(a=a, b=b, ratio=b/a if a > 0 else None, method=method,
                    min_mu=float(mu.min()), rank=2, n_active=int(len(y)), **diagnostics)
    except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
        return _failure(exc)


def fit_loss(x, D, w=None):
    """Centre x with draw weights, then call fit_affine."""
    try:
        loss, xbar = relative_loss(x, w)
    except (ValueError, FloatingPointError) as exc:
        return dict(_failure(exc), xbar=None)
    return dict(fit_affine(loss, D, w), xbar=xbar)


def fit_loss_batch(x, D, w=None):
    """Vectorized *individual* affine fits; never fit the average loss.

    Returns arrays of per-realization coefficients and diagnostics. Difficult
    rows fall back to fit_affine, preserving exactly the scalar objective.
    """
    x, D = np.asarray(x, float), np.asarray(D, float)
    if x.ndim != 2 or D.shape != (x.shape[1],):
        raise ValueError("x must be realizations by observations, matching D")
    nr = len(x)
    weights = _weights(len(D), w)
    active = weights > 0
    xx, dd, ww = x[:,active], D[active], weights[active]
    xbar = np.sum(xx*ww, axis=1)/ww.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        losses = ((xx-xbar[:,None])/xbar[:,None])**2
    a, b, mins = np.full(nr,np.nan), np.full(nr,np.nan), np.full(nr,np.nan)
    converged, iterations = np.zeros(nr,bool), np.zeros(nr,int)
    methods = np.full(nr,"failed",dtype="U48")
    reasons = np.full(nr,"",dtype="U100")
    valid = np.isfinite(losses).all(axis=1) & (xbar > 0)
    scale = np.sum(losses*ww,axis=1)/ww.sum()
    valid &= scale > 0
    z = 1/dd
    lo, span = float(z.min()), float(np.ptp(z))
    if span <= np.finfo(float).eps or np.any(dd <= 0):
        valid[:] = False
    if np.any(valid):
        t = (z-lo)/span
        ww = ww/ww.sum()
        y = np.divide(losses, scale[:,None], out=np.ones_like(losses), where=scale[:,None]>0)
        beta = np.tile([1.,0.],(nr,1))
        pending = valid.copy()
        for iteration in range(1,101):
            ids = np.flatnonzero(pending)
            if not len(ids):
                break
            yy, bb = y[ids], beta[ids]
            mu = bb[:,0,None]+bb[:,1,None]*t
            precision = ww*(mu.min(axis=1)[:,None]/mu)**2
            total = precision.sum(axis=1)
            zbar = np.sum(precision*t,axis=1)/total
            ybar = np.sum(precision*yy,axis=1)/total
            dz = t-zbar[:,None]
            variance = np.sum(precision*dz**2,axis=1)
            with np.errstate(divide="ignore",invalid="ignore"):
                slope = np.sum(precision*dz*(yy-ybar[:,None]),axis=1)/variance
            prop = np.column_stack([ybar-slope*zbar,slope])
            change = (prop[:,0,None]+prop[:,1,None]*t)-mu
            residual = np.max(np.abs(change)/np.maximum(1.,np.abs(mu)),axis=1)
            done = (residual<=1e-9) & (np.min(mu+change,axis=1)>1e-12)
            if done.any():
                finish = ids[done]
                beta[finish] = prop[done]
                converged[finish], iterations[finish], methods[finish] = True, iteration, "batch_damped_fgls"
                pending[finish] = False
            todo = ~done
            if not todo.any():
                continue
            ids, yy, bb, mu, prop, change = ids[todo], yy[todo], bb[todo], mu[todo], prop[todo], change[todo]
            with np.errstate(divide="ignore",invalid="ignore"):
                bound = np.where(change<0,(mu-1e-12)/-change,np.inf)
            step = np.minimum(1.,.99*np.min(bound,axis=1))
            objective = np.sum(ww*(yy/mu+np.log(mu)),axis=1)
            directional = np.sum(ww*(mu-yy)/mu**2*change,axis=1)
            accepted = np.zeros(len(ids),bool)
            for _ in range(40):
                next_mu = mu+step[:,None]*change
                with np.errstate(divide="ignore",invalid="ignore"):
                    next_obj = np.sum(ww*(yy/next_mu+np.log(next_mu)),axis=1)
                good = (~accepted) & (next_mu.min(axis=1)>1e-12) & (next_obj <= objective+1e-4*step*np.minimum(directional,0)+2e-15)
                beta[ids[good]] = bb[good]+step[good,None]*(prop[good]-bb[good])
                accepted |= good
                if accepted.all():
                    break
                step[~accepted] *= .5
            pending[ids[~accepted]] = False  # scalar fallback below
        success = np.flatnonzero(converged)
        b[success] = scale[success]*beta[success,1]/span
        a[success] = scale[success]*beta[success,0]-b[success]*lo
        mins[success] = scale[success]*np.minimum(beta[success,0],beta[success].sum(axis=1))
    for i in np.flatnonzero(~converged):
        f = fit_affine(losses[i], dd, weights[active]) if np.isfinite(losses[i]).all() else _failure("invalid_relative_loss")
        converged[i], iterations[i], methods[i] = f["converged"], f["iterations"], f["method"]
        reasons[i] = f.get("failure") or ""
        if f["converged"]:
            a[i], b[i], mins[i] = f["a"], f["b"], f["min_mu"]
    return dict(a=a,b=b,xbar=xbar,converged=converged,iterations=iterations,
                min_mu=mins,method=methods,failure=reasons)


def bootstrap_weights(frame, rng, block_length=1):
    """State x calendar-week pigeonhole multiplicities in sorted unique order.

    L=1 samples all observed weeks independently. L>1 samples circular moving
    blocks separately within season and contiguous weekly runs (never crosses
    a season boundary or missing-week gap), retaining each run's week count.
    No rows are duplicated and original state IDs are never relabelled.
    """
    if int(block_length) != block_length or block_length < 1:
        raise ValueError("block_length must be a positive integer")
    if len(frame) == 0 or frame[["state", "week_end"]].isna().any().any():
        raise ValueError("Nonempty, nonmissing state/week keys required")
    states, si = np.unique(frame.state.to_numpy(), return_inverse=True)
    dates = pd.to_datetime(frame.week_end).to_numpy()
    weeks, wi = np.unique(dates, return_inverse=True)
    sm = np.bincount(rng.integers(len(states), size=len(states)), minlength=len(states))
    if block_length == 1:
        wm = np.bincount(rng.integers(len(weeks), size=len(weeks)), minlength=len(weeks))
    else:
        if "season" not in frame or frame.season.isna().any():
            raise ValueError("Season labels required for time blocks")
        mapping = frame.assign(_week=dates).groupby("_week").season.agg(lambda s: tuple(s.unique()))
        if any(len(v) != 1 for v in mapping):
            raise ValueError("A calendar week cannot belong to multiple seasons")
        seasons = np.array([mapping.loc[pd.Timestamp(d)][0] for d in weeks])
        wm = np.zeros(len(weeks), dtype=int)
        for season in np.unique(seasons):
            ids = np.flatnonzero(seasons == season)
            cuts = np.flatnonzero(np.diff(weeks[ids]) != np.timedelta64(7, "D")) + 1
            for run in np.split(ids, cuts):
                n = len(run)
                starts = rng.integers(n, size=int(np.ceil(n / block_length)))
                sampled = ((starts[:, None] + np.arange(block_length)) % n).ravel()[:n]
                wm += np.bincount(run[sampled], minlength=len(weeks))
    return (sm[si] * wm[wi]).astype(float), sm, wm


def ci_summary(values, point):
    """Signed percentile CI; missing/undefined draws remain in the denominator.

    nonpositive_fraction is P(draw<=0) over *all requested* draws. It is not a
    denominator-sign diagnostic: callers with ratios must save denominators
    separately. Intervals with missing draws are explicitly conditional on a
    finite estimate and must not be presented as unconditional ratio coverage.
    """
    values = np.asarray([np.nan if v is None else v for v in values], dtype=float)
    finite = np.isfinite(values)
    n = len(values)
    interval = np.quantile(values[finite], [.025, .975]).tolist() if finite.sum() >= 2 else [None, None]
    return dict(value=float(point) if point is not None and np.isfinite(point) else None,
                ci95=interval, finite_fraction=float(finite.mean()) if n else 0.0,
                nonpositive_fraction=float(np.sum(finite & (values <= 0))/n) if n else 0.0,
                n_draws=n, n_finite=int(finite.sum()),
                interval_scope="all_draws" if n and finite.all() else "conditional_on_finite_draws")


def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [json_safe(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (Path, pd.Timestamp)):
        return str(value)
    return value


def write_json(path, result):
    """Explicit output action, called only from runners/main."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json_safe(result), indent=2, allow_nan=False) + "\n", encoding="utf-8")


def save_draws(path, **arrays):
    """Deterministic NPZ bytes (sorted members, fixed ZIP timestamps, no pickle)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name in sorted(arrays):
            buffer = io.BytesIO()
            np.lib.format.write_array(buffer, np.asarray(arrays[name]), allow_pickle=False)
            entry = zipfile.ZipInfo(name + ".npy", date_time=(1980, 1, 1, 0, 0, 0))
            entry.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(entry, buffer.getvalue())


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sample_hash(frame):
    cols = [c for c in ("state", "week_end", "season", "C_t", "C_next", "D_t", "x", "D") if c in frame]
    return hashlib.sha256(frame[cols].to_csv(index=False, float_format="%.17g").encode()).hexdigest()
