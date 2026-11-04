"""P0 known-truth audit on the real 1,748 state-week R_A design (mainline A, Section 4).

Synthetic count panels are generated on the frozen R_A sample: the same rows,
the same state-week structure, the same shared counts (the count of week w+1 is
both C_next of week w and C_t of week w+1, as in the authors' thinning code),
the same statistic R_A, the same relative loss and the same affine FGLS fit
(revision_stats.fit_loss / fit_loss_batch). The weekly growth factor g_i of
every row is known, so the target R_i = wl_R(g_i) is known.

Count alignment (revision 2, spec A1). Inside each maximal run of consecutive
weeks of a state the expected count is carried by the true growth factors; the
level of every run is anchored so that the geometric mean of the expected
count over the run equals the geometric mean of the observed counts of that
run. A global scale lambda in [0.5, 1.5] per generator then sets the mean (over
50 panels) of the median retained analyst count D to 45 +/- 1 (observed: 45).
After lambda, the transmission heterogeneity is recalibrated so that the Monte
Carlo mean of the fitted intercept a equals the observed a (0.014955); lambda
is then re-fitted once with the final heterogeneity.

Generators
  S0 poisson              constant g0; independent Poisson counts around the path
  S1 heterogeneity        iid lognormal g_i; Poisson counts
  S2 nb_reporting         constant g0; latent branching chain with negative-
                          binomial offspring (k=0.3) and binomial reporting rho=0.25
  S3 heterogeneity_nb     S1 heterogeneity and S2 offspring/reporting
  S4 structured           state + national-week heterogeneity, Poisson counts
  S5 structured_nb        S4 heterogeneity with the S3 offspring/reporting law
  S2P branching_poisson   S2 with Poisson offspring (pure branching chain)
  S6 reporting_shock      S3 plus a state-week reporting shock rho_st = 0.25 z_st,
                          log z_st ~ N(-s^2/2, s^2); transmission and shock
                          relative variances each about a_obs/2
  S7a/S7b scale_dep       S3 with sigma^2(C) = sigma0^2 (EC/45)^(-gamma),
                          gamma = 0.25 / 0.5, EC the expected count
Per generator (REPS panels): loss fit a, b, b/a; true-error fit
((x - R_i)/xbar)^2; unthinned b/b_P; within-state slope b_w; week-centred
intercept; affine-minus-power score contrast and the quantile of the observed
contrast; retained-row count quantiles; thinning (THIN_OUTER panels x
THIN_INNER draws per q). Branching baselines: analytic b_Pbr(rho)/b_P on the
observed rows and pure Poisson branching-chain simulations. Grid: b/b_P over
(rho, k) with S1 heterogeneity, count-aligned.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from ma_models import Log, summ, write
from revision_data import load_sample, wl_R
from revision_design import contrast
from revision_empirical import centered_loss
from revision_stats import fit_affine, fit_loss, fit_loss_batch, relative_loss
from s03b_thinning import QS, count_vector, poisson_reference
from s03c_state_fe import fit_group_affine

SEED = 2026110351
MIN_COUNT = 5  # analyst filter of the R_A sample (C_t >= 5)
TARGET_D = 45.0
GI_MEAN, GI_SD = 3.6, 1.6
EXP = Path(__file__).resolve().parents[1]
C03 = EXP.parent / "revision2/checks/c03_baselines_counts_s6.json"


class Design:
    def __init__(self, e):
        self.e = e
        keys, counts, cur, fol = count_vector(e)
        self.keys, self.obs_counts, self.cur, self.fol = keys, counts.astype(float), cur, fol
        nk = len(keys)
        start = np.ones(nk, bool)
        for j in range(1, nk):
            (s0, d0), (s1, d1) = keys[j - 1], keys[j]
            start[j] = not (s0 == s1 and (d1 - d0).days == 7)
        self.start = start
        self.run = np.cumsum(start) - 1
        self.run_n = np.bincount(self.run)
        # row index whose current key is j and following key is j+1 (scored transition)
        self.row_of_transition = np.full(nk, -1)
        for i, (c, f) in enumerate(zip(cur, fol)):
            assert f == c + 1 and not start[f]
            self.row_of_transition[c] = i
        self.next_row = np.array([self.row_of_transition[f] for f in fol])
        self.n = len(e)
        self.D = e.C_t.to_numpy(float)
        self.state = e.state.to_numpy()
        self.week = e.week_end.to_numpy()
        self.last_shock = None

    def transition_g(self, g_rows, g0):
        """Growth factor for each key-to-next-key transition."""
        gt = np.full(len(self.keys), g0)
        has = self.row_of_transition >= 0
        gt[has] = g_rows[self.row_of_transition[has]]
        return gt

    def mean_path(self, g_rows, g0, lam=1.0):
        """Run paths carried by the true growth factors, each run level-anchored to
        the geometric mean of its observed counts, times the global scale lam."""
        gt = self.transition_g(g_rows, g0)
        mu = np.empty(len(self.keys))
        for j in range(len(self.keys)):
            mu[j] = max(self.obs_counts[j], .5) if self.start[j] else mu[j - 1] * gt[j - 1]
        lo, lm = np.log(np.maximum(self.obs_counts, .5)), np.log(mu)
        adj = np.bincount(self.run, lo) / self.run_n - np.bincount(self.run, lm) / self.run_n
        return lam * mu * np.exp(adj[self.run]), gt

    def draw(self, rng, g_rows, g0, reps, k=None, rho=.25, lam=1.0, shock_s=0.0):
        """Return count matrix (reps, n_keys). k=None: independent Poisson counts;
        otherwise a latent branching chain (k=inf: Poisson offspring) reported with
        probability rho, or min(1, rho z) with a lognormal state-week shock z."""
        mu, gt = self.mean_path(g_rows, g0, lam)
        self.last_shock = None
        if k is None:
            return rng.poisson(np.broadcast_to(mu, (reps, len(mu)))).astype(float)
        latent = np.empty((reps, len(mu)))
        for j in range(len(mu)):
            if self.start[j]:
                latent[:, j] = rng.poisson(mu[j] / rho, reps)
            else:
                parent = latent[:, j - 1]
                if np.isinf(k):
                    latent[:, j] = rng.poisson(parent * gt[j - 1])
                else:
                    shape = np.maximum(k * parent, 1e-300)
                    latent[:, j] = np.where(parent > 0, rng.poisson(rng.gamma(shape, gt[j - 1] / k)), 0)
        if shock_s > 0:
            z = np.exp(shock_s * rng.standard_normal(latent.shape) - shock_s ** 2 / 2)
            self.last_shock = z
            return rng.binomial(latent.astype(np.int64), np.minimum(1.0, rho * z)).astype(float)
        return latent if rho >= 1 else rng.binomial(latent.astype(np.int64), rho).astype(float)

    def rows(self, counts):
        return counts[..., self.cur], counts[..., self.fol]


def rel_var(r):
    v = float(np.var(r) / np.mean(r) ** 2)
    return 0.0 if v < 1e-20 else v


def bisect(f, lo, hi, iters=80):
    """Root of an increasing function on [lo, hi]."""
    for _ in range(iters):
        mid = (lo + hi) / 2
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def calibrate_sigma(target, z, g0):
    return bisect(lambda s: rel_var(wl_R(g0 * np.exp(s * z - s ** 2 / 2))) - target, 1e-4, 2.0)


def kappa(g):
    return (GI_MEAN / 7) / (1 + (GI_SD ** 2 / GI_MEAN) * np.log(g) / 7)


def bp_of(ct, cn):
    return poisson_reference(pd.DataFrame({"C_t": ct, "C_next": cn}))


def branching_reference(ct, cn, rho):
    """Analytic plug-in reference for a Poisson branching chain with binomial reporting:
    mean kappa^2 [(1 - rho) + (C_t + .5)/(C_next + .5)] (independent Poisson: 1 + 1/g)."""
    g = (cn + .5) / (ct + .5)
    return float(np.mean(kappa(g) ** 2 * ((1 - rho) + 1 / g)))


def thinning_c(rng, d, counts, inner):
    """Shared-count thinning of one panel, mirroring s03b (full subset).

    As in s03b the thinned fits keep the unthinned analyst rows and regress on
    the unthinned current count.
    """
    ct, cn = d.rows(counts)
    keep = ct >= MIN_COUNT
    D = ct[keep]
    base = fit_loss(wl_R((cn[keep] + .5) / (ct[keep] + .5)), D)
    inv = 1 / QS - 1
    db = []
    for q in QS:
        thin = rng.binomial(np.broadcast_to(counts.astype(np.int64), (inner, len(counts))), q)
        tct, tcn = d.rows(thin)
        ft = fit_loss_batch(wl_R((tcn[:, keep] + .5) / (tct[:, keep] + .5)), D)
        db.append(float(np.nanmean(ft["b"]) - base["b"]))
    db = np.array(db)
    c = float(inv @ db / (inv @ inv))
    bp = bp_of(ct[keep], cn[keep])
    return c, c / bp, db / inv, base["b"] / bp


class Generator:
    """Transmission law (g_rows as a function of a variance multiplier amp and the
    scale lam), offspring/reporting law, and optional reporting shock."""

    def __init__(self, d, g0, z, kind, k=None, rho=.25, half=False, gamma=None, kap=None, g_fixed=None):
        self.d, self.g0, self.z, self.kind, self.k, self.rho = d, g0, z, kind, k, rho
        self.g_fixed = g_fixed
        self.half, self.gamma, self.kap = half, gamma, kap
        self.amp, self.lam, self.sigma, self.shock_s = 1.0, 1.0, 0.0, 0.0

    @property
    def het(self):
        return self.kind != "none"

    def set(self, amp=None, lam=None):
        self.amp = self.amp if amp is None else amp
        self.lam = self.lam if lam is None else lam
        d, g0 = self.d, self.g0
        if not self.het:
            self.g, self.sigma = (np.full(d.n, g0) if self.g_fixed is None else self.g_fixed), 0.0
            return self
        target = self.amp * GLOBAL["a_obs"] * (.5 if self.half else 1.0)
        if self.gamma is None:
            self.sigma = calibrate_sigma(target, self.z, g0)
            self.g = g0 * np.exp(self.sigma * self.z - self.sigma ** 2 / 2)
        else:
            ec = d.mean_path(np.full(d.n, g0), g0, self.lam)[0][d.cur]
            scale = (ec / TARGET_D) ** (-self.gamma / 2)

            def g_of(s0):
                s = s0 * scale
                return g0 * np.exp(s * self.z - s ** 2 / 2)
            self.sigma = bisect(lambda s0: rel_var(wl_R(g_of(s0))) - target, 1e-4, 3.0)
            self.g = g_of(self.sigma)
        if self.half:  # S6: shock variance 2 kappa^2 s^2 equal to the transmission share
            self.shock_s = float(np.sqrt(target / (2 * self.kap ** 2)))
        return self

    def draw(self, rng, reps):
        return self.d.draw(rng, self.g, self.g0, reps, k=self.k, rho=self.rho, lam=self.lam, shock_s=self.shock_s)


GLOBAL = {}


def retained_median(gen, rng, reps):
    ct, _ = gen.d.rows(gen.draw(rng, reps))
    return float(np.mean([np.median(r[r >= MIN_COUNT]) for r in ct]))


def fit_lambda(gen, stream, panels=50):
    """Bisection on lam in [0.5, 1.5] (common random numbers, at most 20 steps)."""
    lo, hi, best = .5, 1.5, None
    for _ in range(20):
        mid = (lo + hi) / 2
        m = retained_median(gen.set(lam=mid), np.random.default_rng([SEED, 50, stream]), panels)
        if best is None or abs(m - TARGET_D) < abs(best[1] - TARGET_D):
            best = (mid, m)
        if abs(m - TARGET_D) <= .25:
            break
        lo, hi = (mid, hi) if m < TARGET_D else (lo, mid)
    gen.set(lam=best[0])
    return {"lambda": best[0], "median_D_50panels": best[1], "at_bound": bool(best[0] < .51 or best[0] > 1.49)}


def mean_fitted_a(gen, stream, panels=40):
    rng = np.random.default_rng([SEED, 51, stream])
    ct, cn = gen.d.rows(gen.draw(rng, panels))
    x = wl_R((cn + .5) / (ct + .5))
    a = []
    for r in range(panels):
        keep = ct[r] >= MIN_COUNT
        f = fit_loss(x[r][keep], ct[r][keep])
        if f["converged"]:
            a.append(f["a"])
    return float(np.mean(a))


def fit_amp(gen, stream, log, tag):
    """Secant on the variance multiplier so the MC mean of fitted a equals a_obs."""
    a_obs = GLOBAL["a_obs"]
    pts = [(1.0, mean_fitted_a(gen.set(amp=1.0), stream))]
    amp = max(.2, min(2.0, a_obs / pts[0][1]))
    for _ in range(6):
        pts.append((amp, mean_fitted_a(gen.set(amp=amp), stream)))
        log(tag, "amp", amp, "mean a", pts[-1][1])
        if abs(pts[-1][1] / a_obs - 1) < .005:
            break
        (x0, y0), (x1, y1) = pts[-2], pts[-1]
        amp = x1 + (a_obs - y1) * (x1 - x0) / (y1 - y0) if y1 != y0 else x1
        amp = max(.1, min(3.0, amp))
    best = min(pts, key=lambda p: abs(p[1] - a_obs))
    gen.set(amp=best[0])
    return {"amp": best[0], "mean_a_40panels": best[1], "iterations": len(pts)}


def scenario_panels(gen, rng, reps, obs_contrast, heavy=True):
    """Analyst view of each synthetic panel: the analyst sees only simulated
    counts, so the scale is the simulated C_t and the R_A filter C_t >= 5 is
    applied to the simulated counts (rows of the real design otherwise)."""
    d = gen.d
    counts = gen.draw(rng, reps)
    shock = d.last_shock
    ct, cn = d.rows(counts)
    x = wl_R((cn + .5) / (ct + .5))
    target = wl_R(gen.g)
    names = ["a", "b", "bP", "n", "D_p10", "D_median", "D_p90", "truth_relvar_kept", "vH_kept", "lag1_corr"]
    if heavy:
        names += ["err_a", "err_b", "bw", "a_week", "a_pooled_c", "contrast"]
    out = {nm: np.full(reps, np.nan) for nm in names}
    for r in range(reps):
        keep = ct[r] >= MIN_COUNT
        D, xr = ct[r][keep], x[r][keep]
        out["n"][r] = keep.sum()
        out["D_p10"][r], out["D_median"][r], out["D_p90"][r] = np.percentile(D, (10, 50, 90))
        out["truth_relvar_kept"][r] = rel_var(target[keep])
        if shock is not None:
            xinf = wl_R(gen.g * shock[r][d.fol] / shock[r][d.cur])
            eh = xinf - target
            out["vH_kept"][r] = float(np.mean((eh[keep] / xinf[keep].mean()) ** 2))
            m = keep & (d.next_row >= 0)
            m[m] &= keep[d.next_row[m]]
            out["lag1_corr"][r] = float(np.corrcoef(eh[m], eh[d.next_row[m]])[0, 1])
        else:
            out["vH_kept"][r] = 0.0
        f = fit_loss(xr, D)
        if f["converged"]:
            out["a"][r], out["b"][r] = f["a"], f["b"]
        out["bP"][r] = bp_of(D, cn[r][keep])
        if not heavy:
            continue
        xbar = xr.mean()
        fe = fit_affine(((xr - target[keep]) / xbar) ** 2, D)
        if fe["converged"]:
            out["err_a"][r], out["err_b"][r] = fe["a"], fe["b"]
        loss, _ = relative_loss(xr)
        fw = fit_group_affine(loss, D, d.state[keep])
        if fw["converged"]:
            out["bw"][r] = fw["b"]
        fr = d.e.loc[keep].assign(x=xr, C_t=D, D=D)
        for mode, key in (("week", "a_week"), ("pooled", "a_pooled_c")):
            L, valid, _ = centered_loss(fr, mode)
            fc = fit_affine(L[valid], D[valid])
            if fc["converged"]:
                out[key][r] = fc["a"]
        out["contrast"][r] = contrast(loss, D)
    out["b_over_bP"] = out["b"] / out["bP"]
    if heavy:
        out["bw_over_bP"] = out["bw"] / out["bP"]
    return out, counts


def ratio_of_means(num, den):
    num, den = np.asarray(num, float), np.asarray(den, float)
    ok = np.isfinite(num) & np.isfinite(den)
    m = float(np.mean(den[ok])) if ok.any() else np.nan
    return float(np.mean(num[ok]) / m) if ok.any() and m > 0 else None


def main(reps, thin_outer, thin_inner, grid_reps, sim_reps, outdir, logdir, log_name):
    log = Log(Path(logdir) / log_name)
    e = load_sample("RA")
    d = Design(e)
    log("design rows", d.n, "unique state-weeks", len(d.keys), "runs", int(d.start.sum()))
    obs = fit_loss(e.x.to_numpy(), d.D)
    a_obs, b_obs = obs["a"], obs["b"]
    GLOBAL["a_obs"] = a_obs
    loss_obs, _ = relative_loss(e.x.to_numpy())
    obs_contrast = contrast(loss_obs, d.D)
    ct_obs, cn_obs = e.C_t.to_numpy(float), e.C_next.to_numpy(float)
    g_obs = (cn_obs + .5) / (ct_obs + .5)
    g0 = float(np.exp(np.mean(np.log(g_obs))))
    bP_obs = poisson_reference(e)
    kap = float(np.mean(kappa(g_obs)))
    fw_obs = fit_group_affine(loss_obs, d.D, d.state)
    Lw, vw, _ = centered_loss(e, "week")
    a_week_obs = fit_affine(Lw[vw], d.D[vw])["a"]
    dq = np.percentile(d.D, (10, 50, 90))
    log("observed a", a_obs, "b", b_obs, "b/a", b_obs / a_obs, "bP", bP_obs, "bw", fw_obs["b"], "contrast", obs_contrast,
        "D quantiles", dq.tolist(), "kappa mean", kap)

    z = np.random.default_rng([SEED, 1]).standard_normal(d.n)
    rs = np.random.default_rng([SEED, 2])
    states, si = np.unique(d.state, return_inverse=True)
    weeks, wi = np.unique(d.week, return_inverse=True)
    us, vw_ = rs.standard_normal(len(states)), rs.standard_normal(len(weeks))
    zs = (us[si] + vw_[wi]) / np.sqrt(2)
    zs = (zs - zs.mean()) / zs.std()
    K = 0.3  # offspring dispersion selected on the k grid in the first audit round
    gens = {"S0_poisson": Generator(d, g0, z, "none"),
            "S1_heterogeneity": Generator(d, g0, z, "iid"),
            "S2_nb_reporting": Generator(d, g0, z, "none", k=K),
            "S3_heterogeneity_nb": Generator(d, g0, z, "iid", k=K),
            "S4_structured": Generator(d, g0, zs, "structured"),
            "S5_structured_nb": Generator(d, g0, zs, "structured", k=K),
            "S2P_branching_poisson": Generator(d, g0, z, "none", k=np.inf),
            "S6_reporting_shock": Generator(d, g0, z, "iid", k=K, half=True, kap=kap),
            "S7a_scale_dep": Generator(d, g0, z, "iid", k=K, gamma=.25),
            "S7b_scale_dep": Generator(d, g0, z, "iid", k=K, gamma=.5)}
    alignment = {}
    for si_, (name, gen) in enumerate(gens.items()):
        gen.set(amp=1.0, lam=1.0)
        l1 = fit_lambda(gen, si_)
        am = fit_amp(gen, si_, log, name) if gen.het else None
        l2 = fit_lambda(gen, si_) if gen.het else l1
        alignment[name] = {"lambda_first": l1, "amp": am, "lambda": l2, "sigma": gen.sigma, "shock_s": gen.shock_s}
        log(name, "lambda", gen.lam, "median D (50 panels)", l2["median_D_50panels"], "amp", gen.amp, "sigma", gen.sigma, "s", gen.shock_s)

    results, draws = {}, {}
    for si_, (name, gen) in enumerate(gens.items()):
        rng = np.random.default_rng([SEED, 10, si_])
        o, _ = scenario_panels(gen, rng, reps, obs_contrast)
        truth_v = rel_var(wl_R(gen.g))
        true_total = o["truth_relvar_kept"] + o["vH_kept"]
        ratio = np.where(o["a"] > 0, o["b"] / o["a"], np.nan)
        res = {"k": None if gen.k is None else (None if np.isinf(gen.k) else gen.k),
               "offspring": None if gen.k is None else ("poisson" if np.isinf(gen.k) else "negative_binomial"),
               "rho": gen.rho if gen.k is not None else None, "gamma": gen.gamma,
               "lambda": gen.lam, "sigma": gen.sigma if gen.het else 0.0, "amp": gen.amp if gen.het else None,
               "alignment": alignment[name], "true_relative_variance": truth_v,
               "a": summ(o["a"]), "b": summ(o["b"]), "b_over_a": summ(ratio),
               "share_a_le_0": float(np.mean(o["a"] <= 0)),
               "a_recovery": summ(o["a"] / truth_v) if truth_v > 0 else None,
               "truth_relvar_kept_rows": summ(o["truth_relvar_kept"]),
               "a_recovery_kept_rows": summ(o["a"] / o["truth_relvar_kept"]) if truth_v > 0 else None,
               "true_total_kept_rows": summ(true_total),
               "a_over_true_total": {"point": ratio_of_means(o["a"], true_total) if np.nanmean(true_total) > 0 else None,
                                     **(summ(o["a"] / true_total) if np.nanmean(true_total) > 0 else {})},
               "n_rows": summ(o["n"]), "n_retained": summ(o["n"]),
               "D_p10": summ(o["D_p10"]), "D_median": summ(o["D_median"]), "D_p90": summ(o["D_p90"]),
               "err_a": summ(o["err_a"]), "err_b": summ(o["err_b"]),
               "truth_ratio_errb_over_v": summ(o["err_b"] / truth_v) if truth_v > 0 else None,
               "bP": summ(o["bP"]), "b_over_bP": summ(o["b_over_bP"]),
               "bw": summ(o["bw"]), "bw_over_bP": summ(o["bw_over_bP"]),
               "pooled_over_within_slope": {"point": ratio_of_means(o["b_over_bP"], o["bw_over_bP"]), **summ(o["b"] / o["bw"])},
               "a_week": summ(o["a_week"]),
               "week_reduction": summ(1 - o["a_week"] / o["a_pooled_c"]),
               "contrast": summ(o["contrast"]),
               "observed_contrast_quantile": float(np.mean(o["contrast"][np.isfinite(o["contrast"])] <= obs_contrast))}
        res["shape_contrast_quantile"] = res["observed_contrast_quantile"]
        if gen.shock_s > 0:
            res["reporting_shock"] = {
                "s": gen.shock_s, "tau2": float(np.exp(gen.shock_s ** 2) - 1), "kappa_mean": kap,
                "vH_analytic_2kappa2s2": float(2 * kap ** 2 * gen.shock_s ** 2),
                "vR_rel": summ(o["truth_relvar_kept"]), "vH_rel": summ(o["vH_kept"]),
                "a_over_vR": {"point": ratio_of_means(o["a"], o["truth_relvar_kept"]), **summ(o["a"] / o["truth_relvar_kept"])},
                "err_a_over_vH": {"point": ratio_of_means(o["err_a"], o["vH_kept"]), **summ(o["err_a"] / o["vH_kept"])},
                "lag1_corr": summ(o["lag1_corr"])}
        log(name, "a", res["a"]["mean"], "truth", truth_v, "b", res["b"]["mean"], "b/a", res["b_over_a"].get("mean"),
            "err_a", res["err_a"]["mean"], "b/bP", res["b_over_bP"]["mean"], "bw/bP", res["bw_over_bP"]["mean"],
            "D p10/50/90", res["D_p10"]["mean"], res["D_median"]["mean"], res["D_p90"]["mean"], "n", res["n_rows"]["mean"],
            "a/true", res["a_over_true_total"]["point"], "contrast_q", res["observed_contrast_quantile"])
        if "reporting_shock" in res:
            rsx = res["reporting_shock"]
            log(name, "vR", rsx["vR_rel"]["mean"], "vH", rsx["vH_rel"]["mean"], "a/vR", rsx["a_over_vR"]["point"],
                "err_a/vH", rsx["err_a_over_vH"]["point"], "lag1", rsx["lag1_corr"]["mean"])
        cs, crs, rates, b1 = [], [], [], []
        for t in range(thin_outer):
            counts = gen.draw(rng, 1)[0]
            c, cr, rate, bb = thinning_c(rng, d, counts, thin_inner)
            cs.append(c); crs.append(cr); rates.append(rate); b1.append(bb)
            if (t + 1) % 25 == 0:
                log(name, "thinning", t + 1, "c/bP mean so far", float(np.mean(crs)))
        rates = np.array(rates)
        res["thinning"] = {"c": summ(cs), "c_over_bP": summ(crs), "unthinned_b_over_bP": summ(b1),
                           "per_q_rate": {f"{q:.1f}": summ(rates[:, i]) for i, q in enumerate(QS)}}
        results[name] = res
        for key in ("a", "b", "err_a", "bw", "contrast", "b_over_bP", "D_median", "n"):
            draws[f"{name}__{key}"] = o[key]
        draws[f"{name}__c_over_bP"] = np.array(crs)
        write(Path(outdir) / "ma_audit_partial.json", {"scenarios": results})
        log(name, "done c/bP", res["thinning"]["c_over_bP"]["mean"])

    # branching-chain Poisson baselines (spec A2): analytic on observed rows; pure simulation, aligned
    baseline = {"independent_poisson_bP": bP_obs, "analytic": {}, "simulation": {}}
    for ri, rho in enumerate((0.1, 0.25, 0.5, 1.0)):
        bpr = branching_reference(ct_obs, cn_obs, rho)
        baseline["analytic"][f"rho{rho:g}"] = {"bPbr": bpr, "over_bP": bpr / bP_obs}
        gen = Generator(d, g0, z, "none", k=np.inf, rho=rho).set(amp=1.0, lam=1.0)
        al = fit_lambda(gen, 30 + ri)
        o, _ = scenario_panels(gen, np.random.default_rng([SEED, 30, ri]), sim_reps, obs_contrast, heavy=False)
        baseline["simulation"][f"rho{rho:g}"] = {"lambda": gen.lam, "alignment": al, "b_over_bP": summ(o["b_over_bP"]),
                                                 "b": summ(o["b"]), "a": summ(o["a"]), "D_median": summ(o["D_median"]),
                                                 "n_rows": summ(o["n"])}
        log("baseline rho", rho, "analytic over_bP", bpr / bP_obs, "sim b/bP", baseline["simulation"][f"rho{rho:g}"]["b_over_bP"]["mean"],
            "lambda", gen.lam)

    g_het = gens["S1_heterogeneity"].g
    grid = []
    for ri, rho in enumerate((0.1, 0.25, 0.5, 1.0)):
        for ki, k in enumerate((0.1, 0.3, 1.0, 3.0, np.inf)):
            gen = Generator(d, g0, z, "none", k=k, rho=rho, g_fixed=g_het).set(lam=1.0)  # S1 heterogeneity, fixed
            al = fit_lambda(gen, 100 + 10 * ri + ki)
            o, _ = scenario_panels(gen, np.random.default_rng([SEED, 20, ri, ki]), grid_reps, obs_contrast, heavy=False)
            grid.append({"rho": rho, "k": None if np.isinf(k) else k, "offspring": "poisson" if np.isinf(k) else "negative_binomial",
                         "lambda": gen.lam, "median_D_50panels": al["median_D_50panels"],
                         "branching_r": baseline["analytic"][f"rho{rho:g}"]["over_bP"],
                         "b": summ(o["b"]), "b_over_bP": summ(o["b_over_bP"]), "a": summ(o["a"]), "D_median": summ(o["D_median"])})
            log("grid rho", rho, "k", k, "lambda", gen.lam, "b/bP", grid[-1]["b_over_bP"]["mean"])

    # self-checks (spec A1(5))
    checks = []
    c03 = json.load(open(C03))["observed"]
    obs_q = {"n": d.n, "p10": float(dq[0]), "median": float(dq[1]), "p90": float(dq[2])}
    checks.append({"check": "observed quantiles equal c03", "value": obs_q,
                   "pass": obs_q["n"] == c03["n"] and all(obs_q[k] == c03["D"][k] for k in ("p10", "median", "p90"))})
    for name, res in results.items():
        if name != "S6_reporting_shock":
            v = res["err_a"]["mean"]
            checks.append({"check": f"{name} |err_a| < 1e-4", "value": v, "pass": abs(v) < 1e-4})
    s6 = results["S6_reporting_shock"]
    v = s6["reporting_shock"]["err_a_over_vH"]["point"]
    checks.append({"check": "S6 err_a_over_vH in [0.85, 1.20]", "value": v, "pass": .85 <= v <= 1.2})
    v = s6["a_over_true_total"]["point"]
    checks.append({"check": "S6 a_over_true_total in [0.9, 1.2]", "value": v, "pass": .9 <= v <= 1.2})
    for c in checks:
        log("self-check", c["check"], c["value"], "PASS" if c["pass"] else "FAIL")

    result = {
        "meta": {"seed": SEED, "reps": reps, "thin_outer": thin_outer, "thin_inner": thin_inner, "grid_reps": grid_reps,
                 "baseline_sim_reps": sim_reps, "qs": QS.tolist(), "rows": d.n, "unique_state_weeks": len(d.keys),
                 "runs": int(d.start.sum()), "offspring_k": K,
                 "design": ("real R_A rows; shared counts; run paths carried by the true growth factors and level-anchored to the "
                            "geometric mean of the observed counts of each run; global scale lambda per generator sets the mean "
                            "median retained D to 45 (50 panels, rng [SEED,50,si]); heterogeneity recalibrated so the MC mean of "
                            "fitted a equals observed a (40 panels, rng [SEED,51,si]); analyst scale D = simulated C_t, filter C_t >= 5"),
                 "superseded": "first-round audit (run-start anchoring, no lambda) kept in analysis/superseded/; not used"},
        "observed": {"a": a_obs, "b": b_obs, "b_over_a": b_obs / a_obs, "bP": bP_obs, "b_over_bP": b_obs / bP_obs,
                     "n": d.n, "D_p10": float(dq[0]), "D_median": float(dq[1]), "D_p90": float(dq[2]),
                     "bw": fw_obs["b"], "bw_over_bP": fw_obs["b"] / bP_obs, "a_week": a_week_obs,
                     "contrast_affine_minus_power": obs_contrast, "g0": g0, "kappa_mean": kap},
        "alignment": alignment, "scenarios": results, "baseline": baseline, "rho_k_grid": grid,
        "self_checks": checks, "self_checks_pass": all(c["pass"] for c in checks)}
    write(Path(outdir) / "ma_audit.json", result)
    np.savez_compressed(Path(outdir) / "ma_audit_draws.npz", **draws)
    log("audit complete; self-checks", "PASS" if result["self_checks_pass"] else "FAIL")
    if not result["self_checks_pass"]:
        sys.exit(3)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--reps", type=int, default=200)
    p.add_argument("--thin-outer", type=int, default=100)
    p.add_argument("--thin-inner", type=int, default=20)
    p.add_argument("--grid-reps", type=int, default=100)
    p.add_argument("--sim-reps", type=int, default=200)
    p.add_argument("--output-dir", type=Path, default=EXP / "analysis")
    p.add_argument("--log-dir", type=Path, default=EXP / "logs")
    p.add_argument("--log-name", default="rev2_audit.log")
    a = p.parse_args()
    if os.environ.get("SMOKE"):
        a.reps, a.thin_outer, a.thin_inner, a.grid_reps, a.sim_reps = 6, 2, 3, 4, 6
        a.output_dir = a.output_dir / "smoke"
        a.log_name = "rev2_audit_smoke.log"
    main(a.reps, a.thin_outer, a.thin_inner, a.grid_reps, a.sim_reps, a.output_dir, a.log_dir, a.log_name)
