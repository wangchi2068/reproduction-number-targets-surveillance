"""Fresh-draw generators (analysis plan Section 5.2, review M4).

The shipped generators (ma_audit.py, ma_audit_national.py) draw the transmission field z once
(seed [SEED, 1]) and hold it fixed across replicates, so their Monte Carlo SD of the fitted intercept
reflects counting and reporting noise only. The fresh-draw generators below keep each shipped law
(sigma, offspring k, reporting share, scale exponent) but redraw z in every replicate:

  FD-S1      S1 law: iid lognormal heterogeneity, independent Poisson counts
  FD-S8      S8 law: iid field, NB(k = 0.3) branching, rho = 0.25, half reporting shock, scale exponent 0.25
  FD-S9-w    S9 law with national calendar-week share w = 0.25, 0.5, 0.75 (post hoc law: S9 was specified
             after the data were seen; labelled as such in the output)

lambda (median retained count 45) and the variance amplitude (MC mean of a = observed a) are calibrated
with ma_audit.calibrate exactly as for the shipped generators (calibration stream 900 + i). 200 panels per
generator, panel stream [20261105, 40, i]. Per panel: two-way CGM sandwich intervals (t, G - 1 df), b/b_P,
within-state slope, week centring, season Wald statistic. Rows are appended to analysis/rev_generators_panels.jsonl.

Usage: python code/rev_generators.py [--reps 200] [--workers 5]
"""
from __future__ import annotations

import argparse
import time
from multiprocessing import Pool

import numpy as np
from scipy import stats

import ma_audit as A
import ma_calibration as MC
import rev_common as R
from ma_models import Log
from revision_empirical import centered_loss
from revision_stats import fit_affine, relative_loss
from s03c_state_fe import fit_group_affine

OUT_JSONL = R.AN / "rev_generators_panels.jsonl"
OUT_JSON = R.AN / "rev_generators.json"
SPECS = {  # name: (index i, offspring k, transmission share, scale exponent, national share, post hoc)
    "FD-S1": (0, None, 1.0, None, 0.0, False),
    "FD-S8": (1, A.K, .5, .25, 0.0, False),
    "FD-S9-25": (2, A.K, .5, .25, .25, True),
    "FD-S9-50": (3, A.K, .5, .25, .50, True),
    "FD-S9-75": (4, A.K, .5, .25, .75, True),
}
SE_OBS = 0.0020316  # two-way CGM SE of the observed intercept (ma_calibration.json)


class FreshDraw(A.Generator):
    """A shipped law whose heterogeneity field is redrawn in every replicate."""

    def __init__(self, d, g0, z, k, share, gamma, kap, nat):
        super().__init__(d, g0, z, "iid", k=k, share=share, gamma=gamma, kap=kap)
        self.nat = nat
        _, self.wi = np.unique(d.week, return_inverse=True)
        self.nw = self.wi.max() + 1
        self.last_g = None

    def field(self, rng):
        z = rng.standard_normal(self.d.n)
        if self.nat <= 0:
            return z
        vw = rng.standard_normal(self.nw)
        vw = (vw - vw.mean()) / vw.std()
        zn = np.sqrt(1 - self.nat) * z + np.sqrt(self.nat) * vw[self.wi]
        return (zn - zn.mean()) / zn.std()

    def g_of(self, z):
        if self.gamma is None:
            s = self.sigma
        else:
            ec = self.d.mean_path(np.full(self.d.n, self.g0), self.g0, self.lam)[0][self.d.cur]
            s = self.sigma * (ec / A.TARGET_D) ** (-self.gamma / 2)
        return self.g0 * np.exp(s * z - s ** 2 / 2)

    def draw(self, rng, reps):
        out, gs, shocks = [], [], []
        for _ in range(reps):
            g = self.g_of(self.field(rng))
            out.append(self.d.draw(rng, g, self.g0, 1, k=self.k, rho=self.rho, lam=self.lam, shock_s=self.shock_s)[0])
            gs.append(g)
            shocks.append(None if self.d.last_shock is None else self.d.last_shock[0])
        self.last_g, self.last_shocks = np.array(gs), shocks
        return np.array(out)


def make(name, d, obs):
    i, k, share, gamma, nat, _ = SPECS[name]
    z0 = np.random.default_rng([A.SEED, 1]).standard_normal(d.n)  # calibration reference field of the shipped laws
    return FreshDraw(d, obs["g0"], z0, k, share, gamma, obs["kappa_mean"], nat)


def panel_row(gen, counts, g, shock, season_all, tcrit_fn):
    d = gen.d
    ct, cn = d.rows(counts)
    keep = ct >= A.MIN_COUNT
    D = ct[keep]
    x = R.wl((cn[keep] + .5) / (ct[keep] + .5))
    st, wk, se = d.state[keep], d.week[keep], season_all[keep]
    target = R.wl(g)
    row = {"n": int(keep.sum()), "truth_relvar_kept": A.rel_var(target[keep]), "vH_kept": 0.0}
    if shock is not None:
        xinf = R.wl(g * shock[d.fol] / shock[d.cur])
        row["vH_kept"] = float(np.mean(((xinf[keep] - target[keep]) / xinf[keep].mean()) ** 2))
    s = MC.sandwich(x, D, st, wk)
    if s is not None:
        G = min(len(np.unique(st)), len(np.unique(wk)))
        iv = MC.ab_intervals(s, tcrit_fn(G - 1))["two_way"]
        row.update(a=float(s["beta"][0]), b=float(s["beta"][1]), se_a=float(iv["se_a"]), a_ci=iv["a_ci"],
                   b_over_a_ci=iv["b_over_a_ci"])
    row["bP"] = A.bp_of(D, cn[keep])
    loss, _ = relative_loss(x)
    fw = fit_group_affine(loss, D, st)
    row["bw"] = float(fw["b"]) if fw["converged"] else None
    fr = d.e.loc[keep].assign(x=x, C_t=D, D=D)
    for mode, key in (("week", "a_week"), ("pooled", "a_pooled_c")):
        L, valid, _ = centered_loss(fr, mode)
        fc = fit_affine(L[valid], D[valid])
        row[key] = float(fc["a"]) if fc["converged"] else None
    gsw = MC.sandwich(x, D, st, wk, groups=se)
    if gsw is not None:
        row["season_theta"], row["season_V2"] = gsw["beta"].tolist(), gsw["V2"].tolist()
    return row


def run_generator(job):
    name, reps = job
    t0 = time.time()
    log = Log(R.LOGS / f"rev_generators_{name}.log")
    e, d, obs, _, _ = R.setup_design()
    gen = make(name, d, obs)
    i = SPECS[name][0]
    cal = A.calibrate(gen, 900 + i, log, name)
    log(name, "calibrated", "lambda", gen.lam, "amp", gen.amp, "sigma", gen.sigma, "shock_s", gen.shock_s)
    rng = np.random.default_rng([R.SEED, 40, i])
    season_all = d.e.season.astype(str).to_numpy()
    tcrit = lambda df: stats.t.ppf(.975, df)
    path = R.AN / f"rev_generators_{name}.jsonl"
    if path.exists():
        path.unlink()
    for r in range(reps):
        counts = gen.draw(rng, 1)[0]
        row = panel_row(gen, counts, gen.last_g[0], gen.last_shocks[0], season_all, tcrit)
        row.update(generator=name, r=r)
        R.append_jsonl(path, row)
        if r % 50 == 49:
            log(name, "panels", r + 1)
    info = {"generator": name, "lambda": gen.lam, "amp": gen.amp, "sigma": gen.sigma, "shock_s": gen.shock_s,
            "calibration": cal, "seconds": round(time.time() - t0, 1)}
    R.append_jsonl(R.AN / "rev_generators_calibration.jsonl", info)
    return name


def cover(ci, t):
    ci = np.asarray([c if c is not None else [np.nan, np.nan] for c in ci], float)
    t = np.broadcast_to(np.asarray(t, float), (len(ci),))
    ok = np.isfinite(ci).all(1) & np.isfinite(t)
    p = float(np.mean((ci[ok, 0] <= t[ok]) & (t[ok] <= ci[ok, 1])))
    return {"coverage": p, "mcse": float(np.sqrt(p * (1 - p) / ok.sum())), "n": int(ok.sum()),
            "below": float(np.mean(t[ok] < ci[ok, 0])), "above": float(np.mean(t[ok] > ci[ok, 1]))}


def summarise(name, rows, info, obs, W_obs):
    f = lambda k: np.array([np.nan if r.get(k) is None else r[k] for r in rows], float)
    a, b, bP, bw = f("a"), f("b"), f("bP"), f("bw")
    vR, vH = f("truth_relvar_kept"), f("vH_kept")
    tot = vR + vH
    ea, eb = np.nanmean(a), np.nanmean(b)
    sd_a = float(np.nanstd(a, ddof=1))
    a_ci = [r.get("a_ci") for r in rows]
    ba_ci = [r.get("b_over_a_ci") for r in rows]
    r_tr = R.ratio_of_means(a, vR)
    r_tot = R.ratio_of_means(a, tot)
    # season Wald null: centre log a_s at the generator's MC mean per season (as in ma_calibration)
    th = [np.array(r["season_theta"]) for r in rows if "season_theta" in r]
    Vs = [np.array(r["season_V2"]) for r in rows if "season_theta" in r]
    la = np.array([np.log(t[:-1]) if np.all(t[:-1] > 0) else np.full(len(t) - 1, np.nan) for t in th])
    centre = np.nanmean(la, 0)
    Wn = np.array([MC.season_wald(t, V, centre) for t, V in zip(th, Vs)])
    Wn = Wn[np.isfinite(Wn)]
    p_season = float((1 + np.sum(Wn >= W_obs)) / (len(Wn) + 1))
    i, k, share, gamma, nat, post_hoc = SPECS[name]
    return {
        "law": {"offspring_k": k, "transmission_share": share, "scale_exponent": gamma, "national_week_share": nat,
                "field": "redrawn every replicate", "post_hoc_law": post_hoc},
        "calibration": {kk: info[kk] for kk in ("lambda", "amp", "sigma", "shock_s")},
        "panels": len(rows),
        "a": R.mc(a), "b": R.mc(b), "b_over_bP": R.mc(b / bP), "true_transmission_relvar": R.mc(vR), "true_total_relvar": R.mc(tot),
        "mc_sd_a": sd_a, "sd_ratio_observed_se_over_mc_sd": SE_OBS / sd_a, "mc_sd_over_observed_se": sd_a / SE_OBS,
        "mean_sandwich_se_over_mc_sd": float(np.nanmean(f("se_a")) / sd_a),
        "coverage_a_ci_of_estimator_mean": cover(a_ci, ea),
        "coverage_a_ci_of_true_total": cover(a_ci, tot),
        "coverage_b_over_a_ci_of_ratio_of_means": cover(ba_ci, eb / ea),
        "a_over_transmission": {"point": r_tr, **R.mc(a / vR)},
        "a_over_true_total": {"point": r_tot, **R.mc(a / tot)},
        "pooled_over_within_slope": {"point": R.ratio_of_means(b, bw), **R.mc(b / bw)},
        "week_reduction": R.mc(1 - f("a_week") / f("a_pooled_c")),
        "implied_transmission_rel_sd": R.implied_sd(obs["a"], OBS_A_CI, r_tr, *np.nanquantile(a / vR, [.025, .975])),
        "implied_total_rel_sd": R.implied_sd(obs["a"], OBS_A_CI, r_tot, *np.nanquantile(a / tot, [.025, .975])),
        "season_test_generator_p": {"p": p_season, "mcse": float(np.sqrt(p_season * (1 - p_season) / len(Wn))), "n": int(len(Wn)),
                                    "W_obs": W_obs, "null": "log a_s centred at the generator MC mean per season"},
    }


OBS_A_CI = None


def main():
    global OBS_A_CI
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--summarise-only", action="store_true")
    a = ap.parse_args()
    log = Log(R.LOGS / "rev_generators.log")
    names = list(SPECS)
    if not a.summarise_only:
        cal_path = R.AN / "rev_generators_calibration.jsonl"
        if cal_path.exists():
            cal_path.unlink()
        log("rev_generators start", names, "reps", a.reps)
        with Pool(min(a.workers, len(names))) as pool:
            for nm in pool.imap_unordered(run_generator, [(n, a.reps) for n in names]):
                log("finished", nm)
    import json
    cal = json.loads((R.AN / "ma_calibration.json").read_text())
    ob = cal["observed"]
    OBS_A_CI = ob["intervals"]["two_way"]["a_ci"]
    W_obs = ob["season_W_sandwich"]
    _, _, obs, _, _ = R.setup_design()
    infos = {r["generator"]: r for r in R.read_jsonl(R.AN / "rev_generators_calibration.jsonl")}
    gens = {nm: summarise(nm, R.read_jsonl(R.AN / f"rev_generators_{nm}.jsonl"), infos[nm], obs, W_obs) for nm in names}
    ratios = {nm: g["mc_sd_over_observed_se"] for nm, g in gens.items()}
    sc5 = {nm: g["sd_ratio_observed_se_over_mc_sd"] for nm, g in gens.items()}
    sc4 = {nm: g["coverage_a_ci_of_estimator_mean"]["coverage"] for nm, g in gens.items()}
    out = {"meta": {"script": "code/rev_generators.py", "plan_section": "5.2", "reps": a.reps,
                    "seed": {"panels": [R.SEED, 40, "i"], "calibration_stream": "ma_audit [SEED, 50/51, 900 + i]"},
                    "observed_two_way_se_a": SE_OBS, "observed_a_ci_two_way": OBS_A_CI, "W_obs": W_obs,
                    "post_hoc_note": "the S9 law (national calendar-week share) and the shipped 30-cell decision rule were specified after the data were seen"},
           "generators": gens,
           "SC5_variance_mismatch_explained": {"rule": "observed SE / MC SD of a in [0.67, 1.5] for at least one fresh-draw generator",
                                               "values": sc5, "pass": any(.67 <= v <= 1.5 for v in sc5.values())},
           "SC4_sandwich_coverage_ge_0.90": {"rule": "coverage of E(a) by the two-way 95% CI >= 0.90 under every fresh-draw generator",
                                             "values": sc4, "pass": all(v >= .9 for v in sc4.values())},
           "mc_sd_over_observed_se": ratios}
    R.write_json(OUT_JSON, out)
    log("rev_generators done", "SC5", out["SC5_variance_mismatch_explained"]["pass"], sc5, "SC4", out["SC4_sandwich_coverage_ge_0.90"]["pass"], sc4)


if __name__ == "__main__":
    main()
