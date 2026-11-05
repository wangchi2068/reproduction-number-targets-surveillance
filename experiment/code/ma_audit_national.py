"""Known-truth audit, third layer: a generator with a national-week component.

S9 is S8 (iid transmission heterogeneity, NB offspring k = 0.3, rho = 0.25,
scale dependence gamma = 0.25, state-week reporting shocks carrying half of the
variance target) in which half of the transmission heterogeneity is a national
calendar-week component shared by all states. The national component is the
week factor of S4 and S5 (same random stream), and the share 0.5 is the S4/S5
design value; neither is tuned to the observed week-centring reduction.

Calibration is exactly as in ma_audit (count scale lambda, then amp so the MC
mean of the fitted intercept equals the observed a; calibration stream 11).
The script reports, on the real R_A design:
  * the main-audit summaries (200 panels): a against transmission and total,
    b/bP, b/a, within-state slope, pooled/within slope ratio, week-centring
    reduction, and thinning c/bP (100 panels x 20 draws);
  * the observed score-contrast quantile over five streams of 200 panels;
  * bootstrap coverage as in ma_coverage (400 panels x 200 pigeonhole draws).
Output: analysis/ma_audit_national.json (parts in analysis/national_parts/).

Run:  python code/ma_audit_national.py [--workers 4]
"""
from __future__ import annotations

import argparse
import os
from multiprocessing import Pool
from pathlib import Path

import numpy as np

import ma_audit as A
import ma_coverage as C
from ma_models import Log, binom_mcse, summ_mc as summ, write
from revision_data import wl_R

NAME, STREAM, NAT_SHARE = "S9_national", 11, 0.5
STREAMS = 5
ST = {}


def national_z(d):
    """Unit-variance draw: sqrt(1 - w) x the iid draw of S1-S8 + sqrt(w) x the S4/S5 national-week factor."""
    z = np.random.default_rng([A.SEED, 1]).standard_normal(d.n)
    rs = np.random.default_rng([A.SEED, 2])
    states, si = np.unique(d.state, return_inverse=True)
    weeks, wi = np.unique(d.week, return_inverse=True)
    rs.standard_normal(len(states))  # state factor of S4 (consumed to keep the stream aligned)
    vw = rs.standard_normal(len(weeks))
    vw = (vw - vw.mean()) / vw.std()
    zn = np.sqrt(1 - NAT_SHARE) * z + np.sqrt(NAT_SHARE) * vw[wi]
    return (zn - zn.mean()) / zn.std()


def make_gen(d, obs):
    return A.Generator(d, obs["g0"], national_z(d), "iid", k=A.K, share=.5, gamma=.25, kap=obs["kappa_mean"])


def init(log_path):
    C.init(log_path)  # loads the design once per worker; C.STATE holds d, obs, gens, log
    d, obs = C.STATE["d"], C.STATE["obs"]
    gen = make_gen(d, obs)
    C.GENS[NAME] = STREAM
    C.STATE["gens"][NAME] = gen
    ST.update(log=C.STATE["log"], d=d, obs=obs, gen=gen, cal=None)


def calibrated():
    if ST["cal"] is None:
        ST["cal"] = A.calibrate(ST["gen"], STREAM, lambda *p: None, NAME)
    return ST["gen"], ST["cal"]


def job_audit(reps, thin_outer, thin_inner, outdir):
    gen, cal = calibrated()
    d, obs, log = ST["d"], ST["obs"], ST["log"]
    rng = np.random.default_rng([A.SEED, 10, STREAM])
    o, _ = A.scenario_panels(gen, rng, reps, obs["contrast_affine_minus_power"])
    true_total = o["truth_relvar_kept"] + o["vH_kept"]
    ratio = np.where(o["a"] > 0, o["b"] / o["a"], np.nan)
    rom = A.ratio_of_means
    a_over_vr = rom(o["a"], o["truth_relvar_kept"])
    res = {"lambda": gen.lam, "amp": gen.amp, "sigma": gen.sigma, "shock_s": gen.shock_s, "national_share": NAT_SHARE,
           "calibration": A.calibration_flags(gen), "alignment": cal,
           "a": summ(o["a"]), "b": summ(o["b"]), "b_over_a": summ(ratio),
           "a_over_transmission": {"point": a_over_vr, **summ(o["a"] / o["truth_relvar_kept"])},
           "a_over_true_total": {"point": rom(o["a"], true_total), **summ(o["a"] / true_total)},
           "err_a_over_vH": {"point": rom(o["err_a"], o["vH_kept"]), **summ(o["err_a"] / o["vH_kept"])},
           "vR_rel": summ(o["truth_relvar_kept"]), "vH_rel": summ(o["vH_kept"]),
           "b_over_bP": summ(o["b_over_bP"]), "bw_over_bP": summ(o["bw_over_bP"]),
           "pooled_over_within_slope": {"point": rom(o["b_over_bP"], o["bw_over_bP"]), **summ(o["b"] / o["bw"])},
           "week_reduction": summ(1 - o["a_week"] / o["a_pooled_c"]),
           "n_retained": summ(o["n"]), "D_p10": summ(o["D_p10"]), "D_median": summ(o["D_median"]), "D_p90": summ(o["D_p90"]),
           "observed_contrast_quantile_200": float(np.mean(o["contrast"][np.isfinite(o["contrast"])] <= obs["contrast_affine_minus_power"]))}
    res["implied_transmission_rel_sd"] = {"point": float(np.sqrt(obs["a"] / a_over_vr)),
                                          "range_from_ratio_limits": [float(np.sqrt(obs["a"] / res["a_over_transmission"]["q975"])),
                                                                      float(np.sqrt(obs["a"] / res["a_over_transmission"]["q025"]))]}
    crs = []
    for t in range(thin_outer):
        counts = gen.draw(rng, 1)[0]
        crs.append(A.thinning_c(rng, d, counts, thin_inner)[1])
    res["thinning"] = {"c_over_bP": summ(crs)}
    write(Path(outdir) / "national_parts" / "audit.json", res)
    log(NAME, "audit a/vR", a_over_vr, "b/bP", res["b_over_bP"]["mean"], "b/a", res["b_over_a"]["mean"],
        "pooled/within", res["pooled_over_within_slope"]["point"], "week red", res["week_reduction"]["mean"],
        "c/bP", res["thinning"]["c_over_bP"]["mean"], "n", res["n_retained"]["mean"])
    return "audit", res


def job_contrast(s, reps, outdir):
    gen, _ = calibrated()
    obs = ST["obs"]
    o, _ = A.scenario_panels(gen, np.random.default_rng([A.SEED, 60, STREAM, s]), reps,
                             obs["contrast_affine_minus_power"], heavy="contrast")
    c = o["contrast"][np.isfinite(o["contrast"])]
    q = float(np.mean(c <= obs["contrast_affine_minus_power"]))
    part = {"stream": s, "observed_contrast_quantile": q, "mcse": binom_mcse(q, len(c)), "n": int(len(c)), "contrast": c.tolist()}
    write(Path(outdir) / "national_parts" / f"contrast_{s}.json", part)
    ST["log"](NAME, "contrast stream", s, "quantile", q)
    return "contrast", part


def run(job):
    if job["kind"] == "audit":
        return job_audit(job["reps"], job["thin_outer"], job["thin_inner"], job["outdir"])
    if job["kind"] == "contrast":
        return job_contrast(job["s"], job["reps"], job["outdir"])
    calibrated()
    return "coverage", C.run(job)[1]


def coverage_summary(rows, B):
    rows = sorted(rows, key=lambda r: r["r"])
    a, b = np.array([r["a"] for r in rows]), np.array([r["b"] for r in rows])
    ea, eb = float(np.nanmean(a)), float(np.nanmean(b))
    tt = np.array([r["true_total"] for r in rows])
    ci = {k: np.array([r[f"{k}_ci"] for r in rows]) for k in ("a", "b", "b_over_a")}
    bsd = {k: np.array([r[f"{k}_bsd"] for r in rows]) for k in ("a", "b", "b_over_a")}
    rat = np.where(a > 0, b / a, np.nan)
    return {"panels": len(rows), "bootstrap_draws": B,
            "a_vs_true_total": C.cover(ci["a"], tt),
            "a_vs_estimator_mean": C.cover(ci["a"], ea), "b_vs_estimator_mean": C.cover(ci["b"], eb),
            "b_over_a_vs_ratio_of_means": C.cover(ci["b_over_a"], eb / ea),
            "estimator_mean": {"a": ea, "b": eb, "b_over_a_ratio_of_means": eb / ea},
            "a_over_true_total": A.ratio_of_means(a, tt),
            "se_ratio_bootstrap_over_mc": {k: float(np.nanmean(bsd[k]) / np.nanstd(v, ddof=1))
                                           for k, v in (("a", a), ("b", b), ("b_over_a", rat))},
            "valid_draw_share": float(np.mean([r["valid_draws"] / B for r in rows]))}


def main(reps, contrast_reps, cov_reps, B, chunk, thin_outer, thin_inner, workers, outdir, logdir, log_name):
    log_path = Path(logdir) / log_name
    log = Log(log_path)
    jobs = [{"kind": "audit", "reps": reps, "thin_outer": thin_outer, "thin_inner": thin_inner, "outdir": str(outdir)}]
    jobs += [{"kind": "contrast", "s": s, "reps": contrast_reps, "outdir": str(outdir)} for s in range(STREAMS)]
    jobs += [{"kind": "coverage", "name": NAME, "lo": lo, "hi": min(lo + chunk, cov_reps), "reps": cov_reps, "B": B,
              "outdir": str(outdir)} for lo in range(0, cov_reps, chunk)]
    log("national audit jobs", len(jobs), "workers", workers, "national share", NAT_SHARE)
    with Pool(workers, initializer=init, initargs=(log_path,)) as pool:
        out = pool.map(run, jobs, chunksize=1)
    audit = next(r for k, r in out if k == "audit")
    parts = sorted((r for k, r in out if k == "contrast"), key=lambda p: p["stream"])
    allc = np.concatenate([np.asarray(p.pop("contrast")) for p in parts])
    init_obs = A.setup(lambda *p: None)[2]
    obs_c = init_obs["contrast_affine_minus_power"]
    qpool = float(np.mean(allc <= obs_c))
    qs = np.array([p["observed_contrast_quantile"] for p in parts])
    cov_rows = [r for k, rr in out if k == "coverage" for r in rr]
    result = {"meta": {"seed": A.SEED, "generator": NAME, "calibration_stream": STREAM, "national_share_of_transmission": NAT_SHARE,
                       "description": "S8 with half of the transmission heterogeneity a national calendar-week component (S4/S5 week factor)",
                       "panels": reps, "contrast_streams": STREAMS, "contrast_panels_per_stream": contrast_reps,
                       "coverage_panels": cov_reps, "bootstrap_draws": B, "thin_outer": thin_outer, "thin_inner": thin_inner,
                       "rng": "audit [SEED,10,11]; contrast [SEED,60,11,s]; coverage panels [SEED,70,11], bootstrap [SEED,71,11,panel]"},
              "observed": {k: init_obs[k] for k in ("a", "b", "b_over_a", "b_over_bP", "bw", "a_week", "contrast_affine_minus_power")},
              "audit": audit,
              "contrast": {"per_stream": parts, "mean": float(qs.mean()), "sd_across_streams": float(qs.std(ddof=1)),
                           "pooled": qpool, "pooled_mcse": binom_mcse(qpool, len(allc)), "summary": summ(allc)},
              "coverage": coverage_summary(cov_rows, B)}
    write(Path(outdir) / "ma_audit_national.json", result)
    log(NAME, "complete; pooled contrast quantile", qpool, "coverage a/E", result["coverage"]["a_vs_estimator_mean"]["coverage"],
        "se ratio", result["coverage"]["se_ratio_bootstrap_over_mc"])


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--reps", type=int, default=200)
    p.add_argument("--contrast-reps", type=int, default=200)
    p.add_argument("--coverage-reps", type=int, default=400)
    p.add_argument("--boot", type=int, default=200)
    p.add_argument("--chunk", type=int, default=50)
    p.add_argument("--thin-outer", type=int, default=100)
    p.add_argument("--thin-inner", type=int, default=20)
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--output-dir", type=Path, default=A.EXP / "analysis")
    p.add_argument("--log-dir", type=Path, default=A.EXP / "logs")
    p.add_argument("--log-name", default="run_audit_national_progress.log")
    a = p.parse_args()
    if os.environ.get("SMOKE"):
        a.reps, a.contrast_reps, a.coverage_reps, a.boot, a.chunk, a.thin_outer, a.thin_inner = 8, 8, 8, 10, 4, 2, 2
        a.output_dir, a.log_name = a.output_dir / "smoke", "run_audit_national_smoke.log"
    main(a.reps, a.contrast_reps, a.coverage_reps, a.boot, a.chunk, a.thin_outer, a.thin_inner, a.workers,
         a.output_dir, a.log_dir, a.log_name)
