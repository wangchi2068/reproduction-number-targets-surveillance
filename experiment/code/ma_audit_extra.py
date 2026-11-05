"""Known-truth audit, second layer (review M2, M5).

1. (gamma, s) grid around the combined generator S8 (iid heterogeneity, NB
   offspring k = 0.3, scale dependence gamma, reporting shocks). The share phi
   of the total relative variance carried by transmission runs over
   {0, .25, .5, .75, 1}; the remaining 1 - phi is the reporting shock
   (2 kappa^2 s^2 = (1 - phi) amp a_obs), so phi indexes the shock scale s.
   gamma runs over {0, .25, .5} (gamma = 0, phi = .5 is S6; gamma = .25,
   phi = .5 is S8; phi = 1 gives S3, S7a, S7b). "calibrated" tunes amp so the MC
   mean of fitted a equals the observed a; "uncalibrated" (gamma = .25 only)
   fixes amp = 1 (total variance target = observed a). For each
   point it reports what the analyst pipeline returns and whether the observed
   a, b/b_P, b/a and score contrast fall inside the simulated 2.5-97.5% range.
2. Observed contrast quantile of every generator S0-S8 across five independent
   random streams (mean, SD across streams, pooled value with binomial MC SE).

Each point is calibrated exactly as in ma_audit (lambda to the observed median
retained count, then amp) and simulated in contrast mode. One JSON part per
job is written as soon as the job finishes; the merged file is
analysis/ma_audit_extra.json.
"""
from __future__ import annotations

import argparse
import json
import os
from multiprocessing import Pool
from pathlib import Path

import numpy as np

import ma_audit as A
from ma_models import Log, binom_mcse, summ_mc as summ, write

SHARES = (0.0, 0.25, 0.5, 0.75, 1.0)
GAMMAS = (0.0, 0.25, 0.5)
STREAMS = 5
STATE = {}


def init(log_path):
    log = Log(log_path)
    e, d, obs, z, gens = A.setup(lambda *p: None)
    STATE.update(log=log, d=d, obs=obs, z=z, gens=gens)


def make(job):
    d, obs, z = STATE["d"], STATE["obs"], STATE["z"]
    if job["kind"] == "share":
        return A.Generator(d, obs["g0"], z, "iid", k=A.K, share=job["share"], gamma=job["gamma"], kap=obs["kappa_mean"])
    return STATE["gens"][job["name"]]


def compat(values, observed):
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    lo, hi = np.quantile(v, (.025, .975))
    return {"observed": float(observed), "sim_q025": float(lo), "sim_q975": float(hi), "inside": bool(lo <= observed <= hi),
            "observed_quantile": float(np.mean(v <= observed))}


def run(job):
    log, obs = STATE["log"], STATE["obs"]
    gen = make(job)
    cal = A.calibrate(gen, job["stream"], log, job["id"], calibrate_amp=job.get("calibrate_amp", True))
    pooled = {k: [] for k in ("a", "b", "b_over_bP", "contrast", "truth", "vH")}
    per_stream = []
    for s in range(STREAMS):
        o, _ = A.scenario_panels(gen, np.random.default_rng([A.SEED, 60, job["stream"], s]), job["reps"],
                                 obs["contrast_affine_minus_power"], heavy="contrast")
        c = o["contrast"][np.isfinite(o["contrast"])]
        q = float(np.mean(c <= obs["contrast_affine_minus_power"]))
        per_stream.append({"stream": s, "observed_contrast_quantile": q, "mcse": binom_mcse(q, len(c)), "n": int(len(c))})
        for k, src in (("a", "a"), ("b", "b"), ("b_over_bP", "b_over_bP"), ("contrast", "contrast"),
                       ("truth", "truth_relvar_kept"), ("vH", "vH_kept")):
            pooled[k].append(o[src])
        log(job["id"], "stream", s, "contrast quantile", q, "mean a", float(np.nanmean(o["a"])))
    P = {k: np.concatenate(v) for k, v in pooled.items()}
    ratio = np.where(P["a"] > 0, P["b"] / P["a"], np.nan)
    qs = np.array([p["observed_contrast_quantile"] for p in per_stream])
    nc = int(np.isfinite(P["contrast"]).sum())
    qpool = float(np.mean(P["contrast"][np.isfinite(P["contrast"])] <= obs["contrast_affine_minus_power"]))
    total = P["truth"] + P["vH"]
    res = {"id": job["id"], "kind": job["kind"], "share": job.get("share"), "gamma": gen.gamma, "name": job.get("name"),
           "calibrated_amp": job.get("calibrate_amp", True), "panels_per_stream": job["reps"], "streams": STREAMS,
           "lambda": gen.lam, "amp": gen.amp if gen.het else None, "sigma": gen.sigma if gen.het else 0.0, "shock_s": gen.shock_s,
           "calibration": A.calibration_flags(gen, job.get("calibrate_amp", True)), "alignment": cal,
           "a": summ(P["a"]), "b": summ(P["b"]), "b_over_a": summ(ratio), "b_over_bP": summ(P["b_over_bP"]),
           "contrast": summ(P["contrast"]),
           "true_transmission_relvar": summ(P["truth"]), "true_total_relvar": summ(total),
           "a_over_true_total": A.ratio_of_means(P["a"], total),
           "a_over_true_transmission": A.ratio_of_means(P["a"], P["truth"]) if np.nanmean(P["truth"]) > 0 else None,
           "observed_contrast_quantile": {"per_stream": per_stream, "mean": float(qs.mean()), "sd_across_streams": float(qs.std(ddof=1)),
                                          "pooled": qpool, "pooled_mcse": binom_mcse(qpool, nc)},
           "compatibility": {"a": compat(P["a"], obs["a"]), "b_over_bP": compat(P["b_over_bP"], obs["b_over_bP"]),
                             "b_over_a": compat(ratio, obs["b_over_a"]), "contrast": compat(P["contrast"], obs["contrast_affine_minus_power"])}}
    res["compatible_all"] = all(v["inside"] for v in res["compatibility"].values())
    write(Path(job["outdir"]) / "extra_parts" / f"{job['id']}.json", res)
    log(job["id"], "done", "a", res["a"]["mean"], "b/bP", res["b_over_bP"]["mean"], "b/a", res["b_over_a"]["mean"],
        "pooled q", qpool, "compatible", res["compatible_all"])
    return res


def main(reps, workers, outdir, logdir, log_name):
    log_path = Path(logdir) / log_name
    log = Log(log_path)
    gen_names = ["S0_poisson", "S1_heterogeneity", "S2_nb_reporting", "S3_heterogeneity_nb", "S4_structured", "S5_structured_nb",
                 "S2P_branching_poisson", "S6_reporting_shock", "S7a_scale_dep", "S7b_scale_dep", "S8_combo"]
    jobs = []
    for gi, gam in enumerate(GAMMAS):
        for ci, cal in enumerate((True, False)):
            if not cal and gam != .25:
                continue
            for i, phi in enumerate(SHARES):
                jobs.append({"id": f"share{phi:g}_g{gam:g}_{'cal' if cal else 'uncal'}", "kind": "share", "share": phi, "gamma": gam,
                             "calibrate_amp": cal, "stream": 200 + 100 * gi + 10 * ci + i, "reps": reps, "outdir": str(outdir)})
    for si_, n in enumerate(gen_names):
        # same calibration stream as ma_audit, so lambda and amp equal the main audit's values
        jobs.append({"id": f"gen_{n}", "kind": "generator", "name": n, "stream": si_, "reps": reps, "outdir": str(outdir)})
    log("jobs", len(jobs), "panels per stream", reps, "streams", STREAMS, "workers", workers)
    with Pool(workers, initializer=init, initargs=(log_path,)) as pool:
        out = pool.map(run, jobs, chunksize=1)
    by = {r["id"]: r for r in out}
    grid = {f"g{gam:g}": [by[f"share{phi:g}_g{gam:g}_cal"] for phi in SHARES] for gam in GAMMAS}
    uncal = [by[f"share{phi:g}_g0.25_uncal"] for phi in SHARES]
    gens = {n: by[f"gen_{n}"] for n in gen_names}
    _, _, obs, _, _ = A.setup(lambda *p: None)
    result = {"meta": {"seed": A.SEED, "streams": STREAMS, "panels_per_stream": reps, "shares": list(SHARES), "gammas": list(GAMMAS),
                       "share_generator": "iid heterogeneity, NB offspring k = 0.3, reporting probability 0.25, scale dependence gamma, reporting shocks",
                       "identities": "gamma 0 / phi 0.5 = S6; gamma 0.25 / phi 0.5 = S8; phi 1 = S3 (gamma 0), S7a (0.25), S7b (0.5)",
                       "rng_grid": "calibration stream 200 + 100 gamma index + 10 variant + share index",
                       "rng": "calibration streams as ma_audit (generators) or 200+ (share grid); simulation rng [SEED, 60, stream, s]",
                       "compatibility": "observed statistic inside the simulated 2.5-97.5% range of the pooled panels"},
              "observed": obs, "share_grid": {"calibrated": grid, "uncalibrated_gamma0.25": uncal},
              "generator_contrast_streams": {n: {"per_stream": g["observed_contrast_quantile"]["per_stream"],
                                                 "mean": g["observed_contrast_quantile"]["mean"],
                                                 "sd_across_streams": g["observed_contrast_quantile"]["sd_across_streams"],
                                                 "pooled": g["observed_contrast_quantile"]["pooled"],
                                                 "pooled_mcse": g["observed_contrast_quantile"]["pooled_mcse"],
                                                 "lambda": g["lambda"], "amp": g["amp"]} for n, g in gens.items()},
              "generators": gens}
    write(Path(outdir) / "ma_audit_extra.json", result)
    log("extra audit complete")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--reps", type=int, default=200)
    p.add_argument("--workers", type=int, default=6)
    p.add_argument("--output-dir", type=Path, default=A.EXP / "analysis")
    p.add_argument("--log-dir", type=Path, default=A.EXP / "logs")
    p.add_argument("--log-name", default="rev3_audit_extra.log")
    a = p.parse_args()
    if os.environ.get("SMOKE"):
        a.reps, a.output_dir, a.log_name = 4, a.output_dir / "smoke", "rev3_audit_extra_smoke.log"
    main(a.reps, a.workers, a.output_dir, a.log_dir, a.log_name)
