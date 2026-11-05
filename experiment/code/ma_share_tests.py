"""Share grid under all six declared tests, with the uncertainty in a propagated (gap-closing round, B-M1).

Grid: scale dependence gamma {0, .25, .5} x transmission share phi
{0, .25, .5, .75, 1} x national-week share of the transmission draw {0, .5}
(30 cells). National share 0 is the ma_audit_extra grid (same calibration
streams 200 + 100 gamma index + share index, same contrast streams); national
share .5 uses the S9 draw of ma_audit_national (calibration streams
220 + 100 gamma index + share index). Each cell is calibrated as in ma_audit
(lambda, then amp so the MC mean of fitted a equals the observed a).

Per cell, the six tests of ma_tests are applied with the same rule (observed
statistic inside the simulated 2.5-97.5% range; contrast: observed quantile
inside [0.025, 0.975]):
  b/b_P, b/a, pooled/within slope b/b_w and the week-centring reduction on 200
  full-measure panels (rng [SEED, 60, stream, 0]); thinning c/b_P on 100 panels
  x 20 thinning draws (rng [SEED, 61, stream]); the score contrast quantile
  pooled over the full-measure stream and contrast streams 1-4 (1000 panels).

Propagation of the uncertainty in a: each cell is recalibrated with the target
set to the lower and upper limits of the calibrated 95% interval for a
(analysis/ma_calibration.json, the interval chosen by its decision rule) on the
same calibration stream, and a / v_R is re-estimated on 200 panels
(rng [SEED, 62, stream, k]). The implied transmission relative SD at each limit
is sqrt(limit / (a / v_R at that limit)); the reported range for a cell runs
from the lower-limit to the upper-limit value, and the point value uses the
observed a.

Output: analysis/ma_share_tests.json (parts in analysis/share_test_parts/).
Run: python code/ma_share_tests.py [--workers 8]
"""
from __future__ import annotations

import argparse
import json
import os
from multiprocessing import Pool
from pathlib import Path

import numpy as np

import ma_audit as A
import ma_audit_national as N
from ma_models import Log, binom_mcse, summ_mc as summ, write

GAMMAS = (0.0, 0.25, 0.5)
SHARES = (0.0, 0.25, 0.5, 0.75, 1.0)
NATIONAL = (0.0, 0.5)
STREAMS = 5
TESTS = ("b_over_bP", "b_over_a", "thinning", "pooled_over_within", "week_reduction", "contrast")
ST = {}


def init(log_path):
    e, d, obs, z, gens = A.setup(lambda *p: None)
    ST.update(log=Log(log_path), d=d, obs=obs, z=z, zn=N.national_z(d), a_obs=A.GLOBAL["a_obs"])


def rng_range(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    lo, hi = np.quantile(v, (.025, .975))
    return float(lo), float(hi)


def run(job):
    d, obs, log = ST["d"], ST["obs"], ST["log"]
    z = ST["zn"] if job["national"] > 0 else ST["z"]
    gen = A.Generator(d, obs["g0"], z, "iid", k=A.K, share=job["share"], gamma=job["gamma"], kap=obs["kappa_mean"])
    stream, oc = job["stream"], obs["contrast_affine_minus_power"]
    A.GLOBAL["a_obs"] = ST["a_obs"]
    cal = A.calibrate(gen, stream, lambda *p: None, job["id"])
    o, _ = A.scenario_panels(gen, np.random.default_rng([A.SEED, 60, stream, 0]), job["reps"], oc, heavy=True)
    contr = [o["contrast"]]
    for s in range(1, STREAMS):
        oc_s, _ = A.scenario_panels(gen, np.random.default_rng([A.SEED, 60, stream, s]), job["reps"], oc, heavy="contrast")
        contr.append(oc_s["contrast"])
    c = np.concatenate(contr)
    c = c[np.isfinite(c)]
    q = float(np.mean(c <= oc))
    rng = np.random.default_rng([A.SEED, 61, stream])
    crs = [A.thinning_c(rng, d, gen.draw(rng, 1)[0], job["thin_inner"])[1] for _ in range(job["thin_outer"])]
    sim = {"b_over_bP": o["b_over_bP"], "b_over_a": np.where(o["a"] > 0, o["b"] / o["a"], np.nan), "thinning": np.asarray(crs, float),
           "pooled_over_within": o["b"] / o["bw"], "week_reduction": 1 - o["a_week"] / o["a_pooled_c"]}
    tests = {}
    for t in TESTS[:-1]:
        lo, hi = rng_range(sim[t])
        ob = job["observed"][t]
        tests[t] = {"observed": ob, "sim_q025": lo, "sim_q975": hi, "pass": bool(lo <= ob <= hi), "sim": summ(sim[t])}
    tests["contrast"] = {"observed_quantile": q, "mcse": binom_mcse(q, len(c)), "n": int(len(c)), "pass": bool(.025 <= q <= .975)}
    rom = A.ratio_of_means
    het = float(np.nanmean(o["truth_relvar_kept"])) > 0
    avr = rom(o["a"], o["truth_relvar_kept"]) if het else None
    res = {"id": job["id"], "gamma": job["gamma"], "share": job["share"], "national": job["national"], "stream": stream,
           "lambda": gen.lam, "amp": gen.amp, "shock_s": gen.shock_s, "alignment": cal, "a": summ(o["a"]),
           "tests": tests, "n_pass": int(sum(v["pass"] for v in tests.values())), "pass_all": all(v["pass"] for v in tests.values()),
           "a_over_vR": avr, "a_over_vR_mc": summ(o["a"] / o["truth_relvar_kept"]) if het else None,
           "implied_sd_point": float(np.sqrt(ST["a_obs"] / avr)) if avr and avr > 0 else None}
    # propagation of the uncertainty in a: recalibrate at the interval limits
    prop = {}
    for k, (lab, tgt) in enumerate((("lower", job["a_limits"][0]), ("upper", job["a_limits"][1]))):
        if not het:
            break
        A.GLOBAL["a_obs"] = tgt
        g2 = A.Generator(d, obs["g0"], z, "iid", k=A.K, share=job["share"], gamma=job["gamma"], kap=obs["kappa_mean"])
        A.calibrate(g2, stream, lambda *p: None, job["id"] + lab)
        o2, _ = A.scenario_panels(g2, np.random.default_rng([A.SEED, 62, stream, k]), job["reps"], oc, heavy=False)
        r2 = rom(o2["a"], o2["truth_relvar_kept"])
        prop[lab] = {"target_a": tgt, "amp": g2.amp, "lambda": g2.lam, "mean_a": float(np.nanmean(o2["a"])), "a_over_vR": r2,
                     "implied_sd": float(np.sqrt(tgt / r2)) if r2 and r2 > 0 else None}
    A.GLOBAL["a_obs"] = ST["a_obs"]
    res["propagation"] = prop
    if prop:
        res["implied_sd_range"] = [prop["lower"]["implied_sd"], prop["upper"]["implied_sd"]]
    write(Path(job["outdir"]) / "share_test_parts" / f"{job['id']}.json", res)
    log(job["id"], "passes", res["n_pass"], {t: v["pass"] for t, v in tests.items()}, "a/vR", avr,
        "implied SD", res["implied_sd_point"], res.get("implied_sd_range"))
    return res


def main(reps, thin_outer, thin_inner, workers, outdir, logdir, log_name):
    log_path = Path(logdir) / log_name
    log = Log(log_path)
    res_all = json.loads((Path(outdir) / "results.json").read_text())
    est = {e["name"]: e for e in res_all["estimates"]}
    observed = {t: est[f"tests/observed/{t}"]["point"] for t in TESTS[:-1]}
    cal = json.loads((Path(outdir) / "ma_calibration.json").read_text())
    if cal["chosen"] == "two_way_sandwich":
        a_lim = cal["applied"]["two_way_sandwich"]["a_ci"]
    else:
        a_lim = cal["applied"]["scaled_bootstrap"]["a_ci_ratio_min"]
    log("observed", observed, "a limits", a_lim, "from", cal["chosen"])
    jobs = []
    for gi, gam in enumerate(GAMMAS):
        for ni, nat in enumerate(NATIONAL):
            for i, phi in enumerate(SHARES):
                jobs.append({"id": f"g{gam:g}_phi{phi:g}_nat{nat:g}", "gamma": gam, "share": phi, "national": nat,
                             "stream": 200 + 100 * gi + 20 * ni + i, "reps": reps, "thin_outer": thin_outer, "thin_inner": thin_inner,
                             "observed": observed, "a_limits": a_lim, "outdir": str(outdir)})
    log("share-test jobs", len(jobs), "panels", reps, "thinning", thin_outer, "x", thin_inner, "workers", workers)
    with Pool(workers, initializer=init, initargs=(log_path,)) as pool:
        out = pool.map(run, jobs, chunksize=1)
    surv = [r for r in out if r["pass_all"]]
    best = max(r["n_pass"] for r in out)
    near = [r["id"] for r in out if r["n_pass"] == best]
    het = [r for r in out if r.get("implied_sd_range")]
    summary = {"cells": len(out), "surviving_all_six": [r["id"] for r in surv], "max_tests_passed": best, "cells_at_max": near,
               "pass_counts_by_test": {t: int(sum(r["tests"][t]["pass"] for r in out)) for t in TESTS},
               "implied_sd_point_range_all_het_cells": [min(r["implied_sd_point"] for r in het), max(r["implied_sd_point"] for r in het)],
               "implied_sd_propagated_range_all_het_cells": [min(r["implied_sd_range"][0] for r in het), max(r["implied_sd_range"][1] for r in het)]}
    if surv:
        sh = [r for r in surv if r.get("implied_sd_range")]
        summary["surviving_shares"] = sorted({r["share"] for r in surv})
        if sh:
            summary["implied_sd_propagated_range_surviving"] = [min(r["implied_sd_range"][0] for r in sh), max(r["implied_sd_range"][1] for r in sh)]
    nm = [r for r in het if r["id"] in near]
    if nm:
        summary["implied_sd_propagated_range_cells_at_max"] = [min(r["implied_sd_range"][0] for r in nm), max(r["implied_sd_range"][1] for r in nm)]
        summary["implied_sd_point_range_cells_at_max"] = [min(r["implied_sd_point"] for r in nm), max(r["implied_sd_point"] for r in nm)]
    write(Path(outdir) / "ma_share_tests.json",
          {"meta": {"seed": A.SEED, "gammas": list(GAMMAS), "shares": list(SHARES), "national": list(NATIONAL), "panels": reps,
                    "contrast_panels": reps * STREAMS, "thinning": [thin_outer, thin_inner], "a_limits": a_lim, "a_limits_source": cal["chosen"],
                    "rule": "observed statistic inside the simulated 2.5-97.5% range; contrast quantile in [0.025, 0.975]",
                    "rng": "calibration 200 + 100 gamma index + 20 national index + share index; panels [SEED, 60, stream, s]; thinning [SEED, 61, stream]; propagation [SEED, 62, stream, k]"},
           "observed": observed, "summary": summary, "cells": out})
    log("share tests complete", summary)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--reps", type=int, default=200)
    p.add_argument("--thin-outer", type=int, default=100)
    p.add_argument("--thin-inner", type=int, default=20)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--output-dir", type=Path, default=A.EXP / "analysis")
    p.add_argument("--log-dir", type=Path, default=A.EXP / "logs")
    p.add_argument("--log-name", default="gap_share_tests.log")
    a = p.parse_args()
    if os.environ.get("SMOKE"):
        a.reps, a.thin_outer, a.thin_inner, a.log_name = 6, 3, 3, "gap_share_tests_smoke.log"
        src = A.EXP / "analysis"
        a.output_dir = a.output_dir / "smoke"
        a.output_dir.mkdir(parents=True, exist_ok=True)
        (a.output_dir / "results.json").write_text((src / "results.json").read_text())
        if not (a.output_dir / "ma_calibration.json").exists():
            raise SystemExit("smoke needs analysis/smoke/ma_calibration.json (run ma_calibration.py smoke first)")
    main(a.reps, a.thin_outer, a.thin_inner, a.workers, a.output_dir, a.log_dir, a.log_name)
