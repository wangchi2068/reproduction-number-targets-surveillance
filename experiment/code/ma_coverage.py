"""Coverage of the state x week pigeonhole bootstrap on known-truth panels (review M6).

For generators S1 (heterogeneity only), S6 (heterogeneity + reporting shocks)
and S8 (combined), calibrated exactly as in ma_audit, each synthetic panel is
analysed as the application is: gamma FGLS affine fit of the relative loss on
the retained rows, and B pigeonhole bootstrap draws (revision_stats.
bootstrap_weights, L = 1) giving percentile 95% intervals for a, b and b/a.

Reported per generator, with binomial Monte Carlo SEs:
  * coverage of the interval for a of the panel's true total relative variance
    (transmission + current-value error on the retained rows);
  * coverage of the intervals for a, b, b/a of the estimator's own expectation
    (MC mean over panels; for b/a the ratio of the means), which isolates the
    interval's calibration from any bias of the estimator;
  * mean bootstrap SD over the MC SD of the point estimates.
Panels are split into chunks across worker processes; one JSON part per chunk.
"""
from __future__ import annotations

import argparse
import os
from multiprocessing import Pool
from pathlib import Path

import numpy as np

import ma_audit as A
from ma_models import Log, binom_mcse, summ_mc as summ, write
from revision_data import wl_R
from revision_stats import bootstrap_weights, fit_loss

GENS = {"S1_heterogeneity": 1, "S6_reporting_shock": 7, "S8_combo": 10}  # calibration stream index in ma_audit
STATE = {}


def init(log_path):
    STATE["log"] = Log(log_path)
    e, d, obs, z, gens = A.setup(lambda *p: None)
    STATE.update(d=d, obs=obs, gens=gens, panels={})


def panels(name, reps):
    if name not in STATE["panels"]:
        gen = STATE["gens"][name]
        A.calibrate(gen, GENS[name], lambda *p: None, name)
        o, counts = A.scenario_panels(gen, np.random.default_rng([A.SEED, 70, GENS[name]]), reps,
                                      STATE["obs"]["contrast_affine_minus_power"], heavy=False)
        STATE["panels"][name] = (gen, o, counts)
    return STATE["panels"][name]


def run(job):
    name, lo, hi, reps, B = job["name"], job["lo"], job["hi"], job["reps"], job["B"]
    gen, o, counts = panels(name, reps)
    d = gen.d
    ct, cn = d.rows(counts)
    rows = []
    for r in range(lo, hi):
        keep = ct[r] >= A.MIN_COUNT
        D, xr = ct[r][keep], wl_R((cn[r][keep] + .5) / (ct[r][keep] + .5))
        frame = d.e.loc[keep]
        rng = np.random.default_rng([A.SEED, 71, GENS[name], r])
        bd = np.full((B, 2), np.nan)
        for j in range(B):
            f = fit_loss(xr, D, w=bootstrap_weights(frame, rng)[0])
            if f["converged"]:
                bd[j] = f["a"], f["b"]
        rat = np.where(bd[:, 0] > 0, bd[:, 1] / bd[:, 0], np.nan)
        row = {"r": r, "a": o["a"][r], "b": o["b"][r], "true_total": o["truth_relvar_kept"][r] + o["vH_kept"][r],
               "valid_draws": int(np.isfinite(bd[:, 0]).sum()), "nonpositive_a_draws": int(np.sum(bd[:, 0] <= 0))}
        for k, v in (("a", bd[:, 0]), ("b", bd[:, 1]), ("b_over_a", rat)):
            g = v[np.isfinite(v)]
            row[f"{k}_ci"] = np.percentile(g, [2.5, 97.5]).tolist() if len(g) > 10 else [np.nan, np.nan]
            row[f"{k}_bsd"] = float(g.std(ddof=1)) if len(g) > 10 else np.nan
        rows.append(row)
    write(Path(job["outdir"]) / "coverage_parts" / f"{name}_{lo:04d}.json", {"rows": rows})
    STATE["log"](name, "panels", lo, hi, "done")
    return name, rows


def cover(ci, target):
    ci = np.asarray(ci, float)
    t = np.broadcast_to(np.asarray(target, float), (len(ci),))
    ok = np.isfinite(ci).all(1) & np.isfinite(t)
    p = float(np.mean((ci[ok, 0] <= t[ok]) & (t[ok] <= ci[ok, 1])))
    return {"coverage": p, "mcse": binom_mcse(p, int(ok.sum())), "n": int(ok.sum()),
            "below": float(np.mean(t[ok] < ci[ok, 0])), "above": float(np.mean(t[ok] > ci[ok, 1]))}


def main(reps, B, chunk, workers, outdir, logdir, log_name):
    log_path = Path(logdir) / log_name
    log = Log(log_path)
    jobs = [{"name": n, "lo": lo, "hi": min(lo + chunk, reps), "reps": reps, "B": B, "outdir": str(outdir)}
            for n in GENS for lo in range(0, reps, chunk)]
    log("coverage jobs", len(jobs), "panels", reps, "bootstrap draws", B, "workers", workers)
    with Pool(workers, initializer=init, initargs=(log_path,)) as pool:
        out = pool.map(run, jobs, chunksize=1)
    res = {}
    for n in GENS:
        rows = sorted((r for nm, rr in out if nm == n for r in rr), key=lambda r: r["r"])
        a, b = np.array([r["a"] for r in rows]), np.array([r["b"] for r in rows])
        ea, eb = float(np.nanmean(a)), float(np.nanmean(b))
        tt = np.array([r["true_total"] for r in rows])
        ci = {k: np.array([r[f"{k}_ci"] for r in rows]) for k in ("a", "b", "b_over_a")}
        bsd = {k: np.array([r[f"{k}_bsd"] for r in rows]) for k in ("a", "b", "b_over_a")}
        rat = np.where(a > 0, b / a, np.nan)
        res[n] = {"panels": len(rows), "bootstrap_draws": B,
                  "a_vs_true_total": cover(ci["a"], tt),
                  "a_vs_estimator_mean": cover(ci["a"], ea), "b_vs_estimator_mean": cover(ci["b"], eb),
                  "b_over_a_vs_ratio_of_means": cover(ci["b_over_a"], eb / ea),
                  "estimator_mean": {"a": ea, "b": eb, "b_over_a_ratio_of_means": eb / ea},
                  "a_over_true_total": A.ratio_of_means(a, tt),
                  "se_ratio_bootstrap_over_mc": {k: float(np.nanmean(bsd[k]) / np.nanstd(v, ddof=1))
                                                 for k, v in (("a", a), ("b", b), ("b_over_a", rat))},
                  "point": {"a": summ(a), "b": summ(b), "b_over_a": summ(rat)},
                  "valid_draw_share": float(np.mean([r["valid_draws"] / B for r in rows])),
                  "nonpositive_a_draw_share": float(np.mean([r["nonpositive_a_draws"] / B for r in rows]))}
        log(n, "coverage a/true", res[n]["a_vs_true_total"]["coverage"], "a/E", res[n]["a_vs_estimator_mean"]["coverage"],
            "b/E", res[n]["b_vs_estimator_mean"]["coverage"], "b/a", res[n]["b_over_a_vs_ratio_of_means"]["coverage"],
            "se ratio", res[n]["se_ratio_bootstrap_over_mc"])
    write(Path(outdir) / "ma_coverage.json",
          {"meta": {"seed": A.SEED, "panels": reps, "bootstrap_draws": B, "bootstrap": "state x calendar-week pigeonhole, L = 1",
                    "interval": "percentile 95%", "rng": "panels [SEED, 70, stream]; bootstrap [SEED, 71, stream, panel]",
                    "generators": list(GENS)}, "results": res})
    log("coverage complete")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--reps", type=int, default=200)
    p.add_argument("--boot", type=int, default=200)
    p.add_argument("--chunk", type=int, default=25)
    p.add_argument("--workers", type=int, default=6)
    p.add_argument("--output-dir", type=Path, default=A.EXP / "analysis")
    p.add_argument("--log-dir", type=Path, default=A.EXP / "logs")
    p.add_argument("--log-name", default="rev3_coverage.log")
    a = p.parse_args()
    if os.environ.get("SMOKE"):
        a.reps, a.boot, a.chunk, a.output_dir, a.log_name = 6, 20, 3, a.output_dir / "smoke", "rev3_coverage_smoke.log"
    main(a.reps, a.boot, a.chunk, a.workers, a.output_dir, a.log_dir, a.log_name)
