"""Reference distribution of the held-out affine-versus-power score contrast (gap-closing round, B-M5).

The observed held-out contrast (ma_extensions.section_validation: leave one
state out, 51 folds; leave one season out, 4 folds; score mean L/mu + log mu
on the held-out rows; positive power-minus-affine favours the affine form) is
compared with its distribution on synthetic panels from every audited
generator (S0-S8 as calibrated in ma_audit, on the same calibration streams as
ma_audit_extra; S9 as in ma_audit_national). Each synthetic panel is analysed
exactly as the data: retained rows C_t >= 5, the same fold design
(fold_scores) and the same weighting (paired_summary).

Per generator and fold design: summary of the simulated held-out contrast
with Monte Carlo SEs, the share of panels in which it favours the affine form,
and the quantile of the observed value (binomial MC SE). Panels:
rng [SEED, 73, stream].
Output: analysis/ma_heldout_sim.json (parts in analysis/heldout_parts/).
Run: python code/ma_heldout_sim.py [--workers 8]
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
from ma_extensions import fold_scores, paired_summary
from ma_models import Log, binom_mcse, summ_mc as summ, write
from revision_data import wl_R

GEN_NAMES = ["S0_poisson", "S1_heterogeneity", "S2_nb_reporting", "S3_heterogeneity_nb", "S4_structured", "S5_structured_nb",
             "S2P_branching_poisson", "S6_reporting_shock", "S7a_scale_dep", "S7b_scale_dep", "S8_combo"]  # ma_audit_extra order = stream
FOLDS = ("state", "season")
ST = {}


def init(log_path):
    e, d, obs, z, gens = A.setup(lambda *p: None)
    gens[N.NAME] = N.make_gen(d, obs)
    ST.update(log=Log(log_path), d=d, obs=obs, gens=gens)


def stream_of(name):
    return N.STREAM if name == N.NAME else GEN_NAMES.index(name)


def run(job):
    name, lo, hi, reps = job["name"], job["lo"], job["hi"], job["reps"]
    d, log = ST["d"], ST["log"]
    gen = ST["gens"][name]
    if ST.get("calibrated") != name:
        A.calibrate(gen, stream_of(name), lambda *p: None, name)
        ST["calibrated"] = name
    counts = gen.draw(np.random.default_rng([A.SEED, 73, stream_of(name)]), reps)
    ct, cn = d.rows(counts)
    rows = []
    for r in range(lo, hi):
        keep = ct[r] >= A.MIN_COUNT
        e = d.e.loc[keep].copy()
        e["x"] = wl_R((cn[r][keep] + .5) / (ct[r][keep] + .5))
        e["D"] = ct[r][keep]
        row = {"r": r}
        for g in FOLDS:
            fr = fold_scores(e, g, np.ones(len(e)))
            if len(fr) > 1:
                p = paired_summary(fr)
                row[g] = p["mean_power_minus_affine"]
                row[f"{g}_folds"] = p["folds"]
                row[f"{g}_affine_better_folds"] = p["affine_better_folds"]
        rows.append(row)
    write(Path(job["outdir"]) / "heldout_parts" / f"{name}_{lo:04d}.json", {"rows": rows})
    log(name, "panels", lo, hi, "done")
    return name, rows


def main(reps, chunk, workers, outdir, logdir, log_name):
    log_path = Path(logdir) / log_name
    log = Log(log_path)
    val = json.loads((Path(outdir) / "ma_extensions.json").read_text())["validation"]
    observed = {g: val[g]["mean_power_minus_affine"] for g in FOLDS}
    names = GEN_NAMES + [N.NAME]
    # one job per generator and chunk; jobs of a generator are contiguous so a worker recalibrates rarely
    jobs = [{"name": n, "lo": lo, "hi": min(lo + chunk, reps), "reps": reps, "outdir": str(outdir)}
            for n in names for lo in range(0, reps, chunk)]
    log("held-out jobs", len(jobs), "panels", reps, "observed", observed, "workers", workers)
    with Pool(workers, initializer=init, initargs=(log_path,)) as pool:
        out = pool.map(run, jobs, chunksize=1)
    res = {}
    for n in names:
        rows = sorted((r for nm, rr in out if nm == n for r in rr), key=lambda r: r["r"])
        g = {"panels": len(rows)}
        for f in FOLDS:
            v = np.array([r.get(f, np.nan) for r in rows], float)
            ok = v[np.isfinite(v)]
            q = float(np.mean(ok <= observed[f]))
            pa = float(np.mean(ok > 0))
            g[f] = {"n": int(len(ok)), "contrast": summ(ok), "share_favouring_affine": pa, "share_mcse": binom_mcse(pa, len(ok)),
                    "observed_quantile": q, "observed_quantile_mcse": binom_mcse(q, len(ok)),
                    "inside_central_95": bool(.025 <= q <= .975)}
        res[n] = g
        log(n, "state q", g["state"]["observed_quantile"], "mean", g["state"]["contrast"]["mean"],
            "season q", g["season"]["observed_quantile"], "mean", g["season"]["contrast"]["mean"])
    write(Path(outdir) / "ma_heldout_sim.json",
          {"meta": {"seed": A.SEED, "panels": reps, "rng": "[SEED, 73, stream]; calibration streams as ma_audit_extra / ma_audit_national",
                    "score": "mean L/mu + log mu on held-out rows; positive power-minus-affine favours affine",
                    "folds": {"state": "leave one state out", "season": "leave one season out"}},
           "observed": observed, "generators": res})
    log("held-out simulation complete")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--reps", type=int, default=200)
    p.add_argument("--chunk", type=int, default=50)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--output-dir", type=Path, default=A.EXP / "analysis")
    p.add_argument("--log-dir", type=Path, default=A.EXP / "logs")
    p.add_argument("--log-name", default="gap_heldout_sim.log")
    a = p.parse_args()
    if os.environ.get("SMOKE"):
        a.reps, a.chunk, a.log_name = 4, 2, "gap_heldout_sim_smoke.log"
        src = A.EXP / "analysis"
        a.output_dir = a.output_dir / "smoke"
        a.output_dir.mkdir(parents=True, exist_ok=True)
        (a.output_dir / "ma_extensions.json").write_text((src / "ma_extensions.json").read_text())
    main(a.reps, a.chunk, a.workers, a.output_dir, a.log_dir, a.log_name)
