"""Decomposition of the pooled vs state fixed-effects b/a gap on one bootstrap stream.

cross-count ratio = (pooled b/a) / (state-FE b/a, mean-intercept convention)
                  = slope ratio (pool_b / state_b) x intercept ratio (a_mean / pool_a)
log_share_slope   = log(slope ratio) / log(cross-count ratio)

L=1 uses B=2000 and rng default_rng([SEED, 1]) on load_sample("RA"): the same
state x week weights as the primary bootstrap in ma_extensions.py, so the pooled
b/a here reproduces primary/RA/b_over_a draw by draw. L=2 and L=3 (circular
moving week blocks within season runs) keep 1000 draws on the same seed.
Per-draw fits come from s03c_state_fe.run_analysis (unchanged); its authors-format
output is written to analysis/fe_rev2/{L1,block2,block3}/ for render_revision.fe_table.

    python ma_fe_decomp.py [--b1 2000] [--b23 1000]
"""
from __future__ import annotations

import argparse
import logging
import os
import time
from pathlib import Path

import numpy as np

from ma_models import boot_ci, write
from revision_data import load_sample
from revision_stats import save_draws, write_json
from s03c_state_fe import _one_draw, run_analysis

SEED = 2026110352
EXP = Path(__file__).resolve().parents[1]
NAMES = ("pooled_a", "pooled_b", "pooled_b_over_a", "state_b", "a_mean", "state_b_over_a",
         "slope_ratio", "intercept_ratio", "crosscount_ratio", "log_share_slope")


def metrics(frame, w, sm):
    v, _ = _one_draw(frame, w, sm)
    pa, pb, sb, am = (v.get(k, np.nan) for k in ("pool_a", "pool_b", "state_b", "a_mean"))
    pa, pb, sb, am = (np.nan if x is None else float(x) for x in (pa, pb, sb, am))
    pr = pb / pa if pa > 0 else np.nan
    sr = pb / sb if sb > 0 else np.nan
    ir = am / pa if pa > 0 and am > 0 else np.nan
    cr = sr * ir
    ls = np.log(sr) / np.log(cr) if np.isfinite(cr) and cr > 0 and abs(np.log(cr)) > 1e-12 and sr > 0 else np.nan
    return [pa, pb, pr, sb, am, sb / am if am > 0 else np.nan, sr, ir, cr, ls]


def run(frame, B, L, log, outdir):
    """Authors' s03c_state_fe.run_analysis on the primary rng stream; the
    decomposition is computed from its stored per-draw values."""
    t0 = time.time()
    result, arrays = run_analysis(frame, B, [SEED, 1], L)
    result["metadata"]["seed"] = [SEED, 1]
    sub = outdir / ("L1" if L == 1 else f"block{L}")
    sub.mkdir(parents=True, exist_ok=True)
    save_draws(sub / "fixed_effects_draws.npz", **arrays)
    write_json(sub / "fixed_effects.json", result)
    names = [str(x) for x in arrays["metric_names"]]
    col = lambda k: arrays["values"][:, names.index(k)]
    states = np.unique(frame.state.to_numpy())
    point = metrics(frame, np.ones(len(frame)), np.ones(len(states), int))
    v = {k: col(k) for k in ("pool_a", "pool_b", "state_b", "a_mean")}
    pa, pb, sb, am = v["pool_a"], v["pool_b"], v["state_b"], v["a_mean"]
    with np.errstate(all="ignore"):
        pr = np.where(pa > 0, pb / pa, np.nan)
        sr = np.where(sb > 0, pb / sb, np.nan)
        ir = np.where((pa > 0) & (am > 0), am / pa, np.nan)
        cr = sr * ir
        ls = np.where(np.isfinite(cr) & (cr > 0) & (np.abs(np.log(cr)) > 1e-12) & (sr > 0), np.log(sr) / np.log(cr), np.nan)
        draws = np.column_stack([pa, pb, pr, sb, am, np.where(am > 0, sb / am, np.nan), sr, ir, cr, ls])
    out = {nm: boot_ci(point[j], draws[:, j], ratio=nm.endswith("ratio") or nm.endswith("b_over_a"))
           for j, nm in enumerate(NAMES)}
    out["meta"] = {"B": B, "block_length": L, "rng": "default_rng([2026110352, 1])", "rows": len(frame),
                   "states": len(states), "seconds": round(time.time() - t0, 1),
                   "authors_format": str((sub / "fixed_effects.json").relative_to(EXP))}
    return out, draws


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--b1", type=int, default=2000)
    p.add_argument("--b23", type=int, default=1000)
    p.add_argument("--output", type=Path, default=EXP / "analysis/ma_fe_decomp.json")
    p.add_argument("--log-name", default=os.environ.get("MA_FE_LOG", "rev2_fe.log"))
    a = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        handlers=[logging.FileHandler(EXP / "logs" / a.log_name), logging.StreamHandler()])
    log = logging.getLogger("fe")
    frame = load_sample("RA")
    res, store = {}, {}
    for L, B in ((1, a.b1), (2, a.b23), (3, a.b23)):
        res[f"L{L}"], store[f"L{L}"] = run(frame, B, L, log, a.output.resolve().parent / "fe_rev2")
        r = res[f"L{L}"]
        log.info("L=%d slope %.3f %s intercept %.3f %s cross %.3f %s log-share %.3f pooled b/a %.3f",
                 L, r["slope_ratio"]["value"], r["slope_ratio"]["ci95"], r["intercept_ratio"]["value"],
                 r["intercept_ratio"]["ci95"], r["crosscount_ratio"]["value"], r["crosscount_ratio"]["ci95"],
                 r["log_share_slope"]["value"], r["pooled_b_over_a"]["value"])
    res["definitions"] = {
        "slope_ratio": "pooled b / state-FE common b",
        "intercept_ratio": "mean state-FE intercept (state-occurrence weighted) / pooled a",
        "crosscount_ratio": "(pooled b/a) / (state-FE b / mean intercept) = slope_ratio x intercept_ratio",
        "log_share_slope": "log(slope_ratio) / log(crosscount_ratio)",
        "interval": "percentile, finite draws (valid_draws reported)"}
    write(a.output, res)
    np.savez_compressed(a.output.with_name("ma_fe_decomp_draws.npz"), names=np.array(NAMES),
                        **{k: v for k, v in store.items()})
    log.info("wrote %s", a.output)


if __name__ == "__main__":
    main()
