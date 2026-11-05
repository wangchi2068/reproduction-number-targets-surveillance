"""Non-parametric Wild Cluster Bootstrap and Permutation Calibration for Season Test.

Addresses reviewer comment C-M1: provides non-parametric calibrations of the two-way
cluster-robust season Wald test W_CGM that operate directly on empirical data, eliminating
reliance on parametric synthetic generators that have lower dependent variance than the data.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import scipy.stats as stats

from ma_calibration import sandwich, season_wald
from ma_models import binom_mcse
from revision_data import load_sample
from revision_stats import fit_loss, relative_loss
from s03c_state_fe import fit_group_affine

SEED = 20261105


def run_wild_cluster_bootstrap(x, D, st, wk, se, B=2000, seed=SEED):
    """Wild Cluster Bootstrap (Webb 6-point) clustered at State level.

    Imposes null of equal season intercepts by obtaining residuals from the pooled
    affine fit, then perturbing cluster residuals with Webb 6-point weights.
    """
    loss, xbar = relative_loss(x)
    n = len(x)
    f0 = fit_loss(x, D)
    if not f0["converged"]:
        raise RuntimeError("Pooled fit failed to converge under null")
    mu0 = f0["a"] + f0["b"] / D
    res0 = loss - mu0

    states = np.unique(st)
    G = len(states)
    state_idx = [np.where(st == s)[0] for s in states]

    # Webb 6-point weights: {-sqrt(3/2), -1, -sqrt(1/2), sqrt(1/2), 1, sqrt(3/2)}
    webb = np.array([-np.sqrt(1.5), -1.0, -np.sqrt(0.5), np.sqrt(0.5), 1.0, np.sqrt(1.5)])
    rng = np.random.default_rng([seed, 1])

    names, codes = np.unique(se, return_inverse=True)
    X = np.column_stack([np.eye(len(names))[codes], 1 / D])
    k_dim = X.shape[1] + 1

    W_boot = []
    for b in range(B):
        v = rng.choice(webb, size=G)
        loss_b = mu0.copy()
        for g in range(G):
            idx = state_idx[g]
            loss_b[idx] += v[g] * res0[idx]
        loss_b = np.maximum(loss_b, 1e-7)

        fb = fit_group_affine(loss_b, D, se)
        if not fb["converged"] or not np.all(fb["active"]) or not np.all(fb["intercepts"] > 0):
            continue

        beta_b = np.r_[np.asarray(fb["intercepts"], float), fb["b"]]
        mu_b = X @ beta_b
        if np.any(mu_b <= 0):
            continue
        r_b = (loss_b - mu_b) / mu_b**2
        psi_b = np.column_stack([x - xbar, X * r_b[:, None]])

        J = np.zeros((k_dim, k_dim))
        J[0, 0] = -n
        u = (x - xbar) / xbar
        J[1:, 0] = X.T @ ((2 * u * (-x / xbar**2)) / mu_b**2)
        J[1:, 1:] = -(X.T * (1 / mu_b**2)) @ X

        try:
            Ji = np.linalg.inv(J)

            def meat(lab):
                _, ci = np.unique(lab, return_inverse=True)
                S_mat = np.zeros((ci.max() + 1, k_dim))
                np.add.at(S_mat, ci, psi_b)
                return S_mat.T @ S_mat, ci.max() + 1

            MS, GS = meat(st)
            MW, GW = meat(wk)
            MR = psi_b.T @ psi_b
            M2 = GS / (GS - 1) * MS + GW / (GW - 1) * MW - MR
            V2 = Ji @ M2 @ Ji.T
            w_eig, U = np.linalg.eigh((V2 + V2.T) / 2)
            V2_proj = (U * np.maximum(w_eig, 0)) @ U.T
            wb = season_wald(beta_b, V2_proj[1:, 1:])
            if np.isfinite(wb) and wb >= 0:
                W_boot.append(wb)
        except Exception:
            continue

    return np.array(W_boot)


def run_permutation_test(x, D, st, wk, se, B=2000, seed=SEED):
    """Within-state season permutation test.

    Under H0 of no season differences in dispersion, season labels within each state
    cluster are exchangeable.
    """
    states = np.unique(st)
    rng = np.random.default_rng([seed, 2])
    W_perm = []

    for b in range(B):
        se_perm = se.copy()
        for s in states:
            idx = np.where(st == s)[0]
            se_perm[idx] = rng.permutation(se[idx])
        g = sandwich(x, D, st, wk, groups=se_perm)
        if g is not None:
            w = season_wald(g["beta"], g["V2"])
            if np.isfinite(w) and w >= 0:
                W_perm.append(w)

    return np.array(W_perm)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--draws", type=int, default=2000, help="Number of bootstrap/permutation draws")
    parser.add_argument("--seed", type=int, default=SEED, help="Random seed")
    args = parser.parse_args()

    e = load_sample("RA")
    x, D = e.x.to_numpy(float), e.D.to_numpy(float)
    st, wk, se = e.state.to_numpy(), e.week_end.to_numpy(), e.season.astype(str).to_numpy()

    # Observed statistic
    g_obs = sandwich(x, D, st, wk, groups=se)
    W_obs = season_wald(g_obs["beta"], g_obs["V2"])
    asymp_p = float(stats.chi2.sf(W_obs, 3))

    t0 = time.time()
    print(f"Running Wild Cluster Bootstrap (B={args.draws}, seed={args.seed})...")
    w_boot = run_wild_cluster_bootstrap(x, D, st, wk, se, B=args.draws, seed=args.seed)
    t1 = time.time()
    print(f"Wild Cluster Bootstrap finished in {t1 - t0:.1f}s ({len(w_boot)} valid draws).")

    print(f"Running Within-State Permutation Test (B={args.draws}, seed={args.seed})...")
    w_perm = run_permutation_test(x, D, st, wk, se, B=args.draws, seed=args.seed)
    t2 = time.time()
    print(f"Permutation Test finished in {t2 - t1:.1f}s ({len(w_perm)} valid draws).")

    p_wcb = float(np.mean(w_boot >= W_obs))
    p_wcb_mcse = float(binom_mcse(p_wcb, len(w_boot)))
    p_perm = float(np.mean(w_perm >= W_obs))
    p_perm_mcse = float(binom_mcse(p_perm, len(w_perm)))

    out = {
        "observed": {
            "W_CGM": float(W_obs),
            "df": 3,
            "asymptotic_chi2_p": asymp_p,
            "season_intercepts": [float(v) for v in g_obs["beta"][:4]],
            "common_slope": float(g_obs["beta"][4]),
        },
        "wild_cluster_bootstrap": {
            "method": "Webb 6-point Wild Cluster Bootstrap (State clusters)",
            "n_draws": len(w_boot),
            "mean": float(np.mean(w_boot)),
            "median": float(np.median(w_boot)),
            "sd": float(np.std(w_boot, ddof=1)),
            "null_q90": float(np.percentile(w_boot, 90)),
            "null_q95": float(np.percentile(w_boot, 95)),
            "null_q99": float(np.percentile(w_boot, 99)),
            "p_value": p_wcb,
            "p_value_mcse": p_wcb_mcse,
        },
        "permutation_test": {
            "method": "Within-state season permutation",
            "n_draws": len(w_perm),
            "mean": float(np.mean(w_perm)),
            "median": float(np.median(w_perm)),
            "sd": float(np.std(w_perm, ddof=1)),
            "null_q90": float(np.percentile(w_perm, 90)),
            "null_q95": float(np.percentile(w_perm, 95)),
            "null_q99": float(np.percentile(w_perm, 99)),
            "p_value": p_perm,
            "p_value_mcse": p_perm_mcse,
        },
    }

    out_file = Path(__file__).resolve().parents[1] / "analysis" / "ma_wild_bootstrap.json"
    out_file.write_text(json.dumps(out, indent=2))
    print(f"Results written to {out_file}")
    print(f"Summary: W_obs = {W_obs:.2f}")
    print(f"  Asymptotic p: {asymp_p:.5f}")
    print(f"  Wild Cluster Bootstrap: null 95% = {out['wild_cluster_bootstrap']['null_q95']:.2f}, p = {p_wcb:.4f} (SE {p_wcb_mcse:.4f})")
    print(f"  Permutation Test:       null 95% = {out['permutation_test']['null_q95']:.2f}, p = {p_perm:.4f} (SE {p_perm_mcse:.4f})")


if __name__ == "__main__":
    main()
