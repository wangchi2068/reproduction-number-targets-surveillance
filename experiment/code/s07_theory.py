"""s07: closed-form theory numbers and curve blocks for Figs 1-2 (mode "analytic"), plus a small
seeded generation-level Monte Carlo for the simulation points of Fig 1B-C (mode "simulated").

Provenance: formulas follow the manuscript (Propositions 1-7, Corollary 2, Remark 1, Eq. 25,
Supplement Eq. S16) and re-use the logic of the review verification script
`.work/case-count-bottleneck/review/checks/theory_checks.py` (A, m_x, v_H, A_H, m_x_H, lever terms,
B_dagger root of g_T = 1, trigamma constant) and of the original Monte Carlo
`experiment/original/code/r4_sim_core.py` (gamma-Poisson representation of NB offspring).

Model (one generation, H1-H4): R_t ~ Gamma(mean Rbar, var v_R); infections I | R_t ~ NB(n k, k/(k+R_t))
(sum of n NB(R_t, k) offspring counts, drawn as Poisson(Gamma(n k, R_t/k))); reported C ~ Bin(I, rho);
Rhat = C / (rho n), m = rho n.
  V(n, rho) = v_R + Rbar/(rho n) + (Rbar^2 + v_R)/(k n) = v_R + A/m,  A = Rbar + rho (Rbar^2 + v_R)/k.
  Fixed n, rho -> 1: V -> v_R + Rbar/n + (Rbar^2 + v_R)/(k n)  (Remark 1).
  Reporting heterogeneity (Prop 3): v_H = tau^2 (Rbar^2+v_R), A_H = Rbar + rho (1+tau^2)(Rbar^2+v_R)/k,
  m_x_H = A_H / (v_R + v_H).
  Budget (Prop 4): m* = sqrt(c0/c1 * m_x).  Prop 7: e(T) = [T(1+phi)/(1-phi) - 2phi(1-phi^T)/(1-phi)^2]/T^2,
  g_T = -T^2 e'(T); T_dagger solves g_T = 1; B_dagger = T_dagger (c0 + c1 m0*), in units of c1.
  Levers (Prop 6, T periods, AR(1) env phi, AR(1) reporting psi):
    env = v_R S(phi,T)/T^2, rep = tau^2 [Rbar^2 S(psi,T) + v_R S(phi psi,T)]/T^2, samp = A_H/(T m),
    S(phi,T) = sum_{h=-(T-1)}^{T-1} (T-|h|) phi^|h|.
  Prop 5(ii): RSE >= [T xi (xi psi'(xi) - 1)]^(-1/2), xi = Rbar^2/v_R (Supplement Eq. S16; T = number
  of periods).  Corollary 2: RSE >= sqrt(2/(S-1)) with the mean estimated.
"""
import numpy as np
from scipy import optimize, special

from cc_common import SEED, Log, Part

log = Log("s07_theory")
P = Part("s07_theory", "s07_theory.py", ["closed form; manuscript Props 1-7"], "analytic", SEED + 700)

RB, VR, K, RHO = 1.15, 0.04, 0.35, 0.25
TAU = 0.15


def A_of(rb=RB, vr=VR, k=K, rho=RHO):
    return rb + rho * (rb ** 2 + vr) / k


def V_m(m, rb=RB, vr=VR, k=K, rho=RHO):
    return vr + A_of(rb, vr, k, rho) / np.asarray(m, float)


def V_nrho(n, rho, rb=RB, vr=VR, k=K):
    return vr + rb / (rho * n) + (rb ** 2 + vr) / (k * n)


def S_ar(phi, T):
    h = np.arange(-(T - 1), T)
    return float(np.sum((T - np.abs(h)) * phi ** np.abs(h)))


def S_ar_cont(phi, T):
    """Closed form of S(phi, T), valid for real T (manuscript Eq. 16)."""
    T = np.asarray(T, float)
    if phi == 0:
        return T
    return T * (1 + phi) / (1 - phi) - 2 * phi * (1 - phi ** T) / (1 - phi) ** 2


def levers(T, m, tau, phi=0.5, psi=0.5, rb=RB, vr=VR, k=K, rho=RHO):
    env = vr * S_ar_cont(phi, T) / T ** 2
    rep = tau ** 2 * (rb ** 2 * S_ar_cont(psi, T) + vr * S_ar_cont(phi * psi, T)) / T ** 2
    AH = rb + rho * (1 + tau ** 2) * (rb ** 2 + vr) / k
    samp = AH / (T * m)
    return {"env": float(env), "rep": float(rep), "samp": float(samp), "total": float(env + rep + samp)}


A = A_of(); MX = A / VR
P.num("A", A, note="A = Rbar + rho(Rbar^2+v_R)/k at Rbar=1.15, v_R=0.04, k=0.35, rho=0.25")
P.num("m_cross", MX, unit="reported cases", note="m_x = A/v_R; manuscript 53.1")
for m in (10, MX, 200):
    P.num(f"doubling_reduction_m{m:.0f}", 1 / (2 * (1 + m / MX)), note="Prop 2(iii): 1/[2(1+m/m_x)]")
P.num("n_over_rho_lever_ratio", A / RB, note="Remark 1: A/Rbar, relative gain of n vs rho per % increase")

# (n, rho) form and the fixed-n floor (Prop 1-2, Remark 1)
for n in (40, 200, 1000):
    P.num(f"V_floor_rho1_n{n}", float(V_nrho(n, 1.0)),
          note="fixed n, rho->1: v_R + Rbar/n + (Rbar^2+v_R)/(k n) (Remark 1); exceeds v_R=0.04")

# Reporting heterogeneity (Prop 3)
vH = TAU ** 2 * (RB ** 2 + VR)
AH = RB + RHO * (1 + TAU ** 2) * (RB ** 2 + VR) / K
P.num("v_H_tau015", vH, note="v_H = tau^2(Rbar^2+v_R); manuscript 0.0307")
P.num("A_H_tau015", AH)
P.num("m_cross_H_tau015", AH / (VR + vH), note="m_x_H = A_H/(v_R+v_H) (Prop 3(ii)); manuscript 30.4")

# Budget: m* and B_dagger (Props 4 and 7)
for mx in (21.3, 49.2, MX):
    for c in (10, 100):
        P.num(f"m_star_mx{mx:.2f}_c{c}", float(np.sqrt(c * mx)), note="m* = sqrt(c0/c1 * m_x)")


def e_T(T, ph):
    return S_ar_cont(ph, T) / np.asarray(T, float) ** 2


def g_T(T, ph, h=1e-6):
    return -T ** 2 * (e_T(T + h, ph) - e_T(T - h, ph)) / (2 * h)


def bdag(ph, c0c1, mx):
    Td = optimize.brentq(lambda t: g_T(t, ph) - 1, 1.0001, 500)
    m0 = np.sqrt(c0c1 * mx)
    return float(Td), float(m0), float(Td * (c0c1 + m0))


Td, m0, B = bdag(0.8, 100, MX)
P.num("T_dagger_phi08", Td, note="root of g_T = 1 at phi=0.8; manuscript 4.48")
P.num("m0_star_baseline_c100", m0, note="manuscript 72.9")
P.num("B_dagger_baseline", B, unit="c1 units", note="phi=0.8, c0/c1=100; manuscript 775")
P.num("B_dagger_mx21.3", bdag(0.8, 100, 21.3)[2], unit="c1 units", note="manuscript range 655-763")
P.num("B_dagger_mx49.2", bdag(0.8, 100, 49.2)[2], unit="c1 units", note="manuscript range 655-763")

# Levers (Prop 6), baseline T=10, m=50, tau=0.15, phi=psi=0.5, units 1e-3
lv = {"base": levers(10, 50, TAU), "m_x10": levers(10, 500, TAU), "T_x4": levers(40, 50, TAU),
      "tau0": levers(10, 50, 0.0)}
for nm, d in lv.items():
    for kk, v in d.items():
        P.num(f"lever_{nm}_{kk}", v * 1e3, unit="1e-3",
              note=("manuscript reports reporting term 7.85 and total 22.5; closed form gives 7.88 and 22.57 "
                    "(rounding/implementation difference in the manuscript, not reproduced)")
              if (nm == "base" and kk in ("rep", "total")) else None)

# Corollary 2 and Prop 5(ii)
for S_ in (3, 4):
    P.num(f"rse_season_S{S_}", float(np.sqrt(2 / (S_ - 1))),
          note=f"sqrt(2/(S-1)), mean estimated; sqrt(2/S) = {np.sqrt(2 / S_):.3f} (the manuscript's 0.82/0.71 "
               "for S=3/4 are sqrt(2/S))")
xi = RB ** 2 / VR
c505 = float(xi * (xi * special.polygamma(1, xi) - 1))
P.num("xi", xi, note="xi = Rbar^2/v_R")
P.num("fisher_const_0505", c505, note="xi(xi psi'(xi) - 1); RSE >= [T * const]^(-1/2), T = number of periods (Supplement Eq. S16)")
for T in (3, 10, 50):
    P.num(f"rse_vR_T{T}", float((c505 * T) ** -0.5), note="Prop 5(ii); Gaussian limit sqrt(2/T)")

# ----------------------------------------------------------------------------- curves (Fig 1)
mgrid = np.logspace(0, 3, 200)
P.put("fig1_V_m", {"m": mgrid, "V": V_m(mgrid), "target1": A / mgrid, "v_R": VR, "m_cross": MX, "A": A})
rgrid = np.linspace(0.02, 1.0, 200)
P.put("fig1_V_rho", {"rho": rgrid, "n": [40, 200, 1000],
                     "V": {str(n): V_nrho(n, rgrid) for n in (40, 200, 1000)},
                     "floor_rho1": {str(n): float(V_nrho(n, 1.0)) for n in (40, 200, 1000)}, "v_R": VR})
vg = np.logspace(np.log10(0.005), np.log10(0.32), 120)
kg = np.logspace(np.log10(0.05), np.log10(5), 120)
VV, KK = np.meshgrid(vg, kg)
P.put("fig1_mx_contour", {"v_R": vg, "k": kg, "m_cross": (RB + RHO * (RB ** 2 + VV) / KK) / VV,
                          "rows": "k", "cols": "v_R", "Rbar": RB, "rho": RHO,
                          "baseline": {"v_R": VR, "k": K, "m_cross": MX},
                          "empirical_scenarios": {"Rbar": 1.6706, "v_R": [0.04, 0.08], "k": [0.2, 0.5]}})

# ----------------------------------------------------------------------------- curves (Fig 2)
taus = [0.0, 0.1, 0.15, 0.2]
blk = {"m": mgrid, "tau": taus, "V": {}, "floor": {}, "m_cross_H": {}}
for t in taus:
    vh = t ** 2 * (RB ** 2 + VR); ah = RB + RHO * (1 + t ** 2) * (RB ** 2 + VR) / K
    blk["V"][str(t)] = VR + vh + ah / mgrid
    blk["floor"][str(t)] = VR + vh
    blk["m_cross_H"][str(t)] = ah / (VR + vh)
P.put("fig2_VH_m", blk)
tg = np.linspace(0, 0.4, 161)
P.put("fig2_mxH_tau", {"tau": tg, "m_cross_H": (RB + RHO * (1 + tg ** 2) * (RB ** 2 + VR) / K) / (VR + tg ** 2 * (RB ** 2 + VR)),
                       "m_cross": MX})
cg = np.logspace(0, 3, 100)
P.put("fig2_mstar", {"c0c1": cg, "m_x": [21.3, 49.2], "m_star": {"21.3": np.sqrt(cg * 21.3), "49.2": np.sqrt(cg * 49.2)}})

# Equal-budget lever comparison (Fig 2D). Budget B = T (c0 + c1 m), c0/c1 = 100, baseline T=10, m=50,
# so B0 = 1500 c1. Volume strategy: T fixed at 10, m = B/(10 c1) - c0/c1. Period strategy: m fixed at
# 50, T = B/(c0 + 50 c1) (real-valued T via the closed form of S). Each at tau = 0.15 and tau = 0
# (reporting standardisation). MSE from Prop 6 with phi = psi = 0.5.
C0C1 = 100.0; T0, M0 = 10, 50
B0 = T0 * (C0C1 + M0)
mult = np.linspace(1, 8, 57)
eq = {"budget_multiple": mult, "B0": B0, "c0c1": C0C1, "T0": T0, "m0": M0, "phi": 0.5, "psi": 0.5, "series": {}}
for tau in (TAU, 0.0):
    vol = [levers(T0, mu * B0 / T0 - C0C1, tau)["total"] * 1e3 for mu in mult]
    per = [levers(mu * B0 / (C0C1 + M0), M0, tau)["total"] * 1e3 for mu in mult]
    eq["series"][f"volume_tau{tau}"] = vol
    eq["series"][f"periods_tau{tau}"] = per
P.put("fig2_equal_budget", eq)
P.num("eqbudget_x2_volume", eq["series"][f"volume_tau{TAU}"][8], unit="1e-3",
      note="MSE at 2x baseline budget, extra budget spent on volume (T=10, m=250)")
P.num("eqbudget_x2_periods", eq["series"][f"periods_tau{TAU}"][8], unit="1e-3",
      note="MSE at 2x baseline budget, extra budget spent on periods (T=20, m=50)")
P.num("eqbudget_x1_standardised", eq["series"]["volume_tau0.0"][0], unit="1e-3",
      note="MSE at baseline budget with tau=0 (reporting standardisation)")

# ----------------------------------------------------------------------------- Monte Carlo points
rng = np.random.default_rng(SEED + 700)
NREP = 20000


def sim(n, rho, reps=NREP):
    sh = RB ** 2 / VR
    R = rng.gamma(sh, RB / sh, reps)
    I = rng.poisson(rng.gamma(n * K, R / K))
    C = rng.binomial(I, rho)
    Rh = C / (rho * n)
    e2, e1 = (Rh - RB) ** 2, (Rh - R) ** 2

    def ci(x):
        se = x.std(ddof=1) / np.sqrt(len(x))
        return float(x.mean()), [float(x.mean() - 1.96 * se), float(x.mean() + 1.96 * se)]
    return ci(e2), ci(e1)


ms = np.round(np.logspace(0, 3, 10), 2)
pts = {"m": ms, "V_sim": [], "V_ci": [], "T1_sim": [], "T1_ci": [], "reps": NREP, "rho": RHO}
for m in ms:
    (v, vc), (t1, tc) = sim(m / RHO, RHO)
    pts["V_sim"].append(v); pts["V_ci"].append(vc); pts["T1_sim"].append(t1); pts["T1_ci"].append(tc)
    log(f"m={m:8.2f} V_sim={v:.5f} theory={V_m(m):.5f} T1={t1:.5f} theory={A / m:.5f}")
P.put("fig1_sim_V_m", pts)
rel = np.array(pts["V_sim"]) / V_m(ms) - 1
P.num("mc_max_rel_dev_V", float(np.max(np.abs(rel))), mode="simulated", seed=SEED + 700,
      note=f"max |V_sim/V_theory - 1| over 10 m values, {NREP} reps each")
rp = {"rho": [0.1, 0.25, 0.5, 0.75, 1.0], "n": [40, 200, 1000], "V_sim": {}, "V_ci": {}, "reps": NREP}
for n in rp["n"]:
    rp["V_sim"][str(n)] = []; rp["V_ci"][str(n)] = []
    for r in rp["rho"]:
        (v, vc), _ = sim(n, r)
        rp["V_sim"][str(n)].append(v); rp["V_ci"][str(n)].append(vc)
    log(n, rp["V_sim"][str(n)])
P.put("fig1_sim_V_rho", rp)
P.num("mc_V_n40_rho1", rp["V_sim"]["40"][-1], mode="simulated", seed=SEED + 700,
      note=f"theory {V_nrho(40, 1.0):.4f}")
log.close()
