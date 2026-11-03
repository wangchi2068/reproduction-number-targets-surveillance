"""Conditional, fixed-design calibration. Importing performs no I/O or simulations.

The variance parameter v is WEEKLY MARGINAL variance of the pre-bound additive R
process, not a per-generation variance. National/state innovation variances are
rho*v*(1-phi_national**2) and (1-rho)*v*(1-phi_state**2). At the iid baseline these
coincide with the component marginal variances. All inputs are stated assumptions;
GI sensitivity values are not posterior draws.

Public integration points: build_design, scenario_grid, simulate_panel,
summarize_panel, calibrate_inference, run_calibration. The latter accepts paired
empirical bootstrap draws; it never runs a second empirical bootstrap. CLI output
is mandatory, caches are not read, and existing calibration outputs are refused.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, replace
import hashlib
import importlib
import json
from pathlib import Path
import time
from typing import Mapping

import numpy as np
import pandas as pd
from scipy.special import gammainc, gammaincinv

VGRID = (0.0, 0.02, 0.04, 0.06, 0.09)
GI_CV = 1.6 / 3.6
DEFAULT_SEED = 20261103
SEED_CHUNK_SIZE = 25  # logical seed partition, independent of execution batch size
PRE_WEEKS = 4
METRICS = ("a", "b", "ratio", "xbar", "n", "VR_week", "r1_state", "r1_national")


def _module(name):
    return importlib.import_module(f"{__package__}.{name}" if __package__ else name)


@dataclass(frozen=True)
class Scenario:
    name: str = "matched36_binomial"
    true_gi_mean: float = 3.6
    analyst_gi_mean: float = 3.6
    gi_cv: float = GI_CV
    true_gi_sd: float | None = None
    analyst_gi_sd: float | None = None
    national_share: float = 0.45
    phi_national: float = 0.0
    phi_state: float = 0.0
    observation: str = "binomial"
    observation_phi: float = 1.0
    k: float = 0.35
    rho_h: float = 0.01
    delay_mean: float = 5.5
    delay_sd: float = 3.0
    reporting_log_sd: float = 0.0
    phi_reporting: float = 0.0
    r_lower: float = 0.05

    def __post_init__(self):
        positive = (self.true_gi_mean, self.analyst_gi_mean, self.gi_cv,
                    self.k, self.delay_mean, self.delay_sd, self.r_lower)
        if not all(np.isfinite(x) and x > 0 for x in positive):
            raise ValueError("GI, dispersion, delay and lower R bound must be positive")
        if any(sd is not None and (not np.isfinite(sd) or sd <= 0) for sd in (self.true_gi_sd, self.analyst_gi_sd)):
            raise ValueError("Explicit GI SD must be positive")
        if not 0 <= self.national_share <= 1 or not 0 < self.rho_h <= 1:
            raise ValueError("National share and admission probability must be in [0,1]")
        if not all(np.isfinite(x) and abs(x) < 1 for x in
                   (self.phi_national, self.phi_state, self.phi_reporting)):
            raise ValueError("Persistence must have absolute value below one")
        if self.observation not in ("binomial", "gamma_poisson"):
            raise ValueError("Unknown observation model")
        if not np.isfinite(self.observation_phi) or self.observation_phi < 1:
            raise ValueError("Observation phi is a variance-to-mean ratio >=1")
        if self.observation == "binomial" and self.observation_phi != 1:
            raise ValueError("observation_phi applies only to gamma_poisson")
        if not np.isfinite(self.reporting_log_sd) or self.reporting_log_sd < 0:
            raise ValueError("Reporting log SD must be nonnegative")


@dataclass(frozen=True)
class SeasonDesign:
    season: str
    start: str
    end: str
    weeks: int
    states: tuple[str, ...]
    initial_counts: tuple[float, ...]
    national_growth: float


@dataclass(frozen=True)
class Design:
    seasons: tuple[SeasonDesign, ...]
    raw_hash: str
    missing_initial: tuple[str, ...] = ()

    @property
    def windows(self):
        return {s.season: (pd.Timestamp(s.start), pd.Timestamp(s.end)) for s in self.seasons}


def scenario_grid():
    """Declared small joint grid plus retained, single-input legacy stresses."""
    base = Scenario()
    result = [base, replace(base, name="matched32_binomial", true_gi_mean=3.2, analyst_gi_mean=3.2),
              replace(base, name="matched32_source_sd21_binomial", true_gi_mean=3.2, analyst_gi_mean=3.2,
                      true_gi_sd=2.1, analyst_gi_sd=2.1)]
    for gi in (2.2, 2.9, 4.3):
        result.append(replace(base, name=f"true{gi:g}_analyst36", true_gi_mean=gi))
    for gi in (3.2, 3.6):
        for persistence in (0.0, 0.25):
            for obs_phi in (1.0, 1.25):
                result.append(replace(base, name=f"joint_gi{gi:g}_pn{persistence:g}_obs{obs_phi:g}",
                                      true_gi_mean=gi, analyst_gi_mean=gi,
                                      phi_national=persistence, observation="gamma_poisson",
                                      observation_phi=obs_phi))
    for name, changes in (
        ("national_share0", {"national_share": 0.0}),
        ("national_share65", {"national_share": 0.65}),
        ("national_persistence25", {"phi_national": 0.25}),
        ("state_persistence50", {"phi_state": 0.5}),
        ("delay_short", {"delay_mean": 3.5, "delay_sd": 2.0}),
        ("delay_long", {"delay_mean": 8.0, "delay_sd": 4.0}),
        ("k01", {"k": 0.1}),
        ("rho_h005", {"rho_h": 0.005}),
        ("rho_h02", {"rho_h": 0.02}),
        ("reporting_shock", {"national_share": 0.0, "reporting_log_sd": 0.12}),
    ):
        result.append(replace(base, name=name, **changes))
    return result


def _frame_hash(frame):
    canonical = frame.sort_values([c for c in ("state", "week_end") if c in frame]).copy()
    return hashlib.sha256(canonical.to_csv(index=False, float_format="%.17g").encode()).hexdigest()


def build_design(raw=None, national=None, windows=None):
    """Freeze actual four window durations and observed first-week state counts.

    Missing initial counts are excluded and reported, never filled with zero.
    Observed zeros remain zero (absorbing without imported infections); no floor is
    imposed. These counts set the deterministic initial exponential intensity,
    not a constraint forcing the first random observed week to match exactly.
    """
    data = _module("revision_data")
    raw = data.load_raw(pathogen="flu") if raw is None else raw.copy()
    if raw.duplicated(["state", "week_end"]).any():
        raise ValueError("Duplicate state-week data")
    raw["week_end"] = pd.to_datetime(raw.week_end)
    national = data.national_series(raw) if national is None else national
    if windows is None:
        windows = data.rule_windows(national)[0]
    if len(windows) != 4:
        raise ValueError("Calibration requires the four fixed observed windows")
    seasons, missing = [], []
    states = sorted(raw.state.unique())
    for season, (first, last) in windows.items():
        first, last = pd.Timestamp(first), pd.Timestamp(last)
        dates = pd.date_range(first, last, freq="7D")
        sub = national.set_index("week_end").reindex(dates)
        growth = sub.g.to_numpy(float)
        if len(dates) < 1 or dates[-1] != last or not np.all(np.isfinite(growth) & (growth > 0)):
            raise ValueError(f"Invalid fixed window/growth data for {season}")
        c0 = raw[raw.week_end == first].set_index("state").C.reindex(states)
        valid = c0.notna() & np.isfinite(c0) & (c0 >= 0)
        missing.extend(f"{season}:{s}" for s in c0.index[~valid])
        c0 = c0[valid]
        if c0.empty:
            raise ValueError(f"No initial states for {season}")
        seasons.append(SeasonDesign(str(season), first.date().isoformat(), last.date().isoformat(),
                                    len(dates), tuple(c0.index), tuple(c0.astype(float)),
                                    float(np.exp(np.mean(np.log(growth))))))
    return Design(tuple(seasons), _frame_hash(raw[["state", "week_end", "C", "coverage"]]), tuple(missing))


def discrete_gamma(mean, sd, *, first_lag=1, quantile=0.999):
    if mean <= 0 or sd <= 0 or first_lag not in (0, 1) or not 0 < quantile < 1:
        raise ValueError("Invalid gamma kernel parameters")
    shape, scale = (mean / sd) ** 2, sd ** 2 / mean
    last = int(np.ceil(gammaincinv(shape, quantile) * scale)) + 1
    lags = np.arange(first_lag, last + 1)
    mass = gammainc(shape, (lags + .5) / scale) - gammainc(shape, np.maximum(lags - .5, 0) / scale)
    mass[0] += gammainc(shape, max(first_lag - .5, 0) / scale)
    return mass / mass.sum()


def _words(value):
    return [int(word) for word in np.frombuffer(hashlib.sha256(str(value).encode()).digest()[:16], dtype="<u4")]


def seed_sequence(seed, scenario, v, replicate, *, season=0, stream=0, chunk=None):
    """Stable, independent scenario/cell/logical-chunk/replicate/season/stream key.

    replicate is the global panel index. Logical chunks have 25 panels, so changing
    execution batch size or iteration order cannot change a panel's draws.
    """
    if seed < 0 or replicate < 0:
        raise ValueError("Seed and replicate must be nonnegative")
    logical_chunk, within = divmod(int(replicate), SEED_CHUNK_SIZE)
    if chunk is not None and chunk != logical_chunk:
        raise ValueError("chunk must equal replicate // SEED_CHUNK_SIZE")
    key = json.dumps(asdict(scenario), sort_keys=True) if isinstance(scenario, Scenario) else str(scenario)
    return np.random.SeedSequence([int(seed), 113, 2, *_words(key), *_words(float(v).hex()),
                                   logical_chunk, within, int(season), int(stream)])


def _ar1(rng, shape, phi):
    z = rng.standard_normal(shape)
    if phi:
        for week in range(1, shape[-1]):
            z[..., week] = phi * z[..., week - 1] + np.sqrt(1 - phi * phi) * z[..., week]
    return z


def weekly_r_process(n_states, n_weeks, rbar, v, scenario, national_rng, state_rng):
    """Stationary-start AR(1) shocks; national shock is one draw shared by states."""
    if not np.isfinite(v) or v < 0 or n_states < 1 or n_weeks < 1:
        raise ValueError("Invalid process dimensions or variance")
    eta = np.sqrt(scenario.national_share * v) * _ar1(national_rng, (n_weeks,), scenario.phi_national)
    eps = np.sqrt((1 - scenario.national_share) * v) * _ar1(state_rng, (n_states, n_weeks), scenario.phi_state)
    unbounded = float(rbar) + eta[None, :] + eps
    return {"R": np.maximum(unbounded, scenario.r_lower), "national": eta, "state": eps,
            "clipped": unbounded < scenario.r_lower}


def reporting_multipliers(n_weeks, scenario, rng):
    if scenario.reporting_log_sd == 0:
        return np.ones(n_weeks)
    sd = scenario.reporting_log_sd
    return np.exp(sd * _ar1(rng, (n_weeks,), scenario.phi_reporting) - .5 * sd * sd)


def _intensity(value):
    value = np.asarray(value, float)
    if not np.all(np.isfinite(value)) or np.any(value < 0):
        raise ValueError("Latent intensity must be finite and nonnegative")
    return value


def gamma_poisson_observation(lam, phi, rng):
    """Draw C|lambda with mean lambda and variance phi*lambda (phi >= 1).

    The caller supplies latent delayed intensity, never observed noisy counts.
    Gamma(shape=lambda/(phi-1), scale=phi-1) mixed with Poisson has these moments.
    Zero intensity is identically zero; phi=1 is the exact Poisson limit.
    """
    lam = _intensity(lam)
    if not np.isfinite(phi) or phi < 1:
        raise ValueError("phi must be finite and >=1")
    if phi == 1:
        return rng.poisson(lam)
    mixture = np.zeros_like(lam)
    positive = lam > 0
    mixture[positive] = rng.gamma(lam[positive] / (phi - 1), phi - 1)
    return rng.poisson(mixture)


def delayed_intensity(infections, delay_kernel, rho_h):
    """Full convolution of latent infection trajectory with admission delay."""
    infections = _intensity(infections)
    kernel = _intensity(delay_kernel)
    if infections.ndim < 1 or kernel.ndim != 1 or not len(kernel) or not np.isclose(kernel.sum(), 1):
        raise ValueError("Delay kernel must be a one-dimensional probability vector")
    if not np.isfinite(rho_h) or not 0 <= rho_h <= 1:
        raise ValueError("rho_h must be a probability")
    out = np.zeros(infections.shape[:-1] + (infections.shape[-1] + len(kernel) - 1,))
    for lag, p in enumerate(kernel):
        out[..., lag:lag + infections.shape[-1]] += rho_h * p * infections
    return out


def binomial_delay_observation(infections, delay_kernel, rho_h, rng, *, reporting=None):
    """Allocate each infection to at most one delayed admission, or no admission.

    Sequential conditional binomials are exactly multinomial allocation over delay
    categories. This avoids a state x day x delay tensor. Reporting, if supplied,
    multiplies the admission probability at the *arrival* day and is common across
    states. Infeasible probabilities fail explicitly; none are silently clipped.
    """
    infections = _intensity(infections)
    kernel = _intensity(delay_kernel)
    if infections.ndim < 1 or np.any(infections != np.floor(infections)):
        raise ValueError("Binomial observation requires integer infection counts")
    if kernel.ndim != 1 or not len(kernel) or not np.isclose(kernel.sum(), 1):
        raise ValueError("Delay kernel must sum to one")
    if not np.isfinite(rho_h) or not 0 <= rho_h <= 1:
        raise ValueError("rho_h must be a probability")
    days, out_days = infections.shape[-1], infections.shape[-1] + len(kernel) - 1
    reporting = np.ones(out_days) if reporting is None else _intensity(reporting)
    if reporting.shape != (out_days,):
        raise ValueError("Reporting must have one multiplier per arrival day")
    probabilities = np.zeros(days)
    for lag, p in enumerate(kernel):
        probabilities += rho_h * p * reporting[lag:lag + days]
    if np.any(probabilities > 1 + 1e-12):
        raise ValueError("Reporting shock makes admission probability exceed one")
    remaining = infections.astype(np.int64)
    mass_left = np.ones(days)
    out = np.zeros(infections.shape[:-1] + (out_days,), dtype=np.int64)
    for lag, p in enumerate(kernel):
        prob = rho_h * p * reporting[lag:lag + days]
        conditional = np.divide(prob, mass_left, out=np.zeros(days), where=mass_left > 1e-15)
        allocated = rng.binomial(remaining, np.clip(conditional, 0, 1))
        out[..., lag:lag + days] += allocated
        remaining -= allocated
        mass_left -= prob
    return out


def simulate_season(design, scenario, v, *, seed=DEFAULT_SEED, replicate=0, season_index=0,
                    meanfield=False, retain_latent=False):
    """Daily renewal vectorized over all states; one loop over causal daily time."""
    def rng(stream):
        return np.random.default_rng(seed_sequence(seed, scenario, v, replicate, season=season_index, stream=stream))

    gk = discrete_gamma(scenario.true_gi_mean, scenario.true_gi_sd if scenario.true_gi_sd is not None else scenario.true_gi_mean * scenario.gi_cv)
    dk = discrete_gamma(scenario.delay_mean, scenario.delay_sd, first_lag=0)
    weeks = PRE_WEEKS + design.weeks + 1
    days, history = 7 * weeks, len(gk) + len(dk)
    growth = np.log(design.national_growth) / 7
    # Discrete Euler-Lotka matches the intended daily growth with the actual kernel.
    rbar = 1 / np.sum(gk * np.exp(-growth * np.arange(1, len(gk) + 1)))
    process = weekly_r_process(len(design.states), weeks, rbar, v, scenario, rng(0), rng(1))
    report = reporting_multipliers(weeks, scenario, rng(3))
    report_daily = np.ones(history + days + len(dk) - 1)
    report_daily[history:history + days] = np.repeat(report, 7)
    mgf_delay = np.sum(dk * np.exp(-growth * np.arange(len(dk))))
    scale = scenario.rho_h * mgf_delay * np.exp(7 * growth * PRE_WEEKS) * np.sum(np.exp(growth * np.arange(7)))
    initial = np.asarray(design.initial_counts) / scale
    infections = np.zeros((len(design.states), history + days))
    infections[:, :history] = initial[:, None] * np.exp(growth * np.arange(-history, 0))
    if not meanfield:
        infections[:, :history] = np.rint(infections[:, :history])
    renewal_rng = rng(2)
    for day in range(days):
        pos = history + day
        infectiousness = np.sum(infections[:, pos - len(gk):pos] * gk[::-1], axis=1)
        R = process["R"][:, day // 7]
        if meanfield:
            infections[:, pos] = R * infectiousness
        else:
            mixed = np.zeros(len(design.states))
            active = infectiousness > 0
            mixed[active] = renewal_rng.gamma(scenario.k * infectiousness[active], R[active] / scenario.k)
            infections[:, pos] = renewal_rng.poisson(mixed)
    # No observation has been drawn before this convolution.
    latent_daily = delayed_intensity(infections, dk, scenario.rho_h) * report_daily
    latent_weekly = latent_daily[:, history:history + days].reshape(len(design.states), weeks, 7).sum(axis=2)
    if meanfield:
        counts = latent_weekly
    elif scenario.observation == "binomial":
        observed_daily = binomial_delay_observation(infections, dk, scenario.rho_h, rng(4), reporting=report_daily)
        counts = observed_daily[:, history:history + days].reshape(len(design.states), weeks, 7).sum(axis=2)
    else:
        counts = gamma_poisson_observation(latent_weekly, scenario.observation_phi, rng(4))
    dates = pd.date_range(pd.Timestamp(design.start) - pd.Timedelta(weeks=PRE_WEEKS), periods=weeks, freq="7D")
    raw = pd.DataFrame({"state": np.repeat(design.states, weeks), "week_end": np.tile(dates, len(design.states)),
                        "C": counts.ravel(), "coverage": 100.0, "R_true": process["R"].ravel(),
                        "latent_admission_mean": latent_weekly.ravel(),
                        "reporting_multiplier": np.tile(report, len(design.states))})
    diagnostics = {"season": design.season, "rbar_target": rbar,
                   "n_states": len(design.states), "observation_weeks": design.weeks,
                   "R_lower_bound_count": int(process["clipped"].sum()),
                   "R_week_draws": int(process["R"].size),
                   "zero_initial_states": int(np.sum(np.asarray(design.initial_counts) == 0)),
                   "weekly_marginal_variance_pre_bound": v,
                   "national_weekly_innovation_variance": scenario.national_share * v * (1 - scenario.phi_national ** 2),
                   "state_weekly_innovation_variance": (1 - scenario.national_share) * v * (1 - scenario.phi_state ** 2)}
    result = {"raw": raw, "diagnostics": diagnostics}
    if retain_latent:
        result.update(infections=infections, latent_weekly=latent_weekly, process=process, reporting=report)
    return result


def simulate_panel(design, scenario, v, *, seed=DEFAULT_SEED, replicate=0, chunk=None, meanfield=False):
    if chunk is not None and chunk != replicate // SEED_CHUNK_SIZE:
        raise ValueError("chunk must equal replicate // SEED_CHUNK_SIZE")
    blocks = [simulate_season(s, scenario, v, seed=seed, replicate=replicate, season_index=i,
                              meanfield=meanfield) for i, s in enumerate(design.seasons)]
    return {"raw": pd.concat([b["raw"] for b in blocks], ignore_index=True),
            "diagnostics": [b["diagnostics"] for b in blocks]}


def _r1(values, blocks, dates):
    frame = pd.DataFrame({"x": values, "block": blocks, "date": pd.to_datetime(dates)})
    frame = frame.sort_values(["block", "date"])
    group = frame.groupby("block").x
    mean, count = group.transform("mean"), group.transform("size")
    dev = np.divide(frame.x, mean, out=np.full(len(frame), np.nan), where=np.asarray(mean) != 0) - 1
    valid = (count.to_numpy() >= 3) & np.isfinite(dev)
    frame = frame.loc[valid].copy()
    frame["dev"] = np.asarray(dev)[valid]
    if frame.empty:
        return np.nan, 0
    adjacent = ((frame.block == frame.block.shift(-1)) &
                ((frame.date.shift(-1) - frame.date).dt.days == 7)).to_numpy()
    d = frame.dev.to_numpy()
    denominator = np.sum(d * d)
    numerator = np.sum(d[:-1][adjacent[:-1]] * d[1:][adjacent[:-1]])
    return (float(numerator / denominator) if denominator > 0 else np.nan), int(adjacent.sum())


def summarize_sample(sample):
    """Shared diagnostics for actual and simulated eligible samples.

    VR_week = sample variance of calendar-week means of signed x/xbar-1,
    divided by sample variance of all signed deviations (both ddof=1).
    r1_state uses state-season mean centering and exact consecutive weeks;
    r1_national uses season-centered weekly cross-state means of the same x.
    Each r1 divides the sum of adjacent products by the full centered sum of
    squares within blocks having at least three observations, as declared here.
    """
    stats = _module("revision_stats")
    x, D = sample.x.to_numpy(float), sample.D.to_numpy(float)
    result = {"n": len(sample), "a": np.nan, "b": np.nan, "ratio": np.nan, "xbar": np.nan,
              "converged": False, "VR_week": np.nan, "r1_state": np.nan, "r1_national": np.nan}
    if len(sample) < 3:
        result["fit_status"] = "too_few_eligible_rows"
        return result
    try:
        fit = stats.fit_loss(x, D)
    except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
        result["fit_status"] = f"{type(exc).__name__}: {exc}"
        fit = {}
    result.update(fit)
    for key in ("a", "b", "ratio", "xbar"):
        if result.get(key) is None:
            result[key] = np.nan
    # Preserve signed b/a for calibration diagnostics, including negative a;
    # revision_stats intentionally leaves its operational ratio undefined there.
    if result.get("converged", False) and np.isfinite(result["a"]) and result["a"] != 0:
        result["ratio"] = result["b"] / result["a"]
    result["n"] = len(sample)
    xb = float(np.mean(x))
    result["xbar"] = xb
    dev = x / xb - 1 if xb != 0 else np.full(len(x), np.nan)
    weekly_dev = pd.Series(dev).groupby(sample.week_end.reset_index(drop=True)).mean().to_numpy()
    total = np.var(dev, ddof=1)
    result["VR_week"] = float(np.var(weekly_dev, ddof=1) / total) if total > 0 and len(weekly_dev) > 1 else np.nan
    blocks = sample.state.astype(str) + "|" + sample.season.astype(str)
    result["r1_state"], result["r1_state_pairs"] = _r1(x, blocks, sample.week_end)
    national = sample.groupby(["season", "week_end"], as_index=False).x.mean()
    result["r1_national"], result["r1_national_pairs"] = _r1(national.x, national.season, national.week_end)
    result["a_nonpositive"] = bool(np.isfinite(result["a"]) and result["a"] <= 0)
    result["D_median"] = float(np.median(D))
    return result


def summarize_panel(raw_panel, design, analyst_gi_mean=3.6, analyst_gi_sd=None):
    data = _module("revision_data")
    raw = raw_panel["raw"] if isinstance(raw_panel, dict) else raw_panel
    panel = data.build_panel(raw)
    out = {}
    for stat in ("RA", "G"):
        sample = data.load_sample(stat=stat, gi_mean=analyst_gi_mean, gi_sd=analyst_gi_sd,
                                  panel=panel, windows=design.windows)
        out[stat] = summarize_sample(sample)
        if "R_true" in raw and len(sample):
            matched = sample[["state", "week_end"]].merge(raw[["state", "week_end", "R_true"]],
                                                         on=["state", "week_end"], validate="one_to_one")
            out[stat]["Rbar_true_retained"] = float(matched.R_true.mean())
    return out


def _clean(value):
    if isinstance(value, dict):
        return {str(k): _clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_clean(v) for v in value]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, (pd.Timestamp, Path)):
        return str(value)
    return value


def _write_json(path, value):
    Path(path).write_text(json.dumps(_clean(value), indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")


def _distribution(values):
    values = np.asarray(values, float)
    finite = values[np.isfinite(values)]
    return {"n_total": len(values), "n_finite": len(finite), "n_nonfinite": int((~np.isfinite(values)).sum()),
            "n_nonpositive": int(np.sum(finite <= 0)), "n_negative": int(np.sum(finite < 0)),
            "n_zero": int(np.sum(finite == 0)), "n_positive": int(np.sum(finite > 0)),
            "mean_finite": float(finite.mean()) if len(finite) else np.nan,
            "sd_finite": float(finite.std(ddof=1)) if len(finite) > 1 else np.nan,
            "q025_50_975_finite": np.quantile(finite, [.025, .5, .975]).tolist() if len(finite) else [np.nan] * 3}


def _records_frame(records):
    if isinstance(records, pd.DataFrame):
        frame = records.copy()
    else:
        frame = pd.json_normalize(records)
    return frame.where(pd.notna(frame), np.nan)


def _valid_fits(cell, stat):
    return cell[f"{stat}.converged"].fillna(False).astype(bool).to_numpy() & np.isfinite(cell[f"{stat}.a"].to_numpy(float))


def summarize_cells(records):
    frame = _records_frame(records)
    result = []
    for (setting, v), cell in frame.groupby(["scenario", "v"], sort=True):
        rec = {"scenario": setting, "v": v, "panels": len(cell)}
        for stat in ("RA", "G"):
            valid = _valid_fits(cell, stat)
            rec[stat] = {"converged": int(valid.sum()), "failed": int((~valid).sum())}
            for metric in METRICS:
                values = cell[f"{stat}.{metric}"].to_numpy(float)
                # Fit-dependent statistics cannot promote failed last iterates.
                if metric in ("a", "b", "ratio"):
                    values = np.where(valid, values, np.nan)
                rec[stat][metric] = _distribution(values)
            a = cell.loc[valid, f"{stat}.a"].to_numpy(float)
            b = cell.loc[valid, f"{stat}.b"].to_numpy(float)
            rec[stat]["a_mean_mc_se"] = float(a.std(ddof=1) / np.sqrt(len(a))) if len(a) > 1 else np.nan
            rec[stat]["ratio_of_mean_b_to_mean_a"] = float(b.mean() / a.mean()) if len(a) and a.mean() != 0 else np.nan
        result.append(rec)
    return result


def _gi_key(gi, sd=None):
    return str(float(gi)) if sd is None or not np.isfinite(sd) else f"{float(gi)}_sd{float(sd)}"


def _gi_prefix(gi, sd=None):
    if np.isclose(gi, 3.2):
        return "RA32source" if sd is not None and np.isfinite(sd) and np.isclose(sd, 2.1) else "RA32"
    return "RA"


def _paired_empirical(draws, gi, sd=None):
    frame = pd.DataFrame(draws).copy()
    source = _gi_prefix(gi, sd)
    if "analyst_gi_mean" in frame:
        frame = frame[np.isclose(frame.analyst_gi_mean, gi)].copy()
        if sd is not None and np.isfinite(sd):
            if "analyst_gi_sd" not in frame:
                raise ValueError("Explicit analyst SD needs GI-SD-tagged or RA32source draws")
            frame = frame[np.isclose(frame.analyst_gi_sd, sd)].copy()
        elif "analyst_gi_sd" in frame:
            frame = frame[frame.analyst_gi_sd.isna()].copy()
        if "RA_a" in frame:
            source = "RA"
    if f"{source}_a" not in frame or "G_a" not in frame:
        raise ValueError(f"Paired empirical draws require {source}_a and G_a for this analyst GI")
    out = pd.DataFrame({"RA_a": frame[f"{source}_a"].to_numpy(float), "G_a": frame.G_a.to_numpy(float)})
    for target, origin in (("RA", source), ("G", "G")):
        for metric in ("b", "converged"):
            if f"{origin}_{metric}" in frame:
                out[f"{target}_{metric}"] = frame[f"{origin}_{metric}"].to_numpy()
        if f"{target}_converged" in out:
            out.loc[~out[f"{target}_converged"].fillna(False).astype(bool), f"{target}_a"] = np.nan
    return out


def _empirical_point_map(points, gi, sd=None):
    if isinstance(points, (str, Path)):
        points = json.loads(Path(points).read_text(encoding="utf-8"))
    if "pooled" in points:
        pooled = points["pooled"]
        return {stat: pooled[prefix].get("point", {}) for stat, prefix in
                (("RA", _gi_prefix(gi, sd)), ("G", "G")) if prefix in pooled}
    key = _gi_key(gi, sd)
    if key in points:
        return points[key]
    if sd is None and gi in points:
        return points[gi]
    return points


def _weighted_line(v, y, se):
    valid = np.isfinite(v) & np.isfinite(y) & np.isfinite(se) & (se > 0)
    if valid.sum() < 3 or np.unique(v[valid]).size < 3:
        return np.full(2, np.nan), valid
    X = np.column_stack((np.ones(valid.sum()), v[valid]))
    coefficient = np.linalg.lstsq(X / se[valid, None], y[valid] / se[valid], rcond=None)[0]
    return coefficient, valid


def calibrate_inference(panel_results, paired_empirical_draws, empirical_points=None, *,
                        bootstrap=2000, seed=DEFAULT_SEED):
    """Paired RA/G empirical + simulation uncertainty, without empirical resampling.

    Paired panel indices are bootstrapped within each v cell; each selected panel's
    RA/G values travel together, preserving their Monte Carlo covariance. A shared
    empirical draw index is then combined with the two calibration lines. Failed,
    nonfinite, nonpositive and out-of-grid inversions remain in raw returned draws.
    Percentiles explicitly refer to finite draws and never a positive-only subset.
    Returns {'summary': ..., 'draws': DataFrame}; does not write files.
    """
    if bootstrap < 0:
        raise ValueError("bootstrap must be nonnegative")
    frame = _records_frame(panel_results)
    summaries, draw_records = {}, []
    for scenario, all_cells in frame.groupby("scenario", sort=True):
        gi = float(all_cells.analyst_gi_mean.iloc[0])
        sd = all_cells.analyst_gi_sd.iloc[0] if "analyst_gi_sd" in all_cells else None
        sd = float(sd) if sd is not None and np.isfinite(sd) else None
        empirical = _paired_empirical(paired_empirical_draws, gi, sd)
        if empirical.empty:
            raise ValueError(f"No paired empirical draws for GI {gi}")
        cells = [c.sort_values("replicate") for _, c in all_cells.groupby("v", sort=True)]
        v = np.array([float(c.v.iloc[0]) for c in cells])
        matrices, mean, se = [], np.full((len(cells), 2), np.nan), np.full((len(cells), 2), np.nan)
        for i, cell in enumerate(cells):
            matrix = np.column_stack([np.where(_valid_fits(cell, stat), cell[f"{stat}.a"], np.nan)
                                      for stat in ("RA", "G")]).astype(float)
            matrices.append(matrix)
            for j in range(2):
                valid = matrix[np.isfinite(matrix[:, j]), j]
                if len(valid) >= 2:
                    mean[i, j], se[i, j] = valid.mean(), valid.std(ddof=1) / np.sqrt(len(valid))
        coefficients, rec = [], {"analyst_gi_mean": gi, "analyst_gi_sd": sd,
                                  "v_grid": v.tolist(), "bootstrap": bootstrap,
                                  "empirical_draws_available": len(empirical)}
        point_map = _empirical_point_map(empirical_points or {}, gi, sd)
        empirical_a = []
        for j, stat in enumerate(("RA", "G")):
            coefficient, valid = _weighted_line(v, mean[:, j], se[:, j])
            coefficients.append(coefficient)
            prefix = _gi_prefix(gi, sd) if stat == "RA" else stat
            point = point_map.get(f"{prefix}_a", point_map.get(stat, {}).get("a", np.nan))
            point = float(point) if point is not None else np.nan
            empirical_a.append(point)
            slope = coefficient[1]
            implied = (point - coefficient[0]) / slope if np.isfinite(slope) and slope != 0 else np.nan
            residual = mean[:, j] - coefficient[0] - slope * v
            rec[stat] = {"c0": coefficient[0], "c1": slope, "cell_a_means": mean[:, j].tolist(),
                         "cell_a_mc_se": se[:, j].tolist(), "line_cells_used": int(valid.sum()),
                         "line_status": "ok" if np.all(np.isfinite(coefficient)) else "insufficient_finite_cells",
                         "max_abs_residual_mc_se": float(np.max(np.abs(residual[valid] / se[valid, j]))) if valid.any() and np.isfinite(slope) else np.nan,
                         "empirical_a": point, "v_implied": implied,
                         "point_nonpositive": bool(np.isfinite(implied) and implied <= 0),
                         "point_outside_grid": bool(np.isfinite(implied) and (implied < v.min() or implied > v.max()))}
        for b in range(bootstrap):
            rng = np.random.default_rng(np.random.SeedSequence([seed, 114, *_words(scenario), b]))
            eidx = int(rng.integers(len(empirical)))
            sampled_mean = np.full_like(mean, np.nan)
            for i, matrix in enumerate(matrices):
                sampled = matrix[rng.integers(len(matrix), size=len(matrix))]
                for j in range(2):
                    finite = sampled[np.isfinite(sampled[:, j]), j]
                    if len(finite) >= 2:
                        sampled_mean[i, j] = finite.mean()
            row = {"scenario": scenario, "analyst_gi_mean": gi, "analyst_gi_sd": sd,
                   "draw": b, "empirical_draw": eidx}
            for j, stat in enumerate(("RA", "G")):
                coefficient, _ = _weighted_line(v, sampled_mean[:, j], se[:, j])
                empirical_value = empirical[f"{stat}_a"].iloc[eidx]
                c0, c1 = coefficient
                implied = (empirical_value - c0) / c1 if np.isfinite(c1) and c1 != 0 else np.nan
                point_c0, point_c1 = coefficients[j]
                row.update({f"{stat}_c0": c0, f"{stat}_c1": c1,
                            f"{stat}_a_empirical": empirical_value, f"{stat}_v": implied,
                            f"{stat}_v_mc_only": (empirical_a[j] - c0) / c1 if np.isfinite(c1) and c1 != 0 else np.nan,
                            f"{stat}_v_empirical_only": (empirical_value - point_c0) / point_c1 if np.isfinite(point_c1) and point_c1 != 0 else np.nan,
                            f"{stat}_slope_nonpositive": bool(np.isfinite(c1) and c1 <= 0),
                            f"{stat}_nonpositive": bool(np.isfinite(implied) and implied <= 0),
                            f"{stat}_outside_grid": bool(np.isfinite(implied) and (implied < v.min() or implied > v.max()))})
            row["difference_RA_minus_G"] = row["RA_v"] - row["G_v"]
            row["ratio_RA_over_G"] = row["RA_v"] / row["G_v"] if np.isfinite(row["G_v"]) and row["G_v"] != 0 else np.nan
            draw_records.append(row)
        summary_draws = pd.DataFrame([r for r in draw_records if r["scenario"] == scenario])
        for stat in ("RA", "G"):
            for suffix in ("v", "v_mc_only", "v_empirical_only"):
                rec[stat][suffix + "_draws"] = _distribution(summary_draws[f"{stat}_{suffix}"] if bootstrap else [])
            rec[stat]["n_nonpositive_slopes"] = int(summary_draws[f"{stat}_slope_nonpositive"].sum()) if bootstrap else 0
            rec[stat]["n_outside_grid"] = int(summary_draws[f"{stat}_outside_grid"].sum()) if bootstrap else 0
        for metric in ("difference_RA_minus_G", "ratio_RA_over_G"):
            rec[metric] = _distribution(summary_draws[metric] if bootstrap else [])
        rec["paired_finite_draws"] = int((np.isfinite(summary_draws.RA_v) & np.isfinite(summary_draws.G_v)).sum()) if bootstrap else 0
        rec["response_checks"] = _response_checks(all_cells, rec, point_map)
        summaries[scenario] = rec
    return {"summary": summaries, "draws": pd.DataFrame(draw_records)}


def _response_checks(cells, calibration, points):
    """Checks at nearest simulated v; extrapolated inversion is flagged, not hidden."""
    result = {}
    grid = np.sort(cells.v.unique())
    for stat in ("RA", "G"):
        implied = calibration[stat]["v_implied"]
        if not np.isfinite(implied):
            result[stat] = {"status": "no_finite_point_inversion"}
            continue
        nearest = float(grid[np.argmin(np.abs(grid - implied))])
        cell = cells[cells.v == nearest]
        valid = _valid_fits(cell, stat)
        empirical = points.get(stat, {})
        rec = {"nearest_grid_v": nearest, "point_outside_grid": calibration[stat]["point_outside_grid"],
               "selection": "nearest declared v cell; diagnostic only, no clamping of inversion", "metrics": {}}
        for metric in METRICS:
            values = cell[f"{stat}.{metric}"].to_numpy(float)
            if metric in ("a", "b", "ratio"):
                values = np.where(valid, values, np.nan)
            target = empirical.get(metric, np.nan)
            target = float(target) if target is not None else np.nan
            finite = values[np.isfinite(values)]
            rec["metrics"][metric] = {"simulated": _distribution(values), "empirical": target,
                                      "empirical_percentile_finite": float(np.mean(finite <= target)) if len(finite) and np.isfinite(target) else np.nan}
        result[stat] = rec
    return result


def read_empirical_draws(path):
    path = Path(path)
    if path.suffix == ".npz":
        with np.load(path, allow_pickle=False) as values:
            return {key: values[key] for key in values.files if np.asarray(values[key]).ndim == 1}
    if path.suffix == ".csv":
        return pd.read_csv(path)
    raise ValueError("Empirical draws must be a paired .npz or .csv file")


def _file_hash(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def run_calibration(output_dir, *, replicates=200, replicate_start=0, bootstrap=2000, empirical_draws=None,
                    empirical_points=None, scenarios=None, vgrid=VGRID, seed=DEFAULT_SEED,
                    design=None, progress=True, simulate_only=False):
    """Run fresh simulations and save every panel before any calibration summary.

    No cache is loaded. An existing calibration_panels.jsonl is refused so an
    interrupted run is identifiable and never silently merged with a new config.
    """
    if replicates < 1 or bootstrap < 0 or replicate_start < 0:
        raise ValueError("replicates must be positive; bootstrap and replicate_start nonnegative")
    scenarios = scenario_grid() if scenarios is None else list(scenarios)
    if not scenarios or len({s.name for s in scenarios}) != len(scenarios):
        raise ValueError("Scenarios must be nonempty with unique names")
    if not vgrid or any(not np.isfinite(v) or v < 0 for v in vgrid) or len(set(vgrid)) != len(vgrid):
        raise ValueError("v grid must contain distinct finite nonnegative values")
    output = Path(output_dir)
    path = output / "calibration_panels.jsonl"
    if path.exists():
        raise FileExistsError(f"Fresh output directory required: {path} already exists")
    design = build_design() if design is None else design
    data = _module("revision_data")
    stats = _module("revision_stats")
    if isinstance(empirical_points, (str, Path)):
        empirical_points = json.loads(Path(empirical_points).read_text(encoding="utf-8"))
    if empirical_points is None:
        candidate = output / "empirical.json"
        empirical_points = json.loads(candidate.read_text(encoding="utf-8")) if candidate.exists() else {}
    # Add diagnostics evaluated on the identical empirical transform; supplied
    # pooled point estimates take precedence when available.
    full_points = {}
    observed_panel = data.build_panel()
    for gi, sd in dict.fromkeys((s.analyst_gi_mean, s.analyst_gi_sd) for s in scenarios):
        supplied = _empirical_point_map(empirical_points, gi, sd)
        full_points[_gi_key(gi, sd)] = {}
        for stat in ("RA", "G"):
            diagnostics = summarize_sample(data.load_sample(stat=stat, gi_mean=gi, gi_sd=sd,
                                                            panel=observed_panel, windows=design.windows))
            diagnostics.update(supplied.get(stat, {}))
            prefix = _gi_prefix(gi, sd) if stat == "RA" else stat
            if f"{prefix}_a" in supplied:
                diagnostics["a"] = supplied[f"{prefix}_a"]
            full_points[_gi_key(gi, sd)][stat] = diagnostics
    empirical_points = full_points
    draws_source = "supplied_in_memory" if empirical_draws is not None else None
    if empirical_draws is None:
        candidate = output / "pooled_draws.npz"
        if candidate.exists():
            empirical_draws, draws_source = read_empirical_draws(candidate), str(candidate)
    elif isinstance(empirical_draws, (str, Path)):
        draws_source = str(empirical_draws)
        empirical_draws = read_empirical_draws(empirical_draws)
    # Validate every requested GI before spending simulation time.
    if empirical_draws is not None:
        for gi, sd in dict.fromkeys((s.analyst_gi_mean, s.analyst_gi_sd) for s in scenarios):
            _paired_empirical(empirical_draws, gi, sd)
    code = Path(__file__).resolve()
    assumption_path = code.parents[1] / "input_assumptions.json"
    config = {"scenarios": [asdict(s) for s in scenarios], "v_grid": list(vgrid), "replicates": replicates,
              "replicate_start": replicate_start, "bootstrap": bootstrap, "seed": seed,
              "seed_logical_chunk_size": SEED_CHUNK_SIZE,
              "pre_weeks": PRE_WEEKS, "variance_parameter": "weekly_marginal_pre_bound",
              "simulate_only": simulate_only}
    manifest = {"status": "running", "config": config, "design": asdict(design),
                "empirical_draws_source": draws_source,
                "hashes": {"config": hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest(),
                           "data": design.raw_hash, "calibration_code": _file_hash(code),
                           "data_code": _file_hash(data.__file__), "stats_code": _file_hash(stats.__file__),
                           "input_assumptions": _file_hash(assumption_path) if assumption_path.exists() else None},
                "cache_policy": "none; fresh run only", "numpy_version": np.__version__, "pandas_version": pd.__version__}
    if empirical_draws is not None:
        manifest["hashes"]["empirical_draws"] = _frame_hash(pd.DataFrame(empirical_draws).assign(state=np.arange(len(pd.DataFrame(empirical_draws)))))
    output.mkdir(parents=True, exist_ok=True)
    _write_json(output / "calibration_manifest.json", manifest)
    _write_json(output / "calibration_empirical.json", empirical_points)
    records = []
    start = time.perf_counter()
    with path.open("x", encoding="utf-8") as handle:
        for scenario in scenarios:
            for v in vgrid:
                cell_start = time.perf_counter()
                for rep in range(replicate_start, replicate_start + replicates):
                    result = simulate_panel(design, scenario, float(v), seed=seed, replicate=rep)
                    summary = summarize_panel(result, design, scenario.analyst_gi_mean, scenario.analyst_gi_sd)
                    row = {"scenario": scenario.name, "analyst_gi_mean": scenario.analyst_gi_mean,
                           "analyst_gi_sd": scenario.analyst_gi_sd,
                           "true_gi_mean": scenario.true_gi_mean, "true_gi_sd": scenario.true_gi_sd,
                           "v": float(v), "replicate": rep,
                           "chunk": rep // SEED_CHUNK_SIZE,
                           "seed_entropy": seed_sequence(seed, scenario, v, rep).entropy,
                           "process": result["diagnostics"], **summary}
                    records.append(row)
                    handle.write(json.dumps(_clean(row), allow_nan=False, ensure_ascii=False) + "\n")
                    handle.flush()
                if progress:
                    print(f"{scenario.name} v={v:g}: {replicates} panels, {time.perf_counter()-cell_start:.2f}s", flush=True)
    simulation_seconds = time.perf_counter() - start
    cells = summarize_cells(records)
    _write_json(output / "calibration_cells.json", cells)
    inference = None
    if empirical_draws is not None and not simulate_only:
        inference = calibrate_inference(records, empirical_draws, empirical_points, bootstrap=bootstrap, seed=seed)
        _write_json(output / "calibration_inference.json", inference["summary"])
        inference["draws"].to_csv(output / "calibration_draws.csv", index=False, float_format="%.17g")
    elapsed = time.perf_counter() - start
    manifest.update(status="completed", panels=len(records), simulation_seconds=simulation_seconds,
                    total_seconds=elapsed, seconds_per_panel=simulation_seconds / len(records),
                    calibration_inference="completed" if inference is not None else ("deferred_simulate_only" if simulate_only else "not_run_no_paired_empirical_draws"),
                    full_declared_grid_panels=len(scenario_grid()) * len(VGRID) * 200,
                    full_declared_grid_seconds_at_pilot_rate=simulation_seconds / len(records) * len(scenario_grid()) * len(VGRID) * 200)
    manifest["hashes"]["raw_panel_results"] = _file_hash(path)
    _write_json(output / "calibration_manifest.json", manifest)
    _write_json(output / "calibration.json", {"meta": manifest, "empirical": empirical_points,
                                             "cells": cells, "inference": inference["summary"] if inference else None})
    return {"manifest": manifest, "cells": cells, "inference": inference, "panel_results": records}


def merge_calibration_batches(batch_dirs, output_dir, *, empirical_draws, empirical_points=None,
                              bootstrap=2000, expected_replicates=200, seed=DEFAULT_SEED):
    """Validate disjoint fresh-run batches and infer on the combined panel pool.

    This is explicit aggregation, not resume caching. Source code/data/assumption
    hashes and scenario configurations must agree. Duplicate panel keys and an
    incomplete per-cell replicate range are errors, so a pilot cannot silently
    enter a production calibration.
    """
    output = Path(output_dir)
    if (output / "calibration_panels.jsonl").exists():
        raise FileExistsError("Merged calibration output already exists")
    records, manifests, inputs = [], [], []
    configs = {}
    for directory in batch_dirs:
        directory = Path(directory)
        manifest = json.loads((directory / "calibration_manifest.json").read_text(encoding="utf-8"))
        if manifest.get("status") != "completed":
            raise ValueError(f"Incomplete batch: {directory}")
        path = directory / "calibration_panels.jsonl"
        if _file_hash(path) != manifest["hashes"]["raw_panel_results"]:
            raise ValueError(f"Changed batch contents: {directory}")
        if manifests:
            for key in ("data", "calibration_code", "data_code", "stats_code", "input_assumptions"):
                if manifest["hashes"][key] != manifests[0]["hashes"][key]:
                    raise ValueError(f"Batch fingerprint mismatch: {key}")
            for key in ("seed", "seed_logical_chunk_size", "pre_weeks", "variance_parameter"):
                if manifest["config"][key] != manifests[0]["config"][key]:
                    raise ValueError(f"Batch configuration mismatch: {key}")
        for scenario in manifest["config"]["scenarios"]:
            if scenario["name"] in configs and configs[scenario["name"]] != scenario:
                raise ValueError("Conflicting scenario configurations")
            configs[scenario["name"]] = scenario
        with path.open(encoding="utf-8") as source:
            records.extend(json.loads(line) for line in source if line.strip())
        manifests.append(manifest)
        inputs.append({"directory": str(directory.resolve()), "panels_sha256": _file_hash(path)})
    if not records:
        raise ValueError("No completed calibration batches supplied")
    frame = _records_frame(records)
    if frame.duplicated(["scenario", "v", "replicate"]).any():
        raise ValueError("Duplicate scenario/cell/replicate across batches")
    expected = set(range(expected_replicates))
    for scenario in configs:
        part = frame[frame.scenario == scenario]
        if set(part.v) != set(VGRID):
            raise ValueError(f"Incomplete variance grid for {scenario}")
        for v, cell in part.groupby("v"):
            if set(cell.replicate) != expected:
                raise ValueError(f"Incomplete replicate range for {scenario}, v={v}")
    if isinstance(empirical_draws, (str, Path)):
        empirical_draws = read_empirical_draws(empirical_draws)
    combined_points = {}
    for directory in batch_dirs:
        combined_points.update(json.loads((Path(directory) / "calibration_empirical.json").read_text(encoding="utf-8")))
    if empirical_points is not None:
        for scenario in configs.values():
            gi, sd = scenario["analyst_gi_mean"], scenario.get("analyst_gi_sd")
            supplied = _empirical_point_map(empirical_points, gi, sd)
            for stat in ("RA", "G"):
                combined_points.setdefault(_gi_key(gi, sd), {}).setdefault(stat, {}).update(supplied.get(stat, {}))
    empirical_points = combined_points
    result = calibrate_inference(records, empirical_draws, empirical_points, bootstrap=bootstrap, seed=seed)
    output.mkdir(parents=True, exist_ok=True)
    with (output / "calibration_panels.jsonl").open("x", encoding="utf-8") as target:
        for row in sorted(records, key=lambda r: (r["scenario"], r["v"], r["replicate"])):
            target.write(json.dumps(_clean(row), allow_nan=False) + "\n")
    _write_json(output / "calibration_cells.json", summarize_cells(records))
    _write_json(output / "calibration_inference.json", result["summary"])
    result["draws"].to_csv(output / "calibration_draws.csv", index=False, float_format="%.17g")
    manifest = {"status": "completed", "aggregation": "validated disjoint batch union", "inputs": inputs,
                "panels": len(records), "replicates_per_cell": expected_replicates, "bootstrap": bootstrap,
                "source_hashes": manifests[0]["hashes"], "aggregation_code_sha256": _file_hash(__file__),
                "scenario_configurations": list(configs.values()),
                "hashes": {"raw_panel_results": _file_hash(output / "calibration_panels.jsonl"),
                           "empirical_draws": _frame_hash(pd.DataFrame(empirical_draws).assign(state=np.arange(len(pd.DataFrame(empirical_draws)))))} }
    _write_json(output / "calibration_manifest.json", manifest)
    _write_json(output / "calibration_empirical.json", empirical_points)
    _write_json(output / "calibration.json", {"meta": manifest, "empirical": empirical_points,
                                             "cells": summarize_cells(records), "inference": result["summary"]})
    return {"manifest": manifest, **result}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--replicate-start", type=int, default=0, help="First global replicate index for independent execution batches")
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--simulate-only", action="store_true", help="Save panel results and defer paired calibration to batch merge")
    parser.add_argument("--empirical-draws", type=Path, help="Paired revised pooled .npz/.csv; auto-detect output-dir/pooled_draws.npz")
    parser.add_argument("--empirical-points", type=Path, help="Revised empirical.json or GI-keyed point summaries")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--scenarios", nargs="+", help="Declared scenario names; omit for full grid")
    parser.add_argument("--v-grid", nargs="+", type=float, default=VGRID)
    parser.add_argument("--merge-batches", nargs="+", type=Path, help="Validate and combine completed batch directories")
    args = parser.parse_args(argv)
    if args.merge_batches:
        if args.empirical_draws is None:
            parser.error("--merge-batches requires --empirical-draws")
        merge_calibration_batches(args.merge_batches, args.output_dir, empirical_draws=args.empirical_draws,
                                  empirical_points=args.empirical_points, bootstrap=args.bootstrap,
                                  expected_replicates=args.replicates, seed=args.seed)
        return 0
    scenarios = scenario_grid()
    if args.scenarios:
        available = {s.name: s for s in scenarios}
        missing = set(args.scenarios) - available.keys()
        if missing:
            parser.error(f"Unknown scenarios {sorted(missing)}; available: {list(available)}")
        scenarios = [available[name] for name in args.scenarios]
    run_calibration(args.output_dir, replicates=args.replicates, replicate_start=args.replicate_start, bootstrap=args.bootstrap,
                    empirical_draws=args.empirical_draws, empirical_points=args.empirical_points,
                    scenarios=scenarios, vgrid=args.v_grid, seed=args.seed, simulate_only=args.simulate_only)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
