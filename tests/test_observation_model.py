"""Moment, limit, pairing, and import contracts for corrected calibration."""
from dataclasses import replace
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

CODE = Path(__file__).resolve().parents[1] / "experiment" / "code"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))
import revision_calibration as cal


@pytest.mark.parametrize("phi", [1.0, 1.25, 2.0])
@pytest.mark.parametrize("mean", [0.2, 3.0, 50.0])
def test_gamma_poisson_conditional_moments(mean, phi):
    n = 120_000
    counts = cal.gamma_poisson_observation(np.full(n, mean), phi, np.random.default_rng(98731))
    variance = phi * mean
    # Six standard errors for mean and variance; NB fourth cumulant is
    # lambda*phi*(6*phi**2 - 6*phi + 1).
    mean_se = np.sqrt(variance / n)
    fourth_cumulant = mean * phi * (6 * phi * phi - 6 * phi + 1)
    variance_se = np.sqrt((fourth_cumulant + 2 * variance * variance) / n)
    assert abs(counts.mean() - mean) < 6 * mean_se
    assert abs(counts.var(ddof=1) - variance) < 6 * variance_se
    assert np.issubdtype(counts.dtype, np.integer)


def test_poisson_exact_limit_and_zero_intensity():
    lam = np.array([0., .2, 1., 20., 1000.])
    observed = cal.gamma_poisson_observation(lam, 1, np.random.default_rng(7))
    expected = np.random.default_rng(7).poisson(lam)
    np.testing.assert_array_equal(observed, expected)
    np.testing.assert_array_equal(cal.gamma_poisson_observation(np.zeros((12, 4)), 1.25, np.random.default_rng(7)), 0)


@pytest.mark.parametrize("phi", [.9, np.nan, np.inf])
def test_observation_rejects_invalid_phi(phi):
    with pytest.raises(ValueError):
        cal.gamma_poisson_observation([1], phi, np.random.default_rng(8))


@pytest.mark.parametrize("lam", [[-1.], [np.nan], [np.inf]])
def test_observation_rejects_invalid_intensity(lam):
    with pytest.raises(ValueError):
        cal.gamma_poisson_observation(lam, 1, np.random.default_rng(8))


def test_delay_binomial_moments_and_allocation_covariance():
    n, infections, rho = 100_000, 80, .08
    delay = np.array([.2, .5, .3])
    trajectory = np.full((n, 1), infections)
    counts = cal.binomial_delay_observation(trajectory, delay, rho, np.random.default_rng(713))
    total = counts.sum(axis=1)
    assert np.all(total <= infections)
    assert abs(total.mean() - infections * rho) < .035
    assert abs(total.var(ddof=1) - infections * rho * (1-rho)) < .12
    np.testing.assert_allclose(counts.mean(axis=0), infections * rho * delay, atol=.025)
    covariance = np.cov(counts.T)
    expected = -infections * (rho * delay[0]) * (rho * delay[1])
    assert abs(covariance[0, 1] - expected) < .025
    np.testing.assert_allclose(cal.delayed_intensity([[infections]], delay, rho),
                               [[infections * rho * p for p in delay]])


def test_delay_allocation_limits_and_reporting_mean():
    infections = np.array([[3, 0, 7], [1, 4, 2]])
    np.testing.assert_array_equal(cal.binomial_delay_observation(infections, [1.], 1, np.random.default_rng(2)), infections)
    np.testing.assert_array_equal(cal.binomial_delay_observation(infections, [0., 1.], 1, np.random.default_rng(2)),
                                  np.pad(infections, ((0, 0), (1, 0))))
    np.testing.assert_array_equal(cal.binomial_delay_observation(infections, [1.], 0, np.random.default_rng(2)), 0)
    counts = cal.binomial_delay_observation(np.full((100_000, 1), 40), [.3, .7], .1,
                                            np.random.default_rng(19), reporting=[.8, 1.2])
    np.testing.assert_allclose(counts.mean(axis=0), [40*.1*.3*.8, 40*.1*.7*1.2], atol=.03)
    with pytest.raises(ValueError, match="exceed one"):
        cal.binomial_delay_observation(infections, [1.], .9, np.random.default_rng(2), reporting=np.full(3, 2.))


def _tiny_season():
    return cal.SeasonDesign("2022-23", "2022-10-01", "2022-10-15", 3,
                            ("AA", "BB", "CC"), (50., 100., 0.), 1.12)


def test_alternative_uses_latent_convolution_before_observation(monkeypatch):
    scenario = replace(cal.Scenario(), observation="gamma_poisson", observation_phi=1.25)
    captured = []
    original = cal.gamma_poisson_observation
    def capture(lam, phi, rng):
        captured.append(lam.copy())
        return original(lam, phi, rng)
    def forbidden(*args, **kwargs):
        raise AssertionError("No admission draw may precede gamma-Poisson observation")
    monkeypatch.setattr(cal, "gamma_poisson_observation", capture)
    monkeypatch.setattr(cal, "binomial_delay_observation", forbidden)
    result = cal.simulate_season(_tiny_season(), scenario, .04, retain_latent=True)
    np.testing.assert_array_equal(captured[0], result["latent_weekly"])
    kernel = cal.discrete_gamma(scenario.delay_mean, scenario.delay_sd, first_lag=0)
    gkernel = cal.discrete_gamma(scenario.true_gi_mean, scenario.true_gi_mean*scenario.gi_cv)
    history = len(kernel)+len(gkernel)
    daily = cal.delayed_intensity(result["infections"], kernel, scenario.rho_h)
    expected = daily[:, history:history+cal.PRE_WEEKS*7+4*7].reshape(3, 8, 7).sum(axis=2)
    np.testing.assert_allclose(captured[0], expected)
    np.testing.assert_array_equal(result["latent_weekly"][2], 0)


def test_common_shocks_reproducible_and_separate_from_state_noise():
    scenario = replace(cal.Scenario(), national_share=1., phi_national=.25)
    def draw(n, rep=0):
        national = np.random.default_rng(cal.seed_sequence(123, scenario, .04, rep, stream=0))
        state = np.random.default_rng(cal.seed_sequence(123, scenario, .04, rep, stream=1))
        return cal.weekly_r_process(n, 40, 1.2, .04, scenario, national, state)
    first, repeat, larger, independent = draw(3), draw(3), draw(8), draw(3, 1)
    np.testing.assert_array_equal(first["R"], repeat["R"])
    np.testing.assert_array_equal(first["national"], larger["national"])
    np.testing.assert_array_equal(first["R"][0], first["R"][1])
    assert not np.array_equal(first["national"], independent["national"])
    np.testing.assert_array_equal(first["state"], 0)


def test_weekly_marginal_variance_separate_persistence_and_zero_limit():
    scenario = replace(cal.Scenario(), national_share=.4, phi_national=.25, phi_state=.6)
    process = cal.weekly_r_process(4, 80_000, 1.2, .04, scenario,
                                   np.random.default_rng(31), np.random.default_rng(32))
    assert process["national"].var() == pytest.approx(.04*.4, rel=.03)
    assert process["state"].var() == pytest.approx(.04*.6, rel=.02)
    assert np.corrcoef(process["national"][:-1], process["national"][1:])[0,1] == pytest.approx(.25, abs=.02)
    assert np.corrcoef(process["state"][0,:-1], process["state"][0,1:])[0,1] == pytest.approx(.6, abs=.02)
    zero = cal.weekly_r_process(3, 10, 1.2, 0, scenario, np.random.default_rng(1), np.random.default_rng(2))
    np.testing.assert_array_equal(zero["R"], 1.2)
    np.testing.assert_array_equal(zero["clipped"], False)


def test_reporting_is_synchronous_independent_stream_and_mean_one():
    scenario = replace(cal.Scenario(), reporting_log_sd=.12)
    seed = cal.seed_sequence(19, scenario, .04, 0, stream=3)
    a = cal.reporting_multipliers(100_000, scenario, np.random.default_rng(seed))
    b = cal.reporting_multipliers(100_000, scenario, np.random.default_rng(seed))
    np.testing.assert_array_equal(a, b)
    assert a.mean() == pytest.approx(1., abs=.002)
    np.testing.assert_array_equal(cal.reporting_multipliers(10, cal.Scenario(), np.random.default_rng(5)), 1)


def test_seed_partition_independent_of_execution_order():
    scenario = cal.Scenario()
    keys = [cal.seed_sequence(1, scenario, .04, r).generate_state(4).tolist() for r in (0, 1, 24, 25, 26)]
    assert len(set(map(tuple, keys))) == 5
    np.testing.assert_array_equal(cal.seed_sequence(1, scenario, .04, 25, chunk=1).generate_state(4), keys[3])
    with pytest.raises(ValueError):
        cal.seed_sequence(1, scenario, .04, 25, chunk=0)
    first = cal.simulate_season(_tiny_season(), scenario, .04, replicate=25)
    cal.simulate_season(_tiny_season(), scenario, .04, replicate=1)
    second = cal.simulate_season(_tiny_season(), scenario, .04, replicate=25)
    pd.testing.assert_frame_equal(first["raw"], second["raw"])


def test_matched_meanfield_growth_uses_actual_discrete_kernel():
    design = _tiny_season()
    scenario = replace(cal.Scenario(), true_gi_mean=3.2, true_gi_sd=2.1)
    result = cal.simulate_season(design, scenario, 0, meanfield=True)
    counts = result["raw"].C.to_numpy().reshape(3, -1)
    np.testing.assert_allclose(counts[:2, cal.PRE_WEEKS], design.initial_counts[:2], rtol=1e-12)
    np.testing.assert_allclose(counts[:2,1:] / counts[:2,:-1], design.national_growth, rtol=1e-12)


def test_shared_empirical_transform_and_fits():
    import revision_data as data
    import revision_stats as stats
    design = cal.build_design()
    assert [s.weeks for s in design.seasons] == [9, 12, 8, 13]
    assert all(len(s.states) == 51 for s in design.seasons)
    raw = cal.simulate_panel(design, cal.Scenario(), .04)["raw"]
    summary = cal.summarize_panel(raw, design)
    for stat in ("RA", "G"):
        sample = data.load_sample(stat, raw=raw, windows=design.windows)
        fit = stats.fit_loss(sample.x, sample.D)
        assert summary[stat]["n"] == len(sample)
        assert summary[stat]["converged"] == fit["converged"]
        if fit["converged"]:
            assert summary[stat]["a"] == fit["a"]
            assert summary[stat]["b"] == fit["b"]
        assert summary[stat]["r1_state_pairs"] > 0
        assert summary[stat]["r1_national_pairs"] == 38


def _paired_panels():
    rows = []
    for v in cal.VGRID:
        for rep in range(12):
            a = .01 + .8*v + (rep - 5.5)*.0003
            record = {"scenario": "paired", "analyst_gi_mean": 3.6, "v": v, "replicate": rep}
            for stat, factor in (("RA", 1), ("G", 2)):
                record[stat] = {"a": factor*a, "b": factor*.2, "ratio": .2/a, "xbar": 1.,
                                "n": 100, "VR_week": .3, "r1_state": -.1, "r1_national": .2,
                                "converged": True}
            rows.append(record)
    return rows


def test_paired_calibration_preserves_mc_covariance_and_boundary_draws():
    draws = np.linspace(-.01, .11, 101)
    empirical = {"RA_a": draws, "G_a": 2*draws}
    points = {"RA": {"a": .04}, "G": {"a": .08}}
    result = cal.calibrate_inference(_paired_panels(), empirical, points, bootstrap=150, seed=723)
    raw = result["draws"]
    np.testing.assert_allclose(raw.RA_v, raw.G_v, atol=1e-12)
    assert (raw.RA_v < 0).any()
    assert (raw.RA_v > max(cal.VGRID)).any()
    assert result["summary"]["paired"]["RA"]["v_draws"]["n_nonpositive"] > 0
    again = cal.calibrate_inference(_paired_panels(), empirical, points, bootstrap=150, seed=723)
    pd.testing.assert_frame_equal(raw, again["draws"])


def test_failed_and_nonfinite_draws_reported_not_replaced():
    panels = _paired_panels()
    panels[0]["RA"].update(a=None, b=None, ratio=None, converged=False)
    empirical = {"RA_a": [np.nan]*10, "G_a": [.04]*10}
    result = cal.calibrate_inference(panels, empirical, bootstrap=20)
    assert result["summary"]["paired"]["RA"]["v_draws"]["n_nonfinite"] == 20
    assert result["summary"]["paired"]["paired_finite_draws"] == 0
    cells = cal.summarize_cells(panels)
    assert cells[0]["RA"]["failed"] == 1
    assert cells[0]["RA"]["a"]["n_nonfinite"] == 1


def test_gi_source_mapping_requires_source_bootstrap():
    draws = {"RA_a": [1,2], "RA32_a": [3,4], "RA32source_a": [5,6], "G_a": [7,8]}
    np.testing.assert_array_equal(cal._paired_empirical(draws, 3.2).RA_a, [3,4])
    np.testing.assert_array_equal(cal._paired_empirical(draws, 3.2, 2.1).RA_a, [5,6])
    with pytest.raises(ValueError):
        cal._paired_empirical({"RA_a": [1], "G_a": [2]}, 3.2)
    source = [s for s in cal.scenario_grid() if s.name == "matched32_source_sd21_binomial"][0]
    assert source.true_gi_sd == source.analyst_gi_sd == 2.1


def test_batch_merge_checks_completeness_duplicates_and_raw_hashes(tmp_path):
    import json
    scenario = cal.Scenario()
    b0, b1 = tmp_path / "b0", tmp_path / "b1"
    for start, output in ((0, b0), (2, b1)):
        cal.run_calibration(output, replicates=2, replicate_start=start, bootstrap=0,
                            scenarios=[scenario], simulate_only=True, progress=False)
    draws = {"RA_a": [.01, .02, .03], "G_a": [.04, .08, .12]}
    with pytest.raises(ValueError, match="Incomplete replicate range"):
        cal.merge_calibration_batches([b0], tmp_path / "incomplete", empirical_draws=draws, expected_replicates=4)
    with pytest.raises(ValueError, match="Duplicate"):
        cal.merge_calibration_batches([b0, b0], tmp_path / "duplicate", empirical_draws=draws, expected_replicates=4)
    result = cal.merge_calibration_batches([b0, b1], tmp_path / "combined", empirical_draws=draws,
                                           expected_replicates=4, bootstrap=10)
    assert result["manifest"]["panels"] == 20
    bundle = json.loads((tmp_path / "combined" / "calibration.json").read_text(encoding="utf-8"))
    assert bundle["meta"]["replicates_per_cell"] == 4
    assert len(bundle["cells"]) == 5
    assert bundle["inference"][scenario.name]["bootstrap"] == 10
    with (b0 / "calibration_panels.jsonl").open("a", encoding="utf-8") as stream:
        stream.write("\n")
    with pytest.raises(ValueError, match="Changed batch"):
        cal.merge_calibration_batches([b0, b1], tmp_path / "tampered", empirical_draws=draws, expected_replicates=4)


def test_import_does_not_touch_data_or_write_files(tmp_path):
    code = """
import builtins, importlib.util, os, pathlib, sys
import numpy, pandas, scipy.special
source = pathlib.Path(sys.argv[1]).read_text(encoding='utf-8')
def forbidden(*args, **kwargs):
    raise AssertionError('Calibration import attempted I/O')
builtins.open = forbidden
pathlib.Path.open = forbidden
pathlib.Path.mkdir = forbidden
pandas.read_csv = forbidden
namespace = {'__name__': 'import_contract', '__file__': sys.argv[1]}
import types
module = types.ModuleType('import_contract')
sys.modules['import_contract'] = module
module.__dict__.update(namespace)
exec(compile(source, sys.argv[1], 'exec'), module.__dict__)
assert 's09_inversion' not in sys.modules
assert 's10_polish' not in sys.modules
assert 'cc_common' not in sys.modules
assert 'revision_data' not in sys.modules
"""
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    result = subprocess.run([sys.executable, "-c", code, str(CODE / "revision_calibration.py")],
                             cwd=tmp_path, env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert not list(tmp_path.iterdir())
