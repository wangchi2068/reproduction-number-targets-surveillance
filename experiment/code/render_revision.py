"""Render the revision's bilingual figures, tables, and numeric macros from JSON.

No historical results, raw data, or fitted models are read. Importing does no I/O.
Production requires empirical, fixed_effects, thinning, calibration, design,
secondary, and theory JSON in --results-dir. Optional block2/block3 FE summaries
are consumed when present. --allow-incomplete is an explicit PILOT mode: missing
panels/tables are labelled, missing numeric macros are empty (never invented),
and the sibling validation_work/core_inference_pilot/thinning.json may be used.

Outputs under --repo-root: paper_{en,zh}/figures/fig_01_targets.pdf through
fig_04_heterogeneity.pdf, tables/tab_s01_pooled.tex through tab_s08_secondary.tex,
generated/numbers.tex, generated/numbers_sources.json, and revision_manifest.json.
The two numbers.tex files are byte-identical. Each manifest asset names the input
JSON SHA256s and this generator's SHA256. NPZ files are not needed or consumed.
PNG previews are written only when --preview-dir is supplied, outside paper dirs.
Tables require booktabs, array, tabularx and longtable in the owning manuscript.

ThinMCStatus is numeric: 1=adequate, 0=tolerance not met, -1=unavailable. ThinMCMax
is the maximum finite metric MC-SE/sampling-SE ratio; undefined metrics remain
explicit in the table and manifest. Positive CV differences favor the affine
model. Percent reduction macros contain numbers only; manuscripts add percent.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import sys
import tempfile
from typing import Any

VERSION = "render-revision-1.0"
PRIMARY = ("empirical", "fixed_effects", "thinning", "calibration", "design", "secondary", "theory")
STATS = ("RA", "G", "RA32", "RA32source")
CENTERS = ("pooled", "season", "season_linear", "week", "week_leave_state_out")
FIGURES = ("fig_01_targets", "fig_02_models", "fig_03_thinning", "fig_04_heterogeneity")
TABLES = ("pooled", "cv", "thinning", "fe", "centering", "calibration", "design", "secondary")
COLORS = ("#0072B2", "#D55E00", "#009E73", "#CC79A7")
_MISSING = object()


class RenderError(ValueError):
    """Missing, incompatible, or internally inconsistent authoritative results."""


def tr(lang, en, zh):
    return zh if lang == "zh" else en


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def at(obj, *keys, default=None):
    for key in keys:
        if isinstance(obj, dict) and key in obj:
            obj = obj[key]
        elif isinstance(obj, (list, tuple)) and isinstance(key, int) and 0 <= key < len(obj):
            obj = obj[key]
        else:
            return default
    return obj


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def number(value, digits=4, scale=1):
    """Plain decimal text, including small signed values; no scientific TeX."""
    if not finite(value):
        return ""
    value *= scale
    if value == 0:
        return "0"
    places = max(0, digits - 1 - int(math.floor(math.log10(abs(value)))))
    text = f"{value:.{min(places, 12)}f}"
    return text.rstrip("0").rstrip(".") if "." in text else text


def interval(bounds, digits=4, scale=1):
    if not isinstance(bounds, (list, tuple)) or len(bounds) != 2 or not all(finite(x) for x in bounds):
        return ""
    if bounds[0] > bounds[1]:
        raise RenderError(f"Reversed interval: {bounds}")
    return "[" + ", ".join(number(x, digits, scale) for x in bounds) + "]"


def shown(value, digits=4, scale=1):
    return number(value, digits, scale) or "--"


def estimate(value, digits=4, scale=1):
    if not isinstance(value, dict):
        return "--"
    point = shown(value.get("value"), digits, scale)
    ci = interval(value.get("ci95"), digits, scale)
    return point + (" " + ci if ci else " [--]")


def draw_counts(value):
    return (at(value, "n_finite", default=at(value, "valid_draws")),
            at(value, "n_draws", default=at(value, "total_draws")))


def draws_text(value):
    valid, total = draw_counts(value)
    return f"{shown(valid)}/{shown(total)}"


def tex_escape(value):
    """Escape ordinary text once, character by character (never re-escape TeX)."""
    replacements = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$",
                    "#": r"\#", "_": r"\_\allowbreak{}", "{": r"\{", "}": r"\}",
                    "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}
    return "".join(replacements.get(c, c) for c in str(value).replace("\n", " "))


def table_tex(caption, label, headers, rows, notes=(), *, long=False, pilot=False):
    """One table renderer; all caller text is escaped, labels alone are validated."""
    if not re.fullmatch(r"tab:[A-Za-z0-9_-]+", label):
        raise RenderError(f"Invalid table label {label!r}")
    n = len(headers)
    if n < 2 or any(len(row) != n for row in rows):
        raise RenderError(f"Table {label}: inconsistent column counts")
    caption = ("PILOT / " if pilot else "") + caption
    row_tex = lambda row: " & ".join(tex_escape(cell) for cell in row) + r" \\"
    if not rows:
        rows = [["Unavailable"] + ["--"] * (n - 1)]
    # All tables are ordinary floats with wrapping paragraph columns. Long tables
    # use a smaller font so that a 48-row table still fits one page.
    widths = ([.24] + [.70 / (n - 1)] * (n - 1)) if n > 3 else ([.35] + [.60 / (n - 1)] * (n - 1))
    spec = "@{}" + "".join(r">{\raggedright\arraybackslash}p{" + f"{w:.4f}" + r"\linewidth}" for w in widths) + "@{}"
    size = r"\scriptsize" if long else r"\small"
    lines = [r"\begin{table}[htbp]", r"\centering", size,
             r"\caption{" + tex_escape(caption) + r"}\label{" + label + "}",
             r"\setlength{\tabcolsep}{2pt}", r"\begin{tabularx}{\linewidth}{" + spec + "}",
             r"\toprule", row_tex(headers), r"\midrule"]
    lines += [row_tex(row) for row in rows]
    lines += [r"\bottomrule", r"\end{tabularx}"]
    if notes:
        lines += [r"\par\smallskip", r"\begin{minipage}{\linewidth}\footnotesize\raggedright",
                  tex_escape(" ".join(str(note) for note in notes)), r"\end{minipage}"]
    lines += [r"\end{table}", ""]
    return "\n".join(lines)


class Inputs:
    def __init__(self, results_dir, repo_root, allow_incomplete=False):
        self.results_dir, self.repo_root = Path(results_dir), Path(repo_root)
        self.allow_incomplete = allow_incomplete
        self.data, self.sources, self.issues = {}, {}, []
        for name in PRIMARY:
            path = self.results_dir / (name + ".json")
            if not path.is_file() and name == "thinning" and allow_incomplete:
                pilot = Path(__file__).resolve().parents[2].parent / "validation_work/core_inference_pilot/thinning.json"
                if pilot.is_file():
                    path = pilot
                    self.issues.append(f"thinning: explicit pilot source {pilot}")
            self.load(name, path, required=True)
        for block in (2, 3):
            self.load(f"block{block}", self.results_dir / f"block{block}/fixed_effects.json", required=False)
        self.validate()

    def issue(self, message):
        if not self.allow_incomplete:
            raise RenderError(message)
        self.issues.append(message)

    def load(self, name, path, required):
        if not path.is_file():
            if required:
                self.issue(f"Missing authoritative input: {path}")
            return
        raw = path.read_bytes()
        try:
            data = json.loads(raw)
        except (UnicodeError, json.JSONDecodeError) as exc:
            raise RenderError(f"Invalid JSON {path}: {exc}") from exc
        if not isinstance(data, dict):
            raise RenderError(f"Expected JSON object in {path}")
        self.data[name] = data
        self.sources[name] = {"path": str(path.resolve()), "sha256": sha256(raw), "bytes": len(raw)}

    def get(self, name, *keys, default=None):
        return at(self.data.get(name), *keys, default=default)

    def need(self, name, *keys, kind=None):
        value = self.get(name, *keys, default=_MISSING)
        if value is _MISSING or (kind and not isinstance(value, kind)):
            self.issue(f"Missing or incompatible {name}.json key: {'.'.join(map(str, keys))}")
            return None
        return value

    def validate(self):
        # Required sections are checked before writing any published asset.
        required = {"empirical": ("pooled", "validation", "centering"),
                    "fixed_effects": ("metadata", "estimates", "point_fits", "bootstrap_diagnostics", "diagnostics"),
                    "thinning": ("metadata", "estimates", "simultaneous_bands", "monte_carlo", "fit_diagnostics"),
                    "calibration": ("meta", "inference", "cells"),
                    "design": ("meta", "generator", "designs"),
                    "secondary": ("meta", "seasons", "participation"),
                    "theory": ("parameters", "checks", "finite_mean_checks")}
        for name, keys in required.items():
            if name in self.data:
                for key in keys:
                    self.need(name, key, kind=(list if key in ("cells", "checks", "finite_mean_checks", "designs") else dict))
        for name in ("fixed_effects", "thinning", "block2", "block3"):
            if name in self.data and self.get(name, "metadata", "status") != "full_budget":
                self.issue(f"{name}.json is not full_budget: {self.get(name, 'metadata', 'status')}")
        if "calibration" in self.data:
            if self.get("calibration", "meta", "status") != "completed" or not self.get("calibration", "inference"):
                self.issue("calibration.json has no completed inference")
            declared = self.get("calibration", "meta", "scenario_configurations", default=self.get("calibration", "meta", "config", "scenarios", default=[]))
            actual = self.get("calibration", "inference", default={}) or {}
            for scenario in declared:
                if scenario.get("name") not in actual:
                    self.issue(f"Missing calibration inference scenario {scenario.get('name')}")
        if "empirical" in self.data:
            for stat in STATS:
                for key in ("point", "a", "b", "ratio", "bins", "power_fit"):
                    self.need("empirical", "pooled", stat, key)
                for axis in ("state", "season"):
                    self.need("empirical", "validation", stat, axis, "folds", kind=list)
            for stat in ("RA", "G"):
                for center in CENTERS:
                    self.need("empirical", "centering", stat, center, "a", kind=dict)
        if "design" in self.data:
            cells = self.get("design", "designs", default=[])
            expected = {(s, m) for s in (4, 8, 12, 20) for m in (1, 2, 4)}
            observed = {(row.get("seasons"), row.get("count_multiplier")) for row in cells}
            if observed != expected or len(cells) != len(expected):
                self.issue("design.json must contain exactly S=4,8,12,20 x count=1,2,4")


def mc_summary(inputs):
    mc = inputs.get("thinning", "monte_carlo", default={})
    metrics = mc.get("metrics", {})
    ratios = [v.get("mc_to_sampling_se") for v in metrics.values()]
    finite_ratios = [v for v in ratios if finite(v)]
    # s03b stores MC SE / |point estimate| under the key mc_to_sampling_se.
    # Also compute the literal MC SE / sampling SE ratio for transparent reporting.
    literal = []
    for value in metrics.values():
        mc_se = at(value, "bootstrap_mean_mc_se")
        sampling_se = at(value, "estimated_sampling_se")
        if finite(mc_se) and finite(sampling_se) and sampling_se > 0:
            literal.append(mc_se / sampling_se)
    return {"maximum": max(finite_ratios) if finite_ratios else None,
            "max_mc_over_sampling": max(literal) if literal else None,
            "undefined_metrics": sum(not finite(v) for v in ratios),
            "n_metrics": len(metrics), "status": 1 if mc.get("adequate") is True else (0 if mc.get("adequate") is False else -1),
            "text": mc.get("status", "unavailable")}


def make_numbers(inputs):
    values, sources, missing = {}, {}, []

    def add(name, source, keys, *, part=None, scale=1, digits=4, optional=False, transform=None):
        raw = inputs.get(source, *keys)
        value = transform(raw) if transform else raw
        text = interval(value, digits, scale) if part == "interval" else number(value, digits, scale)
        if not text:
            missing.append(name)
            if not optional:
                inputs.issue(f"Numeric macro {name} unavailable at {source}.json:{keys}")
        values[name] = text
        sources[name] = {"input": source, "keys": list(keys), "kind": part or "point", "scale": scale,
                         "digits": digits, "value": value, "formatted": text, "optional": optional}
        if transform:
            sources[name]["transform"] = "finite quantiles [0,2]" if part == "interval" else "derived"

    def pair(prefix, source, *keys, scale=1, digits=4, optional=False):
        add(prefix + "Point", source, (*keys, "value"), scale=scale, digits=digits, optional=optional)
        add(prefix + "Interval", source, (*keys, "ci95"), part="interval", scale=scale, digits=digits, optional=optional)

    for stat in ("RA", "G"):
        for metric in ("a", "b", "ratio"):
            pair(stat + metric, "empirical", "pooled", stat, metric)
        pair(stat + "weekReduction", "empirical", "centering", stat, "week_reduction_fraction", scale=100, digits=3)
        for axis in ("state", "season"):
            add(stat + axis + "CVPoint", "empirical", ("validation", stat, axis, "mean_score_difference"))
    for prefix, key in (("ThinAll", "full/full/c"), ("ThinNoq", "full/noq01/c"),
                        ("ThinRatio", "full/full/c_over_bP"), ("ThinHigh", "highCt50/full/c"),
                        ("ThinBoth", "both50/full/c")):
        pair(prefix, "thinning", "estimates", key)
    add("ThinBPPoint", "thinning", ("estimates", "full/bP", "value"))
    for prefix, key in (("FEb", "state_b"), ("FEa", "a_mean"), ("FEratio", "r_mean"),
                        ("FEpaired", "paired_ratio"), ("FEDratio", "state_D_ratio")):
        pair(prefix, "fixed_effects", "estimates", key)
    for block, word in ((2, "Two"), (3, "Three")):
        add(f"FEpairedBlock{word}Interval", f"block{block}", ("estimates", "paired_ratio", "ci95"), part="interval", optional=True)
        add(f"FEBlock{word}Valid", f"block{block}", ("estimates", "paired_ratio", "n_finite"), optional=True)
    for name, stat, key in (("Nra", "RA", "n"), ("Ngs", "G", "n"), ("Nweeks", "RA", "weeks"), ("Nstates", "RA", "states")):
        add(name, "empirical", ("pooled", stat, key), digits=12)
    add("FEvalid", "fixed_effects", ("estimates", "r_mean", "n_finite"), digits=12)
    for name, source, key in (("FEBootstrap", "fixed_effects", "bootstrap"), ("ThinBootstrap", "thinning", "bootstrap"),
                              ("ThinPointReps", "thinning", "point_reps"), ("ThinInnerReps", "thinning", "inner_reps")):
        add(name, source, ("metadata", key), digits=12)
    mc = mc_summary(inputs)
    values["ThinMCMax"] = number(mc["maximum"])
    values["ThinMCStatus"] = str(mc["status"])
    for name, key in (("ThinMCMax", "maximum"), ("ThinMCStatus", "status")):
        sources[name] = {"input": "thinning", "keys": ["monte_carlo", "metrics" if key == "maximum" else "adequate"],
                         "value": mc[key], "formatted": values[name],
                         "transform": "maximum finite mc_to_sampling_se" if key == "maximum" else "true=1, false=0, missing=-1"}
    if not values["ThinMCMax"]:
        missing.append("ThinMCMax")
        inputs.issue("ThinMCMax unavailable: no finite Monte Carlo/sampling SE ratios")
    for prefix, scenario in (("CalSource", "matched32_source_sd21_binomial"), ("CalLegacy", "matched36_binomial")):
        add(prefix + "Point", "calibration", ("inference", scenario, "RA", "v_implied"))
        add(prefix + "Interval", "calibration", ("inference", scenario, "RA", "v_draws", "q025_50_975_finite"),
            part="interval", transform=lambda x: [x[0], x[2]] if isinstance(x, list) and len(x) == 3 else None)
    lines = ["% Auto-generated by experiment/code/render_revision.py; do not edit.",
             "% Numeric-only macros. Empty definitions denote unavailable PILOT or optional inputs.",
             "% ThinMCStatus: 1=adequate; 0=tolerance not met; -1=unavailable."]
    for name, text in values.items():
        source = sources[name]
        lines.append("% " + name + " <- " + source["input"] + ".json " + json.dumps(source["keys"], ensure_ascii=True))
        lines.append("\\newcommand{\\" + name + "}{" + text + "}")
    return "\n".join(lines) + "\n", sources, missing


def failure_notes(records, lang):
    pieces = []
    for name, value in records:
        valid, total = draw_counts(value)
        if finite(valid) and finite(total):
            pieces.append(f"{name}: {int(valid)}/{int(total)} " + tr(lang, "finite draws", "有限重抽样") +
                          (tr(lang, "; interval conditional on finite draws", "；区间以有限估计为条件") if valid < total else ""))
    return "; ".join(pieces) + ("." if pieces else "")


def pooled_table(i, lang):
    rows, failures = [], []
    for stat in STATS:
        d = i.get("empirical", "pooled", stat, default={})
        for metric in ("a", "b", "ratio"):
            rows.append([stat, metric if metric != "ratio" else "b/a", shown(d.get("n")), estimate(d.get(metric)), draws_text(d.get(metric))])
            failures.append((f"{stat} {metric}", d.get(metric)))
    notes = [tr(lang, "Pooled rows: primary bootstrap (2,000 draws); the same state-week weights are used for every statistic and generation interval. Signed percentile intervals; ratios are defined only for positive denominators. Finite-draw counts are retained for every estimate.",
                "合并行：主自助抽样（2000 次）；各统计量与代际间隔使用同一组州周权重。有符号百分位区间；比值只在分母为正时定义。逐项保留有限重抽样次数。"), failure_notes(failures, lang)]
    return tr(lang, "Pooled affine loss fits and generation-interval sensitivities (bootstrap 95% confidence interval in brackets)", "合并仿射损失拟合与代际间隔敏感性（括号内为 bootstrap 95% 置信区间）"), [tr(lang, "Statistic", "统计量"), tr(lang, "Parameter", "参数"), "n", tr(lang, "Estimate [95% CI]", "估计值 [95% CI]"), tr(lang, "Finite/B", "有限/B")], rows, notes, False


def cv_table(i, lang):
    rows = []
    for stat in STATS:
        for axis in ("state", "season"):
            v = i.get("empirical", "validation", stat, axis, default={})
            rows.append([stat + " / " + tr(lang, axis, "州" if axis == "state" else "季"), tr(lang, "Summary", "汇总"),
                         f"{shown(v.get('valid_folds'))}/{shown(v.get('total_folds'))}", "--", "--",
                         shown(v.get("mean_score_difference")), shown(v.get("affine_better_folds"))])
            if axis == "season":
                for f in v.get("folds", []):
                    rows.append([stat, str(f.get("fold", "--")), shown(f.get("test_n")), shown(f.get("score_affine")),
                                 shown(f.get("score_power")), shown(f.get("power_minus_affine")),
                                 tr(lang, "yes", "是") if f.get("converged_and_positive") else tr(lang, "unavailable", "不可用")])
    notes = [tr(lang, "Summary rows show valid/total folds and the number favoring affine in the last column; season rows show test n and convergence. Score is mean L/mu + log(mu), lower is better. Difference = power minus affine, so positive favors affine. The summary difference is the mean across folds weighted by held-out test n; the unweighted fold mean and the count of folds favouring affine are shown separately. Fold fits refit the training center. Training sets overlap; no iid-fold SE or four-season significance test is reported.",
                "汇总行列出有效/总折数及末列偏向仿射模型的折数；季节行列出测试样本量与收敛状态。评分为平均 L/mu + log(mu)，越低越好。差值=幂律减仿射，正值偏向仿射。汇总差值为以测试样本量加权的折均值；同时列出未加权折均值与仿射胜出折数。各折重新拟合训练中心。训练集重叠，不报告独立折标准误或四季显著性检验。")]
    return tr(lang, "Foldwise model validation", "逐折模型验证"), [tr(lang, "Statistic / split", "统计量/划分"), tr(lang, "Fold", "折"), tr(lang, "Valid/total or n", "有效/总数或n"), tr(lang, "Affine", "仿射"), tr(lang, "Power", "幂律"), tr(lang, "Difference", "差值"), tr(lang, "Wins / valid", "胜出/有效")], rows, notes, True


def thinning_table(i, lang):
    rows = []
    estimates = i.get("thinning", "estimates", default={})
    qs = i.get("thinning", "metadata", "config", "qs", default=[])
    for subset in ("full", "highCt50", "both50"):
        for q in qs:
            base = f"{subset}/q{q:g}"
            for metric in ("b", "c", "delta_a", "c_over_bP"):
                value = estimates.get(base + "/" + metric, {})
                band = i.get("thinning", "simultaneous_bands", subset + "/" + metric, default={})
                simultaneous = dict(zip(band.get("metrics", []), band.get("ci95", []))).get(base + "/" + metric)
                rows.append([subset, "q=" + number(q), metric, estimate(value), interval(simultaneous) or "--", draws_text(value)])
        for fit in ("full", "noq01"):
            for metric in ("c", "c_over_bP", "delta_a_rate"):
                value = estimates.get(f"{subset}/{fit}/{metric}", {})
                rows.append([subset, fit, metric, estimate(value), "--", draws_text(value)])
    meta = i.get("thinning", "metadata", default={})
    diag = i.get("thinning", "fit_diagnostics", default={})
    mc = mc_summary(i)
    notes = [tr(lang, "Per-q intervals are pointwise; the separate bands are simultaneous over q within each subset. full uses all q; noq01 excludes q=0.1. highCt50 requires original current count >=50; both50 requires both original adjacent counts >=50. c is the slope increment per (1/q-1); bP is the independent-Poisson reference (plug-in delta value).",
                "逐q区间为点态区间；单列包络为各子集中跨q的同时区间。full使用全部q；noq01排除q=0.1。highCt50要求原当期计数>=50；both50要求原相邻两期计数均>=50。c为每(1/q-1)的斜率增量；bP为独立泊松参照（代入式 delta 值）。"),
             f"B={shown(meta.get('bootstrap'))}; R(point)={shown(meta.get('point_reps'))}; R(inner, actual)={shown(meta.get('inner_reps'))}; "
             f"precision criterion: MC SE at most {shown(i.get('thinning', 'monte_carlo', 'tolerance'))} times the point estimate, applied to the five monitored response summaries; "
             f"max finite MC SE/|point| over all stored metrics={shown(mc['maximum'])}; "
             f"max finite MC SE/sampling SE over all stored metrics={shown(mc['max_mc_over_sampling'])}; "
             f"undefined MC ratios={mc['undefined_metrics']}/{mc['n_metrics']}.",
             f"Point fits {shown(diag.get('point_converged'))}/{shown(diag.get('point_attempted'))}; "
             f"bootstrap fits {shown(diag.get('bootstrap_converged'))}/{shown(diag.get('bootstrap_attempted'))}; "
             f"failed bootstrap averages={shown(diag.get('failed_bootstrap_averages'))}.",
             str(i.get("thinning", "failure_policy", default="Unavailable"))]
    for subset in ("full", "highCt50", "both50"):
        notes.append(f"{subset}: n={shown(i.get('thinning', 'metadata', 'config', 'subsets', subset))}; bP={estimate(estimates.get(subset + '/bP'))}.")
    return tr(lang, "Shared-count binomial thinning of the R_A sample (n=1748): estimates, simultaneous bands and simulation diagnostics", "共享计数稀释估计、同时区间与模拟诊断"), [tr(lang, "Sample", "样本"), tr(lang, "q / fit", "q/拟合"), tr(lang, "Metric", "指标"), tr(lang, "Estimate [95% CI]", "估计值 [95% CI]"), tr(lang, "Simultaneous 95%", "同时95%区间"), tr(lang, "Finite/B", "有限/B")], rows, notes, True


def fe_table(i, lang):
    rows, notes = [], []
    metric_names = (("pool_a", "Pooled a"), ("pool_b", "Pooled b"), ("pool_ratio", "Pooled b/a"),
                    ("state_b", "State FE b"), ("a_mean", "State FE mean a"), ("r_mean", "State FE b/mean a"),
                    ("a_median", "State FE median a"), ("r_median", "State FE b/median a"),
                    ("a_row", "State FE row-weighted a"), ("r_row", "State FE row-weighted ratio"),
                    ("a_case", "State FE count-weighted a"), ("r_case", "State FE count-weighted ratio"),
                    ("state_D_b", "D_t scale b"), ("state_D_mean", "D_t scale mean a"),
                    ("state_D_ratio", "D_t scale ratio"), ("slope_ratio", "Slope ratio (pooled b / state FE b)"),
                    ("intercept_ratio", "Intercept ratio (state FE mean a / pooled a)"),
                    ("paired_ratio", "Cross-count ratio (pooled b/a / state FE b/mean a)"))
    for source in ("fixed_effects", "block2", "block3"):
        if source not in i.data:
            if source != "fixed_effects":
                notes.append(source + ": " + tr(lang, "unavailable (optional)", "不可用（可选）"))
            continue
        block = i.get(source, "metadata", "config", "block_length")
        for key, label in metric_names:
            value = i.get(source, "estimates", key, default={})
            rows.append(["L=" + shown(block), label, estimate(value), draws_text(value), shown(value.get("nonpositive_fraction"))])
        for key, d in i.get(source, "bootstrap_diagnostics", default={}).items():
            notes.append(f"L={shown(block)} {key}: {shown(d.get('converged'))} converged, {shown(d.get('failed'))} failed; " + json.dumps(d.get("failure_reasons", {}), ensure_ascii=False) + ".")
    for key in ("week_only", "state_and_week"):
        d = i.get("fixed_effects", "diagnostics", key, default={})
        rows.append([tr(lang, "Point diagnostic", "点估计诊断"), key, "b=" + shown(d.get("b")) if d.get("converged") else tr(lang, "Unavailable", "不可用"),
                     tr(lang, "No CI", "无区间"), str(d.get("failure") or d.get("method", "unavailable"))])
    notes.insert(0, tr(lang, "Sample: the R_A specification, n=1748 state-weeks in 51 jurisdictions; the G sample is not used here. Intervals retain signed draws and are conditional on finite estimates when finite/B < 1. Last column gives the nonpositive fraction over all requested draws. FE targets the sampled state occurrences; pooled and FE quantities are computed on the same draw. Cross-count ratio = slope ratio x intercept ratio. L denotes the within-season week-block length. Week-only and state-and-week fits are point diagnostics; failed fits have no estimate or interval. L=1 uses the primary bootstrap (2,000 draws, the weights of Table S1). L=2 and L=3 use 1,000 draws each, and their ratios use the week-block bootstrap. Pooled row: primary bootstrap (2,000 draws).",
                         "样本为 $R_A$ 规格的 1748 个州周、51 个辖区；本表不使用 $G$ 样本。区间保留有符号重抽样；有限/B<1时，区间以有限估计为条件。末列为全部请求重抽样中的非正比例。FE以抽中的州出现次数为目标；合并与FE量在同一次重抽样上计算。跨计数比=斜率比×截距比。L为季内周块长度。仅周和州加周拟合为点估计诊断，失败时无估计或区间。L=1 使用主自助抽样（2000 次，与表 S1 权重相同）；L=2、3 各 1000 次，其合并行同样取自主自助抽样（2000 次），比值使用周块自助抽样。"))
    return tr(lang, "Cross-count ratios (pooled b/a over state fixed-effects b/a), slope ratios and intercept ratios for the R_A sample (n=1748), with week-block sensitivity and convergence", "R_A 样本（n=1748）的跨计数比（合并 b/a 与州固定效应 b/a 之比）、斜率比与截距比，含周块敏感性与收敛"), [tr(lang, "Block", "周块"), tr(lang, "Parameter", "参数"), tr(lang, "Estimate [95% CI]", "估计值 [95% CI]"), tr(lang, "Finite/B", "有限/B"), tr(lang, "Nonpositive / diagnostic", "非正比例/诊断")], rows, notes, True


def centering_table(i, lang):
    rows, checks = [], []
    for stat in ("RA", "G"):
        normalizer = i.get("empirical", "pooled", stat, "point", "xbar")
        for center in CENTERS:
            d = i.get("empirical", "centering", stat, center, default={})
            rows.append([stat + " / " + center, shown(normalizer), estimate(d.get("a")), estimate(d.get("b")), draws_text(d.get("a"))])
            checks.append((stat + "/" + center, d.get("a")))
        for key in ("week_change", "week_reduction_fraction"):
            d = i.get("empirical", "centering", stat, key, default={})
            rows.append([stat + " / " + key, shown(normalizer), estimate(d, scale=100 if key.endswith("fraction") else 1), "--", draws_text(d)])
    notes = [tr(lang, "All five centers use the same pooled mean as the relative-loss denominator within each statistic; the denominator is refitted with weights within each bootstrap draw. week_leave_state_out removes the focal state from its weekly center. Reduction is 100*(a_pooled-a_week)/a_pooled, in percent; changes are descriptive and not additive variance components. Pooled-centre rows: primary bootstrap (2,000 draws, as in Table S1); the other centres and the week change use this table's own 1000-draw bootstrap.",
                "同一统计量的五种中心均使用相同合并均值作为相对损失分母；每次重抽样用权重重新拟合该分母。week_leave_state_out的周中心排除本州。减少比例为100*(a_pooled-a_week)/a_pooled，单位为百分比；这些是描述性变化，并非可加总方差分量。合并中心行取自主自助抽样（2000 次，同表 S1）；其余中心与周变化使用本表自身的 1000 次自助抽样。"), failure_notes(checks, lang)]
    return tr(lang, "Centering stages with a common loss normalizer", "采用共同损失分母的中心化阶段"), [tr(lang, "Statistic / center", "统计量/中心"), tr(lang, "Pooled mean", "合并均值"), "a / " + tr(lang, "change [95% CI]", "变化 [95% CI]"), "b [95% CI]", tr(lang, "Finite/B", "有限/B")], rows, notes, True


def calibration_table(i, lang):
    rows = []
    inference = i.get("calibration", "inference", default={}) or {}
    cells = i.get("calibration", "cells", default=[]) or []
    for scenario, rec in sorted(inference.items()):
        for stat in ("RA", "G"):
            d = rec.get(stat, {})
            draws = d.get("v_draws", {})
            q = draws.get("q025_50_975_finite", [])
            ci = interval([q[0], q[2]]) if len(q) == 3 else ""
            fit_cells = [c for c in cells if c.get("scenario") == scenario]
            failed = sum(at(c, stat, "failed", default=0) for c in fit_cells)
            attempted = sum(c.get("panels", 0) for c in fit_cells)
            rows.append([scenario + " / " + stat, shown(d.get("v_implied")) + " " + (ci or "[--]"),
                         f"{shown(draws.get('n_finite'))}/{shown(draws.get('n_total'))}",
                         f"{shown(draws.get('n_nonpositive'))}/{shown(draws.get('n_nonfinite'))}",
                         f"{shown(d.get('n_outside_grid'))}/{shown(d.get('n_nonpositive_slopes'))}",
                         f"{failed}/{attempted}"])
    notes = [tr(lang, "v is a model-conditional weekly marginal variance before the lower bound, not an identified biological variance or upper bound. Intervals are the 2.5th and 97.5th percentiles of all finite signed draws; no positive-only filtering or clamping is applied. Columns report nonpositive/nonfinite inversion draws, out-of-grid/nonpositive-slope draws, and failed/attempted simulation fits.",
                "v为给定模型下、施加下界之前的周边际方差，不是已识别的生物学方差或上界。区间取全部有限有符号重抽样的2.5和97.5百分位；不作仅保留正值的筛选或截断。各列分别列出非正/非有限反演次数、超网格/非正斜率次数，以及失败/尝试的模拟拟合次数。")]
    for scenario, rec in sorted(inference.items()):
        notes.append(scenario + ": grid=" + str(rec.get("v_grid", [])) + "; " + "; ".join(
            stat + " point outside=" + str(at(rec, stat, "point_outside_grid")) +
            ", nonpositive=" + str(at(rec, stat, "point_nonpositive")) +
            ", max |line residual|/MC SE=" + shown(at(rec, stat, "max_abs_residual_mc_se")) for stat in ("RA", "G")) + ".")
    if not inference:
        rows.append([tr(lang, "Calibration unavailable", "校准不可用")] + ["--"] * 5)
    return tr(lang, "Conditional calibration across all declared scenarios", "全部声明情景的条件校准"), [tr(lang, "Scenario / statistic", "情景/统计量"), "v [95% CI]", tr(lang, "Finite/B", "有限/B"), tr(lang, "Nonpositive / nonfinite", "非正/非有限"), tr(lang, "Outside / slope <=0", "超网格/斜率<=0"), tr(lang, "Failed fits / all", "失败拟合/全部")], rows, notes, True


def design_table(i, lang):
    rows = []
    for d in sorted(i.get("design", "designs", default=[]), key=lambda d: (d["seasons"], d["count_multiplier"])):
        rows.append([shown(d.get("seasons")), shown(d.get("count_multiplier")),
                     shown(d.get("size")) + " " + (interval(d.get("size_ci95")) or "[--]"),
                     shown(d.get("power")) + " " + (interval(d.get("power_ci95")) or "[--]"),
                     "/".join(shown(d.get(k)) for k in ("calibration_valid", "null_valid", "alternative_valid")),
                     tr(lang, "yes", "是") if d.get("size_compatible_with_nominal") else tr(lang, "no", "否")])
    notes = [tr(lang, "Conditional discrimination simulations recycle the observed season/count design. Calibration and evaluation samples are independent. Valid counts are calibration/null/alternative. The final column states whether the size interval contains nominal alpha; all cells are reported, including failed size checks. These simulations do not imply a universal season threshold or a universal 80% power claim.",
                "条件区分模拟复用观测季节/计数设计，校准与评估样本独立。有效次数依次为校准/零假设/备择假设。末列表示检验大小区间是否包含名义alpha；保留全部单元，包括大小检验不通过的单元。这些模拟不意味着普适季节门槛或普适80%功效。"),
             f"alpha={shown(i.get('design', 'generator', 'nominal_alpha'))}; "
             f"replicates/stage={shown(i.get('design', 'meta', 'replicates_per_stage'))}."]
    return tr(lang, "Model-conditional design size and power", "给定模型下的设计检验大小与功效"), [tr(lang, "Seasons", "季数"), tr(lang, "Count multiplier", "计数倍数"), tr(lang, "Size [95% CI]", "大小 [95% CI]"), tr(lang, "Power [95% CI]", "功效 [95% CI]"), tr(lang, "Valid C/N/A", "有效C/N/A"), tr(lang, "Size compatible", "大小相容")], rows, notes, False


def secondary_table(i, lang):
    rows, failures = [], []
    for stat in ("RA", "G"):
        for season, d in i.get("secondary", "seasons", stat, default={}).items():
            rows.append([stat + " / " + season, shown(d.get("n")), estimate(d.get("a")),
                         shown(at(d, "point", "b")), shown(at(d, "point", "ratio")), draws_text(d.get("a"))])
            failures.append((stat + "/" + season, d.get("a")))
        d = i.get("secondary", "participation", stat, default={})
        rows.append([stat + " / " + tr(lang, "reporting adjustment", "报告调整"), shown(d.get("n")),
                     estimate(d.get("delta_original_minus_adjusted")), shown(d.get("original_a")),
                     shown(d.get("adjusted_a")), draws_text(d.get("delta_original_minus_adjusted"))])
    notes = [tr(lang, "Season rows report a with its interval and b, b/a as point estimates only. Reporting-adjustment rows use the a column for original-minus-adjusted a [95% CI], b for original a, and b/a for adjusted a; the original normalizer is fixed. Season b and b/a come from separate per-season fits without a stored bootstrap; season b/a with intervals from a common-slope model are given in the season table. Summer rows are reported in a separate table.",
                "季节行报告a及其区间，b与b/a仅为点估计。报告调整行的a列为原值减调整后a [95% CI]，b列为原a，b/a列为调整后a；保持原分母不变。季节 b 与 b/a 来自不存储自助抽样的分季拟合；带区间的季节 b/a（共同斜率模型）见季节表。夏季行另列一表。"), failure_notes(failures, lang)]
    return tr(lang, "Season-specific fits and reporting participation sensitivity", "分季拟合与报告参与率敏感性"), [tr(lang, "Statistic / period", "统计量/时期"), "n", "a / " + tr(lang, "change [95% CI]", "变化 [95% CI]"), "b / " + tr(lang, "original a", "原a"), "b/a / " + tr(lang, "adjusted a", "调整a"), tr(lang, "Finite/B", "有限/B")], rows, notes, True


def summer_table(i, lang):
    rows = []
    for year, d in i.get("secondary", "summer", default={}).items():
        if isinstance(d, dict):
            rows.append([tr(lang, "Summer ", "夏季 ") + year, shown(d.get("n")), shown(d.get("median_count_scale")),
                         shown(d.get("coverage_median")), shown(d.get("iqr_logG"))])
    notes = [tr(lang, "Descriptive summaries of May-October state-weeks with D_t >= 5: rows, median count scale D_t, median reporting participation (%) and IQR of log G. No model is fitted, so no interval is given. Reporting was voluntary in the 2024 summer only; that contrast is descriptive, not a causal or independently replicated intervention.",
                "5-10 月且 D_t>=5 的州周的描述性汇总：行数、计数尺度 D_t 中位数、报告参与率中位数(%)与 log G 的四分位距。未拟合模型，故不给区间。仅 2024 年夏季为自愿报告；该对比仅作描述，不作因果或独立重复干预解释。")]
    return tr(lang, "Off-season summer descriptors", "非流行季夏季描述"), [tr(lang, "Period", "时期"), "n", tr(lang, "Median count scale", "计数尺度中位数"), tr(lang, "Median participation (%)", "参与率中位数(%)"), "IQR(log G)"], rows, notes, False


TABLE_BUILDERS = (pooled_table, cv_table, thinning_table, fe_table, centering_table, calibration_table, design_table, secondary_table)
TABLE_INPUTS = (("empirical",), ("empirical",), ("thinning",), ("fixed_effects", "block2", "block3"),
                ("empirical",), ("calibration",), ("design",), ("secondary",))


def font_config(lang):
    from matplotlib import font_manager
    installed = {f.name for f in font_manager.fontManager.ttflist}
    latin = next((name for name in ("Times New Roman", "Times", "Liberation Serif", "DejaVu Serif") if name in installed), "DejaVu Serif")
    cjk = next((name for name in ("SimSun", "Noto Serif CJK SC", "Noto Serif SC", "Source Han Serif SC") if name in installed), None)
    if lang == "zh" and cjk is None:
        raise RenderError("Chinese PDF rendering requires SimSun or Noto Serif CJK SC (discovered through matplotlib.font_manager)")
    families = [latin] + ([cjk] if lang == "zh" else [])
    return {"font.family": families, "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 9,
            "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
            "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": .7,
            "lines.linewidth": 1.3, "lines.markersize": 3.5, "pdf.fonttype": 42, "ps.fonttype": 42,
            "mathtext.fontset": "stix", "text.usetex": False, "axes.unicode_minus": False,
            "savefig.facecolor": "white"}, {"latin": latin, "cjk": cjk}


def unavailable(ax, lang, reason):
    ax.text(.5, .5, tr(lang, "Unavailable", "不可用") + "\n" + reason, transform=ax.transAxes,
            ha="center", va="center", wrap=True, color="#555555")
    ax.set_axis_off()


def panel_title(ax, letter, title):
    ax.set_title(letter + "  " + title, loc="left", pad=8)


def points_ci(ax, positions, records, *, horizontal=False, color=COLORS[0], marker="o", label=None):
    # Plot endpoints directly: percentile CIs need not contain the point estimate.
    used_label = False
    for pos, record in zip(positions, records):
        point, bounds = at(record, "value"), at(record, "ci95")
        if isinstance(bounds, (list, tuple)) and len(bounds) == 2 and all(finite(x) for x in bounds):
            if horizontal:
                ax.plot(bounds, [pos, pos], color=color, linewidth=1)
                ax.plot(bounds, [pos, pos], "|", color=color)
            else:
                ax.plot([pos, pos], bounds, color=color, linewidth=1)
                ax.plot([pos, pos], bounds, "_", color=color)
        if finite(point):
            ax.plot(point if horizontal else pos, pos if horizontal else point, marker=marker, color=color,
                    linestyle="none", label=label if not used_label else None)
            used_label = True


def fig_targets(i, lang, plt, np):
    fig, axes = plt.subplots(2, 2, figsize=(7.15, 5.35), constrained_layout=True)
    checks = i.get("theory", "checks", default=[])
    finite_checks = i.get("theory", "finite_mean_checks", default=[])
    titles = [tr(lang, "Inferential targets", "推断目标"), tr(lang, "Reporting variation", "报告波动"),
              tr(lang, "Exact and Monte Carlo", "解析值与蒙特卡洛"), tr(lang, "Realised mean", "有限期实现均值")]
    for ax, letter, title in zip(axes.flat, "ABCD", titles): panel_title(ax, letter, title)
    if not checks or not finite_checks:
        for ax in axes.flat: unavailable(ax, lang, "theory.json")
        return fig
    tau0 = min(row["tau"] for row in checks)
    current = sorted((r for r in checks if r["target"] == "current" and r["tau"] == tau0), key=lambda r: r["m"])
    process = {r["m"]: r for r in checks if r["target"] == "process" and r["tau"] == tau0}
    T = finite_checks[0]["T"]
    ax = axes[0, 0]
    x = [r["m"] for r in current]
    ax.loglog(x, [r["analytic"] for r in current], color=COLORS[0], label=tr(lang, "Current value", "当期值"))
    ax.loglog(x, [r["analytic"] / T for r in current], "--", color=COLORS[1], label=tr(lang, "Realised mean", "有限期均值") + f" ($T={T}$)")
    ax.loglog(x, [process[r["m"]]["analytic"] / T for r in current], "-.", color=COLORS[2], label=tr(lang, "Process mean", "过程均值") + f" ($T={T}$)")
    ax.set_xlabel(tr(lang, "Expected count $m$", "期望计数 m")); ax.set_ylabel("MSE"); ax.legend(frameon=False)
    ax = axes[0, 1]
    for j, tau in enumerate(sorted({r["tau"] for r in checks})):
        rows = sorted((r for r in checks if r["target"] == "current" and r["tau"] == tau), key=lambda r: r["m"])
        ax.loglog([r["m"] for r in rows], [r["analytic"] for r in rows], linestyle=("-", "--", "-.")[j % 3],
                  color=COLORS[j % len(COLORS)], label=r"$\tau=$" + number(tau))
    ax.set_xlabel(tr(lang, "Expected count $m$", "期望计数 m")); ax.set_ylabel(tr(lang, "Current-target MSE", "当期目标 MSE")); ax.legend(frameon=False)
    ax = axes[1, 0]
    values = []
    for j, target in enumerate(("current", "process")):
        rows = [r for r in checks if r["target"] == target]
        exact = np.array([r["analytic"] for r in rows]); simulated = np.array([r["monte_carlo"] for r in rows])
        se = np.array([r["mc_se"] for r in rows])
        ax.errorbar(exact, simulated, yerr=1.96 * se, fmt=("o", "s")[j], ms=2.5, alpha=.75, color=COLORS[j],
                    label=tr(lang, "Current value" if target == "current" else "Process mean", "当期" if target == "current" else "过程"), elinewidth=.6)
        values.extend(exact); values.extend(simulated)
    ax.plot([min(values), max(values)], [min(values), max(values)], ":", color="#555555")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel(tr(lang, "Exact MSE", "解析 MSE")); ax.set_ylabel(tr(lang, "MC MSE (±1.96 MC SE)", "MC MSE (±1.96 MC标准误)")); ax.legend(frameon=False)
    ax = axes[1, 1]
    for j, tau in enumerate(sorted({r["tau"] for r in finite_checks})):
        rows = sorted((r for r in finite_checks if r["tau"] == tau), key=lambda r: r["m"])
        ax.loglog([r["m"] for r in rows], [r["analytic"] for r in rows], color=COLORS[j], linestyle=("-", "--")[j % 2], label=r"$\tau=$" + number(tau))
        ax.errorbar([r["m"] for r in rows], [r["monte_carlo"] for r in rows], yerr=[1.96 * r["mc_se"] for r in rows], fmt="o", color=COLORS[j], ms=3)
    ax.set_xlabel(tr(lang, "Expected count $m$", "期望计数 m")); ax.set_ylabel(f"MSE ($T={T}$)"); ax.legend(frameon=False)
    return fig


def fig_models(i, lang, plt, np):
    fig, axes = plt.subplots(2, 2, figsize=(7.15, 5.5), constrained_layout=True)
    for j, stat in enumerate(("RA", "G")):
        ax = axes[0, j]; panel_title(ax, "AB"[j], {"RA": "$R_A$", "G": "$G$"}[stat] + " / " + tr(lang, "observed support", "观测支持域"))
        d = i.get("empirical", "pooled", stat, default={}); bins = d.get("bins", [])
        if not bins:
            unavailable(ax, lang, "empirical.json pooled." + stat); continue
        x = np.array([b["D"] for b in bins]); y = np.array([b["L"] for b in bins])
        ax.scatter(x, y, s=np.sqrt([b["n"] for b in bins]) * 2, color="#333333", label=tr(lang, "Binned loss", "分箱损失"), zorder=3)
        grid = np.geomspace(x.min(), x.max(), 200)
        fit, power = d.get("point", {}), d.get("power_fit", {})
        if fit.get("converged") and finite(fit.get("a")) and finite(fit.get("b")):
            ax.plot(grid, fit["a"] + fit["b"] / grid, color=COLORS[0], label=tr(lang, "Affine", "仿射"))
        if power.get("converged") and finite(power.get("coefficient")) and finite(power.get("exponent")):
            ax.plot(grid, power["coefficient"] * grid ** (-power["exponent"]), "--", color=COLORS[1], label=tr(lang, "Power", "幂律"))
        ax.set_xscale("log"); ax.set_xlabel(tr(lang, "Count scale $D$ (bin means)", "计数尺度 D（箱均值）")); ax.set_ylabel(tr(lang, "Relative squared loss", "相对平方损失")); ax.legend(frameon=False)
    for j, axis in enumerate(("season", "state")):
        ax = axes[1, j]; panel_title(ax, "CD"[j], tr(lang, "Leave-one-" + axis + " validation", "留一" + ("季" if axis == "season" else "州") + "验证"))
        groups = sorted({str(f["fold"]) for stat in ("RA", "G") for f in i.get("empirical", "validation", stat, axis, "folds", default=[])})
        for k, stat in enumerate(("RA", "G")):
            by_fold = {str(f["fold"]): f for f in i.get("empirical", "validation", stat, axis, "folds", default=[])}
            pos, diffs = [], []
            for n, group in enumerate(groups):
                f = by_fold.get(group, {})
                if f.get("converged_and_positive") and finite(f.get("power_minus_affine")):
                    pos.append(n + (k - .5) * .17); diffs.append(f["power_minus_affine"])
            ax.plot(pos, diffs, ("o", "s")[k], color=COLORS[k], label={"RA": "$R_A$", "G": "$G$"}[stat], ms=3)
        ax.axhline(0, color="#555555", linewidth=.8)
        ax.set_ylabel(tr(lang, "Power − affine score", "幂律 − 仿射评分")); ax.set_xlabel(tr(lang, "Held-out " + axis + (" (alphabetical; positive favours affine)" if axis == "state" else " (positive favours affine)"), "留出" + ("季" if axis == "season" else "州") + "（正值偏向仿射）"))
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(groups if axis == "season" else [], fontsize=7)
        ax.legend(frameon=False)
    return fig


def fig_thinning(i, lang, plt, np):
    fig, axes = plt.subplots(2, 2, figsize=(7.15, 5.4), constrained_layout=True)
    qs = i.get("thinning", "metadata", "config", "qs", default=[])
    titles = [tr(lang, "Slope by thinning probability", "逐稀释概率的斜率"), tr(lang, "Intercept change: simultaneous bands", "截距变化：同时区间"),
              tr(lang, "Slope-increment sensitivity", "斜率增量敏感性"), tr(lang, "Monte Carlo precision", "蒙特卡洛精度")]
    for ax, letter, title in zip(axes.flat, "ABCD", titles): panel_title(ax, letter, title)
    if not qs:
        for ax in axes.flat: unavailable(ax, lang, "thinning.json")
        return fig
    order = sorted(range(len(qs)), key=lambda k: 1 / qs[k] - 1)
    qvals = [qs[k] for k in order]; x = [1 / q - 1 for q in qvals]
    ax = axes[0, 0]
    records = [i.get("thinning", "estimates", f"full/q{q:g}/b", default={}) for q in qvals]
    points_ci(ax, x, records, label=tr(lang, "b(q), pointwise 95%", "b(q)，点态95%区间"))
    bp = i.get("thinning", "estimates", "full/bP", "value"); base = i.get("thinning", "estimates", "full/unthinned_b", "value")
    if finite(bp) and finite(base): ax.plot(x, [base + bp * u for u in x], "--", color=COLORS[1], label=tr(lang, "Independent-Poisson reference", "独立泊松参照"))
    ax.set_xlabel("1/q − 1"); ax.set_ylabel("b(q)"); ax.legend(frameon=False)
    ax = axes[0, 1]
    band = i.get("thinning", "simultaneous_bands", "full/delta_a", default={})
    bands = dict(zip(band.get("metrics", []), band.get("ci95", [])))
    records = [{"value": i.get("thinning", "estimates", f"full/q{q:g}/delta_a", "value"), "ci95": bands.get(f"full/q{q:g}/delta_a")} for q in qvals]
    points_ci(ax, x, records, color=COLORS[1]); ax.axhline(0, color="#555555", linewidth=.8)
    ax.set_xlabel("1/q − 1"); ax.set_ylabel("delta a(q)")
    ax = axes[1, 0]
    keys = ("full/full/c", "full/noq01/c", "highCt50/full/c", "both50/full/c")
    labels = [tr(lang, "All q", "全部q"), tr(lang, "Exclude q=0.1", "排除q=0.1"), tr(lang, "Current count >=50", "当期计数>=50"), tr(lang, "Both counts >=50", "两期计数>=50")]
    records = [i.get("thinning", "estimates", key, default={}) for key in keys]
    points_ci(ax, range(len(keys)), records, horizontal=True, color=COLORS[2]); ax.set_yticks(range(len(keys)), labels); ax.invert_yaxis()
    if finite(bp): ax.axvline(bp, color=COLORS[1], linestyle="--", label="bP")
    ax.set_xlabel(tr(lang, "c [pointwise 95% CI]", "c [点态95%区间]")); ax.legend(frameon=False)
    ax = axes[1, 1]
    metrics = i.get("thinning", "monte_carlo", "metrics", default={})
    for j, metric in enumerate(("c", "delta_a")):
        vals = [at(metrics, f"full/q{q:g}/{metric}", "mc_to_sampling_se") for q in qvals]
        ax.plot([u for u, v in zip(x, vals) if finite(v)], [v for v in vals if finite(v)], ("o-", "s--")[j], color=COLORS[j], label=metric)
    tolerance = i.get("thinning", "monte_carlo", "tolerance")
    if finite(tolerance): ax.axhline(tolerance, color="#555555", linestyle=":", label=tr(lang, "0.1 reference", "0.1参考线"))
    ax.set_xlabel("1/q − 1"); ax.set_ylabel("MC SE / point estimate"); ax.legend(frameon=False)
    fig.suptitle(tr(lang, "Thinning", "稀释") + f"  B={shown(i.get('thinning', 'metadata', 'bootstrap'))}, "
                 f"R={shown(i.get('thinning', 'metadata', 'inner_reps'))}", fontsize=9)
    return fig


def fig_heterogeneity(i, lang, plt, np):
    fig, axes = plt.subplots(2, 2, figsize=(7.15, 5.6), constrained_layout=True)
    titles = [tr(lang, "Pooled and within-state ratios", "合并与州内比值"), tr(lang, "Cross-count ratios", "跨计数比"),
              tr(lang, "$R_A$ centring stages", "RA中心化阶段"), tr(lang, "$G$ centring stages", "G中心化阶段")]
    for ax, letter, title in zip(axes.flat, "ABCD", titles): panel_title(ax, letter, title)
    ax = axes[0, 0]
    keys = ("pool_ratio", "r_mean", "state_D_ratio")
    records = [i.get("fixed_effects", "estimates", key, default={}) for key in keys]
    points_ci(ax, range(3), records, horizontal=True)
    ax.set_yticks(range(3), [tr(lang, "Pooled", "合并"), tr(lang, "State FE", "州FE"), tr(lang, "State FE, $D_t$ scale", "州FE，D_t尺度")]); ax.invert_yaxis()
    ax.set_xlabel(tr(lang, "$b/a$ [95% CI]", "b/a [95%区间]"))
    ax = axes[0, 1]
    sources = [s for s in ("fixed_effects", "block2", "block3") if s in i.data]
    records = [i.get(s, "estimates", "paired_ratio", default={}) for s in sources]
    points_ci(ax, range(len(sources)), records, horizontal=True, color=COLORS[2])
    ax.set_yticks(range(len(sources)), ["$L=" + shown(i.get(s, "metadata", "config", "block_length")) + "$; " + draws_text(r) for s, r in zip(sources, records)]); ax.invert_yaxis()
    ax.axvline(1, color="#555555", linestyle=":"); ax.set_xlabel(tr(lang, "Cross-count ratio (pooled / state FE)", "跨计数比（合并/州FE）"))
    labels = [tr(lang, "Pooled", "合并"), tr(lang, "Season", "季节"), tr(lang, "Season trend", "季内趋势"), tr(lang, "Week", "周"), tr(lang, "Week, leave state out", "周，排除本州")]
    for j, stat in enumerate(("RA", "G")):
        ax = axes[1, j]
        records = [i.get("empirical", "centering", stat, center, "a", default={}) for center in CENTERS]
        points_ci(ax, range(len(CENTERS)), records, horizontal=True, color=COLORS[j])
        ax.set_yticks(range(len(CENTERS)), labels); ax.invert_yaxis(); ax.set_xlabel(tr(lang, "Intercept $a$ [95% CI], common normaliser", "a [95%区间]，共同分母"))
        ax.axvline(0, color="#555555", linewidth=.6)
    return fig


FIGURE_BUILDERS = (fig_targets, fig_models, fig_thinning, fig_heterogeneity)
FIGURE_INPUTS = (("theory",), ("empirical",), ("thinning",), ("empirical", "fixed_effects", "block2", "block3"))


def render_revision(results_dir=None, repo_root=None, allow_incomplete=False, preview_dir=None):
    """Render once from a frozen in-memory JSON snapshot; return per-language manifests."""
    repo_root = Path(repo_root or Path(__file__).resolve().parents[2]).resolve()
    results_dir = Path(results_dir or repo_root / "experiment/analysis/revised").resolve()
    preview_dir = Path(preview_dir).resolve() if preview_dir else None
    if preview_dir and any(preview_dir == repo_root / f"paper_{lang}" or repo_root / f"paper_{lang}" in preview_dir.parents for lang in ("en", "zh")):
        raise RenderError("--preview-dir must be outside the published paper_en/paper_zh directories")
    inputs = Inputs(results_dir, repo_root, allow_incomplete)
    numbers, macro_sources, missing_macros = make_numbers(inputs)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    fonts = {lang: font_config(lang) for lang in ("en", "zh")}
    generator = {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__).read_bytes()), "version": VERSION}
    timestamp = datetime.now(timezone.utc).isoformat()
    manifests = {}
    # Render into a staging directory first. Schema/font/render errors leave papers intact.
    with tempfile.TemporaryDirectory(prefix="render_revision_") as temp:
        stage = Path(temp)
        for lang in ("en", "zh"):
            root = stage / f"paper_{lang}"
            assets = []
            def record(relative, dependencies, kind, **extra):
                path = root / relative
                assets.append({"path": relative, "kind": kind, "sha256": sha256(path.read_bytes()),
                               "inputs": {name: inputs.sources[name] for name in dependencies if name in inputs.sources},
                               "missing_inputs": [name for name in dependencies if name not in inputs.sources],
                               "generator": generator, **extra})
            def write(relative, text):
                path = root / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(text, encoding="utf-8", newline="\n")
            write("generated/numbers.tex", numbers)
            write("generated/numbers_sources.json", json.dumps(macro_sources, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
            record("generated/numbers.tex", tuple(inputs.sources), "macros", macro_names=list(macro_sources))
            record("generated/numbers_sources.json", tuple(inputs.sources), "macro_sources")
            for index, (name, builder, dependencies) in enumerate(zip(TABLES, TABLE_BUILDERS, TABLE_INPUTS), 1):
                caption, headers, rows, notes, long = builder(inputs, lang)
                relative = f"tables/tab_s{index:02d}_{name}.tex"
                write(relative, table_tex(caption, "tab:" + name, headers, rows, notes, long=long, pilot=allow_incomplete))
                record(relative, dependencies, "table", label="tab:" + name, caption=caption, rows=len(rows), environment="longtable" if long else "tabularx")
            with matplotlib.rc_context(fonts[lang][0]):
                for name, builder, dependencies in zip(FIGURES, FIGURE_BUILDERS, FIGURE_INPUTS):
                    fig = builder(inputs, lang, plt, np)
                    if allow_incomplete:
                        fig.text(.5, .995, tr(lang, "PILOT — not for publication", "试运行 — 不用于发表"), ha="center", va="top", color="#9A3B00", fontsize=8)
                    path = root / "figures" / (name + ".pdf")
                    path.parent.mkdir(parents=True, exist_ok=True)
                    fig.savefig(path, format="pdf", bbox_inches="tight", metadata={"Creator": VERSION, "CreationDate": None, "ModDate": None})
                    if preview_dir:
                        target = stage / "previews" / lang / (name + ".png")
                        target.parent.mkdir(parents=True, exist_ok=True)
                        fig.savefig(target, dpi=160, bbox_inches="tight")
                    plt.close(fig)
                    record(f"figures/{name}.pdf", dependencies, "figure", font_selection=fonts[lang][1])
            manifest = {"schema_version": 1, "generator": generator, "generated_at": timestamp, "language": lang,
                        "mode": "pilot" if allow_incomplete else "production", "results_dir": str(results_dir),
                        "inputs": inputs.sources, "issues": inputs.issues, "missing_macros": missing_macros,
                        "macros": macro_sources, "monte_carlo": mc_summary(inputs), "assets": assets,
                        "table_files": [a["path"] for a in assets if a["kind"] == "table"],
                        "figure_files": [a["path"] for a in assets if a["kind"] == "figure"],
                        "npz_policy": "No NPZ consumed; plotted arrays and intervals come from hashed JSON summaries."}
            write("revision_manifest.json", json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
            manifests[lang] = manifest
        # Atomic replacement per file; manifest is installed last as the run marker.
        for lang in ("en", "zh"):
            root = stage / f"paper_{lang}"
            paths = sorted((p for p in root.rglob("*") if p.is_file()), key=lambda p: (p.name == "revision_manifest.json", str(p)))
            for source in paths:
                target = repo_root / f"paper_{lang}" / source.relative_to(root)
                target.parent.mkdir(parents=True, exist_ok=True)
                temporary = target.with_suffix(target.suffix + ".rendering")
                temporary.write_bytes(source.read_bytes()); temporary.replace(target)
        if preview_dir:
            for source in (stage / "previews").rglob("*.png"):
                target = preview_dir / source.relative_to(stage / "previews")
                target.parent.mkdir(parents=True, exist_ok=True); target.write_bytes(source.read_bytes())
    return manifests


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results-dir", type=Path, help="Authoritative revised JSON directory (default: repo-root/experiment/analysis/revised)")
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2], help="Repository/output root; default is derived from this file")
    parser.add_argument("--allow-incomplete", action="store_true", help="Explicit PILOT only: label outputs and allow missing or pilot inputs")
    parser.add_argument("--preview-dir", type=Path, help="Optional external directory for PNG previews; default writes PDF only")
    args = parser.parse_args(argv)
    try:
        manifests = render_revision(args.results_dir, args.repo_root, args.allow_incomplete, args.preview_dir)
    except (RenderError, OSError, KeyError, TypeError) as exc:
        print(f"render_revision: {exc}", file=sys.stderr)
        return 2
    print(json.dumps({"mode": manifests["en"]["mode"], "macros": list(manifests["en"]["macros"]),
                      "missing_macros": manifests["en"]["missing_macros"],
                      "manifests": [str(args.repo_root.resolve() / f"paper_{lang}/revision_manifest.json") for lang in manifests]}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
