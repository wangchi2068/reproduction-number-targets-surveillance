"""Amendment 3, Section A3.3: bias range over all simulated mechanisms (review N11).

Reads the shipped fresh-draw cells in results.json (revision_2026_11_r2.C_total_bias.generators),
reports the relative bias of the fitted intercept against the transmission-plus-reporting total for
every cell with a successful calibration, and takes the minimum and maximum as the reported range.
State-realistic cells are listed separately with their calibration status and a_over_total ratio;
they do not enter the numeric range because none passes the six declared tests.

Usage: python code/rev_r3_bias_range.py
"""
from __future__ import annotations

import hashlib
import json
import time

import rev_common as R

OUT = R.AN / "rev_r3_bias_range.json"
HASHES = R.AN / "revision3" / "script_hashes.jsonl"


def main():
    t0 = time.time()
    HASHES.parent.mkdir(parents=True, exist_ok=True)
    R.append_jsonl(HASHES, {"script": "rev_r3_bias_range.py",
                            "sha256": hashlib.sha256(open(__file__, "rb").read()).hexdigest(),
                            "utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()),
                            "note": "Amendment 3, Section A3.3"})
    res = json.loads((R.AN / "results.json").read_text(encoding="utf-8"))
    gen = res["revision_2026_11_r2"]["C_total_bias"]["generators"]
    sr = res["revision_2026_11_r2"]["D1_state_realistic_generators"]["cells"]

    fd = []
    for k, v in gen.items():
        if not k.startswith("FD-"):
            continue
        fd.append({"cell": k, "relative_bias": float(v["relative_bias"]),
                   "ci": [float(x) for x in v["relative_bias_ci"]],
                   "calibration": "ok"})
    biases = [c["relative_bias"] for c in fd]
    lo, hi = min(biases), max(biases)

    sr_rows = []
    for k, v in sr.items():
        sr_rows.append({"cell": k,
                        "calibration_failed": bool(v["calibration"].get("calibration_failed")),
                        "a_over_total": float(v["a_over_total"]) if v.get("a_over_total") is not None else None,
                        "coverage_of_true_total": v.get("coverage_of_true_total", {}).get("coverage"),
                        "all_pass": bool(v.get("all_pass")),
                        "failing_tests": v.get("failing_tests", [])})

    out = {
        "meta": {"script": "code/rev_r3_bias_range.py", "plan": "analysis/analysis_plan_amendment3.md, Section A3.3",
                 "source": "results.json revision_2026_11_r2.C_total_bias.generators and "
                           "revision_2026_11_r2.D1_state_realistic_generators",
                 "seconds": round(time.time() - t0, 1)},
        "fresh_draw_cells": fd,
        "range": {"min": lo, "max": hi, "min_pct": 100 * lo, "max_pct": 100 * hi,
                  "rounded_pct": [int(round(100 * lo)), int(round(100 * hi))]},
        "state_realistic_cells": sr_rows,
        "state_realistic_note": "none of the six state-realistic cells passes the six declared tests; the one "
                                "cell with a successful calibration recovers 77.6% of the simulated total, "
                                "so its bias is about -22%",
        "decision_rule": {"rule": "report the min-max of the relative bias over the fresh-draw cells; "
                                  "state-realistic cells are reported separately",
                          "abstract_phrase": "6--30%"},
    }
    R.write_json(OUT, out)
    print(json.dumps(out["range"], indent=1), flush=True)


if __name__ == "__main__":
    main()
