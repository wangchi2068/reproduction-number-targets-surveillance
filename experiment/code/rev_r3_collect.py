"""Amendment 3, Section A3.6: collect the revision-3 results into `results.json`.

Reads the five section outputs written by the revision-3 lane and adds a single `revision_2026_11_r3`
block with the numbers the manuscript cites, the decision-rule verdicts and the script hashes. No
estimation happens here; every value is copied from its section file.

Usage: python code/rev_r3_collect.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import rev_common as R

PLAN = "analysis/analysis_plan_amendment3.md"
HASHES = R.AN / "revision3" / "script_hashes.jsonl"
SECTIONS = {
    "A3.1_population": "rev_r3_population.json",
    "A3.2_reference": "rev_r3_reference.json",
    "A3.3_bias_range": "rev_r3_bias_range.json",
    "A3.5_tables": "rev_r3_tables.json",
}


def main():
    ap = argparse.ArgumentParser()
    ap.parse_args()
    t0 = time.time()
    R.append_jsonl(HASHES, {"script": "rev_r3_collect.py",
                            "sha256": hashlib.sha256(open(__file__, "rb").read()).hexdigest(),
                            "utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()),
                            "note": "Amendment 3, Section A3.6"})

    src = {k: json.loads((R.AN / v).read_text(encoding="utf-8")) for k, v in SECTIONS.items()}

    pop = src["A3.1_population"]
    ref = src["A3.2_reference"]
    bias = src["A3.3_bias_range"]
    tabs = src["A3.5_tables"]

    block = {
        "generated_utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()),
        "analysis_plan": {"amendment3": PLAN, "lock": "analysis/analysis_plan_amendment3.lock"},
        "data": "frozen NHSN 2026-10-20 vintage",
        "A3.1_population": {
            "M1_a1": pop["M1_affine_population"]["a1"],
            "M2_a1": pop["M2_free_exponent_population"]["a1"],
            "M2_gamma": pop["M2_free_exponent_population"]["gamma"],
            "M1_a0": pop["M1_affine_population"]["a0"],
            "M1_ratio_p90_over_p10": pop["M1_affine_population"]["ratio_p90_over_p10"]["point"],
            "M2_ratio_p90_over_p10": pop["M2_free_exponent_population"]["ratio_p90_over_p10"]["point"],
            "bootstrap": {k: {"ci": v["ci"], "n": v["n"]} for k, v in pop["bootstrap"].items()},
            "matched_count_table": pop["matched_count_table"],
            "verdict": pop["decision_rule"]["verdict"],
        },
        "A3.2_reference": {
            "point": ref["point"],
            "interval": ref["interval"],
            "point_range": ref["point_range"],
            "lower_limit_exceeds_1_all_references": ref["lower_limit_exceeds_1_all_references"],
            "B": ref["meta"]["B"],
            "verdict": ("lower limit exceeds 1 under all references"
                        if ref["lower_limit_exceeds_1_all_references"]
                        else "whether b/b_P exceeds 1 depends on the reference"),
        },
        "A3.3_bias_range": {
            "range_pct": [bias["range"]["min_pct"], bias["range"]["max_pct"]],
            "rounded_pct": bias["range"]["rounded_pct"],
            "fresh_draw_cells": bias["fresh_draw_cells"],
            "state_realistic_note": bias["state_realistic_note"],
        },
        "A3.5_tables": {
            "counts": tabs["counts"],
            "s_a5_note": tabs["s_a5_note"],
            "bootstrap_t_note": tabs["bootstrap_t_note"],
            "table_split": tabs["table_split"],
        },
        "scripts": [json.loads(line) for line in HASHES.read_text(encoding="utf-8").splitlines() if line.strip()]
                   if HASHES.exists() else [],
    }

    path = R.AN / "results.json"
    res = json.loads(path.read_text(encoding="utf-8"))
    res["revision_2026_11_r3"] = block
    path.write_text(json.dumps(res, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"r3_written": True, "A3.1": pop["decision_rule"]["verdict"],
                      "A3.2_lower_all": ref["lower_limit_exceeds_1_all_references"],
                      "A3.3_range": bias["range"]["rounded_pct"]}, ensure_ascii=False, indent=1))
    print("seconds", round(time.time() - t0, 1))


if __name__ == "__main__":
    main()
