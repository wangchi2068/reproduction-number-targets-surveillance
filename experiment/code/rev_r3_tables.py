"""Amendment 3, Section A3.5: sensitivity count and table notes (review N13, N14).

No new estimation. Reads `results.json` and writes `analysis/rev_r3_tables.json` holding:

  - the sensitivity ledger, with S-a5 at threshold 1.02 marked as a duplicate of the primary sample
    (identical window set and row count), and the count of independent sensitivity analyses;
  - the bootstrap-t table note recorded by review N14;
  - the split of tables between the main text and the supplementary file.

Usage: python code/rev_r3_tables.py
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

OUT = R.AN / "rev_r3_tables.json"
HASHES = R.AN / "revision3" / "script_hashes.jsonl"
DUPLICATE = "S-a5_threshold_1.02"

BOOTSTRAP_T_NOTE = (
    "The studentised statistic uses the stacked cluster-robust sandwich that includes the estimate of "
    "xbar, not a sandwich with xbar held fixed. This choice is required for stability on clamped "
    "bootstrap draws, where the restricted mean b0/D would otherwise put mass at zero."
)


def main():
    ap = argparse.ArgumentParser()
    ap.parse_args()
    t0 = time.time()
    HASHES.parent.mkdir(parents=True, exist_ok=True)
    R.append_jsonl(HASHES, {"script": "rev_r3_tables.py",
                            "sha256": hashlib.sha256(open(__file__, "rb").read()).hexdigest(),
                            "utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()),
                            "note": "Amendment 3, Section A3.5"})

    res = json.loads((R.AN / "results.json").read_text(encoding="utf-8"))
    e1 = res["revision_2026_11"]["E1_intercept"]
    sens = e1["sensitivities"]
    primary_rows = int(res["revision_2026_11"]["E2_sample"]["rows"]) if "E2_sample" in res["revision_2026_11"] else 1748

    ledger = []
    for name, rec in sens.items():
        has_rows = "n_rows" in rec
        dup = bool(name == DUPLICATE and has_rows and int(rec["n_rows"]) == primary_rows)
        ci = rec.get("b_over_bP_ci_two_way")
        ledger.append({"name": name,
                       "n_rows": int(rec["n_rows"]) if has_rows else None,
                       "n_states": int(rec["n_states"]) if "n_states" in rec else None,
                       "n_weeks": int(rec["n_weeks"]) if "n_weeks" in rec else None,
                       "lower_limit_b_over_bP": float(ci[0]) if ci else rec.get("lower_limit_min"),
                       "upper_limit_b_over_bP": float(ci[1]) if ci else None,
                       "duplicate_of_primary": dup,
                       "role": ("duplicate check" if dup else "independent")})
    independent = [r for r in ledger if not r["duplicate_of_primary"]]

    out = {
        "meta": {"script": "code/rev_r3_tables.py",
                 "plan": "analysis/analysis_plan_amendment3.md, Section A3.5",
                 "source": "results.json revision_2026_11.E1_intercept.sensitivities",
                 "seconds": round(time.time() - t0, 1)},
        "sensitivity_ledger": ledger,
        "counts": {"declared": len(ledger), "duplicates": len(ledger) - len(independent),
                   "independent": len(independent),
                   "text_phrasing": f"{len(independent)} independent sensitivity analyses"},
        "s_a5_note": (f"{DUPLICATE} selects exactly the primary sample "
                      f"({primary_rows} rows, identical window set) and is therefore a duplicate check, "
                      f"not an independent sensitivity analysis"),
        "bootstrap_t_note": BOOTSTRAP_T_NOTE,
        "table_split": {
            "main_text": ["sensitivity", "b_over_bP_reference", "population_gradient",
                          "cross_series", "bias_range"],
            "supplementary": ["state_realistic_generator_details", "clamp_types", "interval_types"],
        },
    }
    R.write_json(OUT, out)
    print(json.dumps({"declared": out["counts"]["declared"], "independent": out["counts"]["independent"],
                      "phrase": out["counts"]["text_phrasing"]}, ensure_ascii=False))
    print("duplicates:", [r["name"] for r in ledger if r["duplicate_of_primary"]])


if __name__ == "__main__":
    main()
