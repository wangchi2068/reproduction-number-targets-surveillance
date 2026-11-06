"""Measured reporting noise from the NHSN revision history (ObsAwareData).

Each state-week enters every CDC release with the value that release published. The spread of those
values around the final value is a direct measurement of the reporting/recording error in a settled
count: it contains no transmission variation, because the week and the jurisdiction are fixed. This
is the quantity the model calls reporting noise, and Proposition 2 says it cannot be recovered from a
single release of the data.

For every state-week with at least two published values we record:

    v_rev = mean over releases of ((C_release - C_final) / C_final)^2

and the affine diagnostic's count is C_final. The first release (earliest as_of_date) is used for the
"provisional" variant, which is what a real-time analyst would have seen.

Output: analysis/revision3/measured_revision.json
Usage: python code/rev_measure_revision.py [--root DIR] [--pathogen influenza]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import statistics
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import rev_common as R

FEATURE = "hospital_admission_count"
VALID = {"observed_valid", "observed_zero"}


def read_state(path, pathogen):
    """Yield (state, week, as_of, value) for valid admission counts."""
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row["pathogen"] != pathogen or row["feature_name"] != FEATURE:
                continue
            if row["observation_state"] not in VALID:
                continue
            try:
                value = float(row["weekly_feature"])
            except (TypeError, ValueError):
                continue
            yield (row["raw_geography_id"], row["observation_week"], row["as_of_date"], value)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(R.EXP / "data" / "obsaware"))
    ap.add_argument("--pathogen", default="influenza")
    args = ap.parse_args()
    t0 = time.time()
    root = Path(args.root)

    series = defaultdict(dict)          # (state, week) -> {as_of: value}
    states = 0
    for path in sorted(root.glob("*.csv.gz")):
        states += 1
        for state, week, as_of, value in read_state(path, args.pathogen):
            series[(state, week)][as_of] = value

    rows = []
    for (state, week), releases in series.items():
        if len(releases) < 2:
            continue
        ordered = sorted(releases.items())
        final = ordered[-1][1]
        first = ordered[0][1]
        if final <= 0:
            continue
        sq = [(v - final) ** 2 / final ** 2 for _, v in ordered]
        rows.append({"state": state, "week": week, "n_releases": len(releases),
                     "final": final, "first": first,
                     "v_rev": statistics.fmean(sq),
                     "first_rel": (first - final) / final,
                     "all_equal": all(abs(v - final) < 1e-12 for _, v in ordered)})

    out = {
        "meta": {"script": "code/rev_measure_revision.py",
                 "source": "ObsAwareData release history, NHSN hospital_admission_count",
                 "pathogen": args.pathogen, "states_read": states,
                 "states_with_data": len({r["state"] for r in rows}),
                 "state_weeks": len(rows),
                 "definition": "v_rev = mean_r ((C_r - C_final)/C_final)^2 over releases r",
                 "seconds": round(time.time() - t0, 1)},
        "summary": {
            "share_never_revised": (sum(1 for r in rows if r["all_equal"]) / len(rows)) if rows else None,
            "median_v_rev": statistics.median([r["v_rev"] for r in rows]) if rows else None,
            "mean_v_rev": statistics.fmean([r["v_rev"] for r in rows]) if rows else None,
            "rms_relative_revision": (statistics.fmean([r["v_rev"] for r in rows]) ** 0.5) if rows else None,
            "median_n_releases": statistics.median([r["n_releases"] for r in rows]) if rows else None,
        },
        "rows": rows,
    }
    dest = R.AN / "revision3" / "measured_revision.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps(out["meta"], indent=1))
    print(json.dumps(out["summary"], indent=1))


if __name__ == "__main__":
    main()
