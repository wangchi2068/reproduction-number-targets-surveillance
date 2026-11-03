"""s99: merge analysis/parts/*.json into analysis/results.json with per-number provenance.

Every number carries: value, optional ci95, script, inputs, seed (number-specific if recorded, else the
part's base seed), mode (measured = computed from NHSN/Delphi data; simulated = Monte Carlo; analytic =
closed form or arithmetic on fitted values). Keys are '<part>.<key>'. Detailed blocks stay in the part
files; results.json lists them under 'blocks' by name so figures and tables can be traced.
"""
import hashlib
import json
import time

from cc_common import ANALYSIS, DATA, PARTS, SEED, sha256

out = {"meta": {"generated": time.strftime("%Y-%m-%dT%H:%M:%S"), "base_seed": SEED,
                "data": {"nhsn_vintage": {"file": "data/external/nhsn_hrd_vintage_2026-10-20.csv.gz",
                                          "sha256": sha256(DATA / "external/nhsn_hrd_vintage_2026-10-20.csv.gz")},
                         "delphi_fluview_clinical": {"file": "data/external/delphi_fluview_clinical_raw.json",
                                                     "sha256": sha256(DATA / "external/delphi_fluview_clinical_raw.json")}},
                "interval": "95% percentile interval of the two-way (state x week) pigeonhole bootstrap unless a "
                            "number's note says otherwise",
                "parts": {}},
       "numbers": {}, "blocks": {}}
R3_PARTS = ("inversion", "multiscale", "reporting_regime", "returns", "parag", "polish")  # round 3 (s09) and final polish (s10)
_hash = {}


def input_hashes(paths):
    out = {}
    for pth in paths or []:
        fp = DATA.parent / pth
        if fp.is_file():
            _hash.setdefault(pth, sha256(fp))
            out[pth] = _hash[pth]
    return out


for f in sorted(PARTS.glob("s*.json")) + [PARTS / f"{n}.json" for n in R3_PARTS if (PARTS / f"{n}.json").exists()]:
    d = json.loads(f.read_text())
    name = f.stem
    meta = d.get("meta", {})
    out["meta"]["parts"][name] = {k: meta.get(k) for k in ("script", "inputs", "mode", "seed", "updated")}
    out["meta"]["parts"][name]["sha256_part"] = hashlib.sha256(f.read_bytes()).hexdigest()
    for k, rec in d.get("numbers", {}).items():
        r = dict(rec)
        r.setdefault("mode", meta.get("mode"))
        r.setdefault("seed", meta.get("seed"))
        r["script"] = meta.get("script")
        r["inputs"] = meta.get("inputs")
        if name in R3_PARTS:
            r["inputs_sha256"] = input_hashes(meta.get("inputs"))
        out["numbers"][f"{name}.{k}"] = r
    out["blocks"][name] = sorted(d.get("blocks", {}))
(ANALYSIS / "results.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
print(len(out["numbers"]), "numbers from", len(out["meta"]["parts"]), "parts")
