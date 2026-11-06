"""Verify the frozen inputs against the SHA-256 values recorded in data/DATA_MANIFEST.md (run first by run_all.sh).

Inputs read by the pipeline must match exactly (exit 4 otherwise). Derived or legacy files are reported but not
enforced: the China extract and its script now live under legacy/china_covid19/ and no pipeline step reads them.
Usage: python code/rev_verify_data.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

D = Path(__file__).resolve().parents[1] / "data"
REQUIRED = {
    "external/nhsn_hrd_vintage_2026-10-20.csv.gz": "9685ddd4b0035694754334e530d5b356bd3503f57b82821f7ee7883c4f5fb94b",
    "external/census_NST-EST2025-ALLDATA.csv": "92188e29cb0a67dcf95afa7d6c47359409782f086478b70ea4128eb70e223ca9",
    "external/delphi_nhsn/MANIFEST.json": "fab6e1742518f37453bdf96e523359997d2285ab2a8db73d47a217a4378dbd9a",
}
REPORTED = {
    "external/jhu_covid19_china_cases.csv": "4004ab1c36c56804b9ff4d0ceaef94042f17426d3f4df07f0430c147b83316f8",
    "external/delphi_fluview_clinical_raw.json": "7c14547acead56bf52333fad18449efc13485d67b0db3f635cb573ae7283651a",
    "panel_nhsn_vintage_2026-10-20.csv.gz": "244a025f9833c54bd07a3fbb1d122e19dbe696e89f8e33aa022dfb6fa9d86c7a",
    "national_flu.csv": "12acdc02e0b27754ac2d5ed5ddc258770870a9a5e5861ac9f86ae26a47b91c03",
}


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    bad = 0
    for group, table in (("required", REQUIRED), ("reported", REPORTED)):
        for name, want in table.items():
            p = D / name
            got = sha(p) if p.exists() else "missing"
            ok = got == want
            print(f"{group:8s} {'OK  ' if ok else 'DIFF'} {name} {got}")
            bad += (not ok) and group == "required"
    man = D / "external/delphi_nhsn/MANIFEST.json"
    if man.exists():  # the 51 Delphi issue-history files are checked against their (checksummed) manifest
        files = json.loads(man.read_text())["files"]
        diff = [n for n, e in files.items() if not (man.parent / n).exists() or sha(man.parent / n) != e["sha256"]]
        print(f"required {'OK  ' if not diff else 'DIFF'} external/delphi_nhsn/*.json.gz ({len(files)} files; {len(diff)} differ)")
        bad += bool(diff)
    print("verify_data", "PASS" if not bad else f"FAIL ({bad} required file(s) differ)")
    sys.exit(4 if bad else 0)


if __name__ == "__main__":
    main()
