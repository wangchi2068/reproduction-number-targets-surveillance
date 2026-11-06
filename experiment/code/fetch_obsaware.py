"""Fetch the ObsAwareData release history used by rev_measure_revision.py.

Downloads one `<ST>_long_all_versions.csv.gz` per jurisdiction (51 files, about 0.9 GB in total) from
the DMA-PRIME/ObsAwareData GitHub releases, verifying each file against the SHA-256 listed in the
repository's metadata/release_assets.csv.

Source: https://github.com/DMA-PRIME/ObsAwareData (CC BY 4.0)

Usage: python code/fetch_obsaware.py [--dest DIR] [--states AK,AL,...]
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import rev_common as R

META = "https://raw.githubusercontent.com/DMA-PRIME/ObsAwareData/main/metadata/release_assets.csv"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dest", default=str(R.EXP / "data" / "obsaware"))
    ap.add_argument("--states", default="", help="comma-separated postal codes; default all")
    args = ap.parse_args()
    dest = Path(args.dest)
    dest.mkdir(parents=True, exist_ok=True)

    with urllib.request.urlopen(META, timeout=120) as fh:
        rows = list(csv.DictReader(io.TextIOWrapper(fh, encoding="utf-8")))
    full = [r for r in rows if r["file"].endswith("_long_all_versions.csv.gz") and "added" not in r["file"]]
    want = {s.strip().upper() for s in args.states.split(",") if s.strip()}
    if want:
        full = [r for r in full if r["state"] in want]
    print(f"{len(full)} file(s) to fetch")

    for r in full:
        out = dest / f"{r['state']}.csv.gz"
        if out.exists() and sha256(out) == r["sha256"]:
            print(f"{r['state']} cached")
            continue
        with urllib.request.urlopen(r["url"], timeout=900) as fh, open(out, "wb") as w:
            w.write(fh.read())
        got = sha256(out)
        status = "ok" if got == r["sha256"] else "HASH MISMATCH"
        print(f"{r['state']} {out.stat().st_size} bytes {status}")
        if got != r["sha256"]:
            raise SystemExit(f"{r['state']}: expected {r['sha256']}, got {got}")


if __name__ == "__main__":
    main()
