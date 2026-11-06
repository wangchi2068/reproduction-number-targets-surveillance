# ObsAwareData revision history

Source: https://github.com/DMA-PRIME/ObsAwareData (CC BY 4.0)

Per-state release history of NHSN hospital admissions. Each state-week appears once per CDC release,
with the value that release published, so the spread of values measures reporting/recording error
directly (no transmission variation, because the week and jurisdiction are fixed).

Not redistributed in this repository (887 MB). To fetch:

    python code/fetch_obsaware.py

Files land in this directory as `<ST>.csv.gz`. Consumers: `code/rev_measure_revision.py`.
