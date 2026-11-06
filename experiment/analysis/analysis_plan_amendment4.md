# Analysis plan, Amendment 4, 2026-11-06

Paper: "What more cases buy: inferential targets for reproduction-number estimates and what cross-sectional dispersion diagnostics measure".

This amendment belongs to `analysis_plan.md`, `analysis_plan_amendment2.md` and `analysis_plan_amendment3.md`. It changes the scope of the shipped paper, not any estimate. Its sha256 and UTC time are recorded in `analysis_plan_amendment4.lock`.

## A4.1 Withdrawal of A3.4 (cross-series replication)

Section A3.4 of Amendment 3 ran the same diagnostic on NHSN COVID-19 and RSV admissions and reported the results as a replication check. Those results were computed and are recorded in `analysis/rev_r3_crosspathogen.json`; the figures they produced are:

- COVID-19 passed the window rule in all four seasons (1,614 rows) and had affine intercept 0.0056 (0.0035--0.0077), model-free floor 7.1% (5.8--8.1%) and population coefficient $-0.0012$ ($-0.0023$ to $-0.0002$).
- RSV passed the window rule in two seasons (720 rows) and had affine intercept 0.0113 (0.0027--0.0199), model-free floor 9.6% (5.2--12.2%) and population coefficient $-0.0029$ ($-0.0067$ to 0.0010).

The paper is now restricted to a single empirical series, NHSN influenza hospital admissions, and Section A3.4 is withdrawn from the manuscript. The reasons are recorded here rather than in the paper:

1. The replication criterion of A3.4 fixed the floor as replicated only if its point estimate lay inside the influenza 95% interval. Under the alternative criterion "affine intercept positive with an interval excluding zero", both COVID-19 and RSV would have counted as replicating the floor. The criterion, not the data, decided the verdict, and a single-series paper does not need to carry that choice.
2. The RSV sample is 720 rows over two seasons, with 2022--23 carrying no national counts and 2023--24 failing the reporting-coverage floor. It is too small to separate non-replication from low power.
3. The paper's argument concerns what a cross-sectional dispersion diagnostic measures, not whether a particular pathogen replicates. A cross-series section that partly succeeds and partly fails invites the reading that the diagnostic is pathogen-specific, which the single series cannot settle either way.

## A4.2 What is unaffected

Every quantity in the influenza analysis is unchanged. The primary test (E1), the sensitivity set (16 independent analyses), the population gradient under the free count exponent, the bias range over the simulated mechanisms, the reference-conditioned interval for $b/b_P$ and the season test are computed from the influenza data only and are therefore untouched by this withdrawal.

## A4.3 Record keeping

The code `code/rev_r3_crosspathogen.py`, the output `analysis/rev_r3_crosspathogen.json` and its log remain in the repository as part of the record of what was run and what it produced. They are not run by `run_all.sh` and no number from them appears in the paper. The `revision_2026_11_r3` block of `results.json` drops the `A3.4_crosspathogen` key; the key `A3.5_tables` and lower keep their names.

## A4.4 Decision rule, fixed here

If a later revision reinstates a cross-series check, it must state the replication criterion before the run, report the verdict under both the interval-membership and the interval-excludes-zero criteria, and state the sample size of every series that fails the window rule. No analysis in this amendment is confirmatory.
