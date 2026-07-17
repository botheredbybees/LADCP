# Task brief: bulk validation runs for I7N 2018 and A16N 2013

**Audience:** an autonomous Claude session (any model; written for Haiku).
**Mission:** download the remaining raw LADCP casts for two cruises from
NCEI, run the existing validation harness over every cast, and produce a
results report. This is data-gathering and harness-running ONLY.

## Ground rules (binding)

1. Do **NOT** modify anything under `src/`, `scripts/`, `tests/`, `docs/`,
   or `octave_harness/`. No exceptions, including "small fixes".
2. If the harness errors on a cast, **record the error text and move on**
   to the next cast. Do not debug, do not work around, do not retry more
   than once. Casts that error are a *result*, not a problem to solve.
3. If something structural blocks all progress (site down, disk full,
   harness won't start at all), STOP and write what you observed to
   `BULK_VALIDATION_REPORT.md`. Do not improvise.
4. All new files go under `test_data/2018_I7N/`, `test_data/2013_A16N/`,
   or `BULK_VALIDATION_REPORT.md` at the repo root. Nothing else.
5. Everything is resumable: before downloading any file, check whether it
   already exists with size > 100 KB — if so, skip the download. Before
   running a batch, check whether its CSV already exists — if so, skip
   the batch. A later session can pick up where a killed one stopped.
6. RuntimeWarnings ("Mean of empty slice" etc.) and CTD-lag warnings in
   harness output are normal. Only `ERROR` rows in the per-cast table
   count as failures.

## Before starting

- Working directory: `C:\Users\peter_sha\Documents\sourcecode\Nuyina\LADCP`
- Check free disk space; need ~10 GB. If less, STOP and report.
- Downloads use plain `curl -sf --retry 3 -o <dest> <url>`. If curl exits
  nonzero, DELETE the partial destination file (truncated files corrupt
  results silently — see `test_data/2018_I7N/DOWNLOAD_NOTES.md`, cast 020),
  record the cast as `DOWNLOAD-FAILED`, and continue.
- Run all long commands in the background and check output files; each
  cast takes ~2–3 minutes in the harness, so a full-cruise batch sequence
  runs for hours. That is expected. Work serially; no parallelism needed.

## Cruise 1: I7N 2018 (124 reference casts)

Reference NCs: `test_data/2018_I7N/processed_uv/001.nc … 124.nc`.
Raw archive root (HTTPS, browsable):

    https://www.ncei.noaa.gov/data/oceans/archive/arc0163/0222105/1.1/data/0-data/2018_I07N/raw-level0/

Per cast `NNN` (3-digit, e.g. `047`), three files go into `test_data/2018_I7N/raw/`:

| file | from |
|---|---|
| `NNNDL000.000` | `LADCP_raw/NNN/NNNDL000.000` |
| `NNNUL000.000` | `LADCP_raw/NNN/NNNUL000.000` |
| `NNN_01.cnv`   | `CTD_24Hz/NNN_01.cnv` |

If a URL 404s, first fetch the directory listing
(`curl -s <root>/LADCP_raw/NNN/`) and use whichever file matches
`*DL*.000` / `*UL*.000` (alphabetically first if several); same idea for
`CTD_24Hz/` with `NNN_*.cnv`. If still nothing, record the cast as
`MISSING-ON-SERVER` and continue.

Casts 003 and 010 are already on disk. Cast 020's CTD is truncated —
delete `test_data/2018_I7N/raw/020_01.cnv` and re-download it.

**Batching:** process casts in batches of 10 (001–010, 011–020, …,
121–124). For each batch:

1. Download the batch's missing files (skip-if-exists rule above).
2. Run (single command, from the repo root):

       TEST_DATA_DIR=test_data uv run python scripts/validate_multicast.py --raw-dir test_data/2018_I7N/raw --ref-dir test_data/2018_I7N/processed_uv --casts <comma-separated batch casts> --out test_data/2018_I7N/bulk/batch_NN.csv > test_data/2018_I7N/bulk/batch_NN.log 2>&1

   (create `test_data/2018_I7N/bulk/` first; NN = batch number, zero-padded)

After all batches: concatenate the batch CSVs into
`test_data/2018_I7N/validate_all.csv` (keep only the first header line).

## Cruise 2: A16N 2013 (casts 001–095 only)

Casts 096+ have no uplooker and CANNOT run in this harness — do not
download them. Reference NCs: `test_data/2013_A16N/processed_nc/`.
Raw archive root:

    https://www.ncei.noaa.gov/data/oceans/archive/arc0147/0205839/1.1/data/0-data/2013.rb1304.A16N/raw-level0/

Exact filenames are already listed on disk — use these instead of
scraping (each line is a filename):

- `test_data/2013_A16N/listing_ladcp_raw_downlooker_WH150.txt` → prefix URL path `ladcp_raw/downlooker/WH150/`
- `test_data/2013_A16N/listing_ladcp_raw_uplooker_WH300.txt` → `ladcp_raw/uplooker/WH300/`
- `test_data/2013_A16N/listing_ctd_timeseries.txt` → `ctd_timeseries/`

For each cast `NNN` in 001–095, the three files contain `00NNN_` in their
names (5-digit station). Download into `test_data/2013_A16N/raw/`.
Casts 003, 010, 030, 060, 090 are already on disk.

**Batching:** batches of 10 as above, output to
`test_data/2013_A16N/bulk/batch_NN.csv`, using:

       TEST_DATA_DIR=test_data uv run python scripts/validate_multicast.py --raw-dir test_data/2013_A16N/raw --ref-dir test_data/2013_A16N/processed_nc --casts <batch> --dl-glob "DLWH*00{cast}_*.PD0" --ul-glob "ULWH*00{cast}_*.PD0" --ctd-glob "ctd_timeseries_00{cast}_*_gps.txt" --out test_data/2013_A16N/bulk/batch_NN.csv > test_data/2013_A16N/bulk/batch_NN.log 2>&1

After all batches: concatenate into `test_data/2013_A16N/validate_all.csv`.

## Deliverable: BULK_VALIDATION_REPORT.md (repo root)

Write (or update, if resuming) a report containing, per cruise:

1. Counts: casts attempted / succeeded / ERROR / SKIP-missing-files /
   DOWNLOAD-FAILED / MISSING-ON-SERVER.
2. Pass counts: how many succeeded casts have `u_rmse < 0.05` and
   `v_rmse < 0.05` (both), and how many pass u only.
3. Distribution: median / p90 / worst u_rmse and v_rmse over succeeded casts.
4. A table of every ERROR cast with the one-line exception text from the log.
5. A table of the 10 worst casts by u_rmse (cast, n, u_rmse, v_rmse, r_u, max_depth).
6. Anything anomalous you noticed, stated as observations only — no fixes,
   no code suggestions.

Do NOT commit anything to git. Leave that to the main session.

## Expected scale (sanity reference)

- I7N: ~122 casts to download ≈ 5.5 GB, ≈ 5–6 h of harness time.
- A16N: ~90 casts ≈ 1.2 GB, ≈ 3–4 h of harness time.
- Prior spot-check results, for calibration: I7N 003/010 pass u (<0.05);
  A16N 003/010 pass both; A16N deep casts (>4 km) are KNOWN to fail with
  u ~0.1–0.45 — when you see those, that is the expected result, not a
  mistake in your setup.
