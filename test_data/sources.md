For LADCP specifically, there are a few excellent sources where you can get **both raw and processed** data (plus ancillary CTD/SADCP) to test a Python toolkit.

## Best all‑round source: GO‑SHIP / NCEI

1. **GO‑SHIP LADCP portal (processed profiles + figures)**  
   - Site: GO‑SHIP LADCP viewer at UH. [currents.soest.hawaii](https://currents.soest.hawaii.edu/go-ship/ladcp/)
   - What you get:  
     - Processed full‑depth profiles for many WOCE/CLIVAR/GO‑SHIP cruises.  
     - Per‑station plots, diagnostics tables, and NetCDF/Matlab/ASCII downloads.  
   - How it helps:  
     - Great for validating that your Python pipeline reproduces velocities and diagnostics for specific stations.  
     - Good “targets” to aim at before touching raw data.  

2. **NCEI GO‑SHIP LADCP datasets (raw + processed + ancillary)**  
   - Dataset overview: NCEI “GOSHIP‑LADCP” collection. [ncei.noaa](https://www.ncei.noaa.gov/metadata/granule/geoportal/rest/metadata/item/GOSHIP-LADCP.0221195/html)
   - Data structure: each dataset includes:  
     - Level‑0: raw LADCP (up/down cast) + CTD, navigation, shipboard ADCP. [accession.nodc.noaa](https://accession.nodc.noaa.gov/WOCE-LADCP)
     - Level‑1: processed LADCP (down, up, averaged) + ancillary. [accession.nodc.noaa](https://accession.nodc.noaa.gov/WOCE-LADCP)
     - Level‑2: science‑ready subsets of selected parameters. [accession.nodc.noaa](https://accession.nodc.noaa.gov/WOCE-LADCP)
   - How it helps:  
     - Exactly what you asked for: **paired raw and processed** data plus CTD and navigation for end‑to‑end testing.  
     - Ideal for checking file handling, transforms, inversion, etc.

3. **WOCE/CLIVAR archive index (cruise‑level LADCP packages)**  
   - UH LADCP data page: lists WOCE/CLIVAR cruises and points to NCEI. [currents.soest.hawaii](https://currents.soest.hawaii.edu/clivar/ladcp/)
   - Example: CLIVAR P02 cruise report describes a directory layout including `raw`, `CTD`, `SADCP`, `processed`, and `processed_noedit` for each cruise. [ldeo.columbia](https://www.ldeo.columbia.edu/~bhuber/LADCP/refs/reports/CLIVAR_PO2_Cruise_Report.pdf)
   - How it helps:  
     - Gives you “canonical” cruise structures you can mirror in your tests.  
     - `processed_noedit` vs. `processed` is perfect for exploring editing/QA decisions.

## Other useful public datasets

4. **Research Data Australia – Aurora Australis LADCP**  
   - Dataset: “LADCP current velocity data for CTD stations from the Australis 2006” (or similar). [researchdata.edu](https://researchdata.edu.au/ladcp-current-velocity-australis-2006/700644)
   - Notes: upward and downward `.adp` files, processed with MATLAB software. [researchdata.edu](https://researchdata.edu.au/ladcp-current-velocity-australis-2006/700644)
   - How it helps:  
     - Australian context, with documented processing; good to see how another group structured things.

5. **Manual / best‑practice cross‑checks**  
   - “A Manual for Acquiring Lowered Doppler Current Profiler (LADCP) Data” in Ocean Best Practices gives context and references to public software/datasets. [repository.oceanbestpractices](https://repository.oceanbestpractices.org/handle/11329/375)
   - Not a dataset, but useful for understanding how the raw/ancillary streams are expected to look.

## Where to start, concretely

If you want a minimal, realistic test set for a future Python pipeline:

1. Pick one GO‑SHIP cruise (e.g. a well‑documented P02 occupation) from the GO‑SHIP LADCP viewer. [currents.soest.hawaii](https://currents.soest.hawaii.edu/go-ship/ladcp/)
2. Grab the matching raw+processed package for that cruise from NCEI’s GOSHIP‑LADCP collection. [ncei.noaa](https://www.ncei.noaa.gov/metadata/granule/geoportal/rest/metadata/item/GOSHIP-LADCP.0221195/html)
3. Use the CLIVAR/WOCE cruise report to understand the directory conventions (`raw`, `CTD`, `SADCP`, `processed`). [ldeo.columbia](https://www.ldeo.columbia.edu/~bhuber/LADCP/refs/reports/CLIVAR_PO2_Cruise_Report.pdf)

That gives you a **closed system**:  
- Raw LADCP + CTD + navigation + SADCP.  
- Official processed outputs and diagnostics as your benchmark.  

Would you prefer something closer to your historical region (e.g. Indian/Southern Ocean sections) or are you happy to grab any well‑documented GO‑SHIP line as a first test set?

---

## Task 1 archive inspection (2026-06-20)

Three tgz archives were inspected in `test_data/cruise_data/`:

1. **processed_uv_20181105.tgz** (269 MB)
   - Contains: MATLAB files (001.mat – 030.mat), ASCII text files, and diagnostic PDFs (001_01.pdf – 030_14.pdf)
   - Raw PD0 files: NOT FOUND

2. **processed_uv_netcdf_20181105.tgz** (11 MB)
   - Contains: NetCDF files (001.nc – 030.nc)
   - Raw PD0 files: NOT FOUND

3. **processed_w_20181230.tgz** (3.0 MB)
   - Contains: Vertical velocity NetCDF profiles and PDFs (e.g., 001_VKEprof.nc, 001_prof.nc, 001_VKE.pdf)
   - Raw PD0 files: NOT FOUND

**Conclusion:** These archives contain **only processed outputs** (MATLAB, NetCDF, PDF) from cruise I7N (2018 Indian Ocean). Raw PD0 binary files (002DL000.000, 002UL000.000) are **not present** in any archive. Unit tests do not require raw files and will proceed with mock/synthetic data as needed.