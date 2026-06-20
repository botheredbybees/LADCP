# Project Proposal: A Modern Python Toolkit for LADCP Processing

## Overview

This proposal outlines a future software project to build a modern, open, Python-based toolkit for processing Lowered Acoustic Doppler Current Profiler (LADCP) data. The need remains clear because LADCP is still one of the main ways to obtain direct, full-water-column current profiles during CTD operations, while post-processing remains technically demanding and heavily dependent on legacy workflows (Thurnherr, n.d.; University of Hawaiʻi at Mānoa, 2015).

The project would focus on the parts of the workflow that have historically been awkward for operators and analysts: coordinate transforms, heading and tilt correction, synchronization with CTD and navigation data, inverse solutions, diagnostic plotting, and reproducible provenance. The goal is not to replace the scientific basis of established methods, but to make those methods easier to inspect, reproduce, validate, and extend in modern Python-first environments (Thurnherr, n.d.; University of Hawaiʻi at Mānoa, 2015).

## Background

Standard Argo floats are extremely valuable for temperature and salinity profiling, but they do not provide the same kind of direct full-depth current profile that a well-processed LADCP cast can provide. Standard Argo missions infer park-depth currents from float displacement, whereas LADCP systems can estimate velocity structure through the full water column during a cast (Integrated Marine Observing System, 2025; International Pacific Research Center, n.d.; NOAA Atlantic Oceanographic and Meteorological Laboratory, 2024).

That difference is one reason LADCP remains important for ship-based deep ocean programs. At the same time, LADCP processing is still widely regarded as non-trivial, with mature but complex software lineages that depend on MATLAB, cruise-specific conventions, and operator experience (British Oceanographic Data Centre, n.d.; Thurnherr, n.d.; University of Hawaiʻi at Mānoa, 2015).

## Note on Perl in the workflow

The recollection of using Perl is well supported by the historical tooling around LADCP operations. Andreas Thurnherr notes that the portable LDEO acquisition system used public-domain tools, including *expect* scripts and a Perl-based communication program called *bbabble* to talk to one or two RDI ADCPs, while other cruise documentation describes UNIX shell-command based acquisition wrappers around onboard LADCP operations (Thurnherr, n.d.; Cruise Report P02E, 2022).

That detail matters because the historical toolchain was often a hybrid of MATLAB, shell scripts, Perl-based communication layers, and cruise-specific helper utilities rather than a single coherent software package. A modern proposal should therefore account for both acquisition-time and post-processing-time complexity, including interfaces with inherited scripts and file conventions (Thurnherr, n.d.; Cruise Report P02E, 2022).

## Problem statement

The core practical problem is not that no software exists; it is that the available software is fragmented, expert-centric, and difficult to integrate into modern reproducible scientific computing workflows. The official and semi-official LADCP processing pathways are scientifically credible, but they are distributed across aging documentation, MATLAB scripts, GitHub repositories, and institution-specific operational knowledge (Thurnherr, n.d.; British Oceanographic Data Centre, n.d.; GitHub, 2014).

This is especially awkward for users working in Linux, containers, Python, version control, and automated data pipelines. It also creates unnecessary friction for long-term reprocessing, auditability, team handover, and method comparison across cruises or institutions (Thurnherr, n.d.; GitHub, 2014).

## Project goals

The project should aim to deliver the following capabilities:

- A Python package for ingesting LADCP raw files together with associated CTD, GPS, and navigation inputs.
- Explicit implementations of beam-to-instrument, instrument-to-ship, and ship-to-earth coordinate transforms (apaloczy, 2018; SourceForge, 2021).
- Modular correction stages for heading, tilt, rotation, sound-speed handling, and timing offsets (Thurnherr, n.d.; University of Hawaiʻi at Mānoa, 2015).
- Support for both classic shear-based approaches and inverse or velocity-based solution methods used in established workflows (Thurnherr, n.d.; JAMSTEC, n.d.).
- Diagnostic plots, cast-level quality-control summaries, and machine-readable processing provenance (University of Hawaiʻi at Mānoa, 2015).
- Containerized execution suitable for laptops, shipboard Linux systems, and institutional servers.

## Proposed architecture

### 1. Data ingestion

This layer would standardize handling of raw LADCP files, CTD exports, station metadata, GPS fixes, and optional shipboard ADCP references. Reliable ingestion is a large part of the operational burden, and cleaning up this stage would remove one of the recurring sources of downstream confusion (Cruise Report P02E, 2022; SourceForge, 2021).

### 2. Core transformations

This layer would expose orientation and geometry explicitly, including tilt, heading, and frame-rotation assumptions. Existing tools already perform these operations, but often in ways that are harder to inspect or reuse in a transparent Python pipeline (Thurnherr, n.d.; University of Hawaiʻi at Mānoa, 2015; apaloczy, 2018).

### 3. Solution engine

This layer would provide multiple solution paths: a shear-based workflow, an inverse workflow, and comparison modes that show where the methods diverge. Preserving both methodological families is important for continuity with existing practice and for scientific trust (Thurnherr, n.d.; JAMSTEC, n.d.).

### 4. QA and diagnostics

This layer would generate the plots and summaries operators actually need, including tilt and heading diagnostics, residual checks, bottom-track diagnostics where relevant, and cast summary reports. Existing documentation makes clear that interpretability is central to reliable LADCP processing, not a cosmetic extra (University of Hawaiʻi at Mānoa, 2015).

### 5. Reproducible deployment

This layer would package the toolkit as both a Python API and a command-line application, with Docker images for repeatable execution. That would support routine reprocessing, batch operation, and deployment on servers without relying on ad hoc desktop setups (GitHub, 2014; SourceForge, 2021).

## Validation plan

The safest path would be staged validation against established workflows rather than an all-at-once rewrite.

1. Benchmark the Python implementation against the LDEO MATLAB workflow on a small set of trusted reference casts (Thurnherr, n.d.; University of Hawaiʻi at Mānoa, 2015).
2. Replicate key diagnostics and final velocity profiles within predefined tolerance bands.
3. Test sensitivity to heading, tilt, timing offsets, and sound-speed assumptions.
4. Compare constrained and unconstrained solutions where bottom track, shipboard ADCP, or GPS references are available (NOAA Atlantic Oceanographic and Meteorological Laboratory, 2024; Thurnherr, 2010).
5. Freeze a validated minimum viable workflow before broadening scope.

This approach treats the legacy system as the scientific reference implementation while allowing the new toolkit to earn confidence through side-by-side reproducibility.

## Deliverables

| Deliverable | Description |
|---|---|
| Python package | Core LADCP processing library with tests and documented interfaces. |
| CLI workflow | Batch and single-cast processing from the command line. |
| Container image | Docker-based runtime for reproducible deployment. |
| QC outputs | Standard plots, summaries, and logs for operator review. |
| Validation dataset | A curated set of casts and expected outputs for regression testing. |
| Documentation | User guide, processing assumptions, troubleshooting notes, and examples. |

## Risks and constraints

The main technical risk is scientific equivalence. Small choices in synchronization, heading treatment, inversion constraints, or ancillary data use can materially change the resulting profile, so validation has to be treated as a central deliverable rather than an afterthought (Thurnherr, n.d.; Thurnherr, 2010).

A second risk is scope creep. It would be easy to let the project broaden from “modern LADCP processing” into a general-purpose ADCP ecosystem, which would dilute effort and delay useful results (SourceForge, 2021; IAHR, n.d.).

A third risk is maintenance load. This proposal is most realistic as a future, phased project with an intentionally narrow first milestone rather than as an attempt to solve every workflow issue at once.

## Recommended first phase

A sensible first phase would be deliberately modest:

- Support one instrument family and one cruise data convention first.
- Reproduce a limited set of LDEO outputs before adding new features (Thurnherr, n.d.; University of Hawaiʻi at Mānoa, 2015).
- Focus on diagnostics, provenance, and repeatability before interface polish.
- Defer GUI development until the scientific core is stable.
- Preserve compatibility with inherited shell and Perl-era conventions where practical, even if the main implementation is Python (Thurnherr, n.d.; Cruise Report P02E, 2022).

## Why this proposal is worth keeping

This project is worth preserving because it addresses a real, technically interesting, and still underserved gap. The current ecosystem demonstrates that the science is mature, but the software stack remains awkwardly split between MATLAB code, shell wrappers, Perl-based utilities, and institution-specific practice (Thurnherr, n.d.; Cruise Report P02E, 2022; GitHub, 2014).

A modern Python implementation would align far better with current research computing practice and with long-term goals around reproducibility, automation, and maintainable scientific software. Even if development is deferred, the proposal captures a meaningful and tractable future contribution.

## References

Except where noted, the following documents and code have been downloaded and copied into the docs/legacy directory

apaloczy. (2018). *ADCPtools: Beam to Earth coordinate transformations and other utilities for ADCP data* [GitHub repository]. GitHub. https://github.com/apaloczy/ADCPtools

British Oceanographic Data Centre. (n.d.). *How to process LADCP data with the LDEO software*. https://www.bodc.ac.uk/data/documents/nodb/pdf/ladcp_ldeo_processing_IX.7_IX.10.pdf

Cruise Report P02E. (2022). *LADCP*. https://cruise-report-2022-p02e.readthedocs.io/en/latest/ladcp.html

GitHub. (2014). *jgrelet/ladcp: Matlab LADCP-processing system based on IFM-GEOMAR/LDEO* [GitHub repository]. https://github.com/jgrelet/ladcp

IAHR. (n.d.). *ADCPtool: An open source ADCP postprocessing framework*. https://www.iahr.org/library/infor?pid=15031

Integrated Marine Observing System. (2025). *Argo floats*. https://imos.org.au/facility/argo-floats

International Pacific Research Center. (n.d.). *Argo deep ocean currents*. https://iprc.soest.hawaii.edu/newsletters/newsletter_sections/iprc_climate_vol6_1/argo_deep_ocean_currents.pdf

JAMSTEC. (n.d.). *5.21 LADCP*. https://www.jamstec.go.jp/cindy/obs/CruiseReport_MR11-07_partial/5.21_LADCP.pdf

NOAA Atlantic Oceanographic and Meteorological Laboratory. (2024, June 6). *Live Science Update: Acoustic measurements determine water column dynamics with LADCP technology*. https://www.aoml.noaa.gov/ladcp-technology-i07/

SourceForge. (2021). *adcptools download*. https://sourceforge.net/projects/adcptools/ - no longer available online

Thurnherr, A. M. (2010). A practical assessment of the errors associated with full-depth LADCP profiles obtained using Teledyne RDI Workhorse acoustic Doppler current profilers. *Journal of Atmospheric and Oceanic Technology, 27*(7), 1215-1227. https://www2.whoi.edu/site/ladder/wp-content/uploads/sites/59/2020/02/thurnherr_2010_JAOT_84904.pdf

Thurnherr, A. M. (n.d.). *Acquisition and processing of LADCP data*. Lamont-Doherty Earth Observatory. https://www.ldeo.columbia.edu/~ant/LADCP.html

University of Hawaiʻi at Mānoa. (2015). *Processing LADCP data with LDEO Matlab code*. https://currents.soest.hawaii.edu/docs/ladcp_doc/Computer/ProcessingPC/ldeo_matlab_processing.html