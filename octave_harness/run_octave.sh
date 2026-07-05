#!/usr/bin/env bash
# Run octave-cli inside the gnuoctave/octave:9.2.0 container, with the LADCP
# repo root mounted at /work. Works around MSYS/Git-Bash path mangling on
# Windows hosts (MSYS_NO_PATHCONV=1) -- see CONTINUATION_PLAN.md "Getting
# Octave" section.
#
# Usage:
#   ./run_octave.sh --eval "disp(1+1)"
#   ./run_octave.sh /work/octave_harness/some_script.m
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."
REPO_ROOT="$(pwd -W 2>/dev/null || pwd)"

MSYS_NO_PATHCONV=1 docker run --rm \
  -v "${REPO_ROOT}:/work" \
  -w /work \
  -e OCTAVE_HARNESS_BEGIN_STEP \
  docker.io/gnuoctave/octave:9.2.0 \
  octave-cli --no-gui "$@"
