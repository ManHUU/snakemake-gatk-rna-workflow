#!/bin/bash
# =============================================================================
# run_smoke_test.sh  (reviewer comments C-07 / C-08 / G-02)
#
# One documented command to verify a working installation on the MINIATURE
# fixture: it builds the fixture if missing, runs the real workflow + containers
# through run.sh with config/smoke/config.yaml, then checks the expected
# outputs. This is a fast installation check (design target ~15-30 min of
# execution on the documented HPC allocation) - NOT the full-genome tutorial.
#
# Usage:
#   bash workflow/scripts/run_smoke_test.sh              # build (if needed) + run + check
#   bash workflow/scripts/run_smoke_test.sh --dry-run    # plan only
#
# Report separately (do not fold into the execution time): scheduler queue
# time, first-time container pulls, and fixture-download/build time. Record the
# measured wall time, peak RAM, disk use and allocation in the manuscript
# (currently [MEASURE] placeholders) from the first successful run.
# =============================================================================

set -euo pipefail
cd "$(cd "$(dirname "$0")/../.." && pwd)"

SMOKE_CFG="config/smoke/config.yaml"
SMOKE_OUT="results_smoke"
SMOKE_FASTA="resources/smoke/smoke.fa"

if [[ ! -e "$SMOKE_FASTA" || ! -d data/smoke ]]; then
    echo "Smoke fixture not found - building it first ..."
    bash workflow/scripts/make_smoke_fixture.sh
fi

echo "Running smoke test (config: $SMOKE_CFG) ..."
start=$(date +%s)
bash run.sh --configfile "$SMOKE_CFG" "$@"
rc=$?
end=$(date +%s)
echo "run.sh returned $rc; wall time: $(( end - start )) s (execution only; excludes queue/pull/build)."

# Skip the output check on a dry run.
for a in "$@"; do [[ "$a" == "--dry-run" || "$a" == "-n" ]] && exit 0; done

echo "Checking smoke-test outputs ..."
bash workflow/scripts/check_smoke_outputs.sh "$SMOKE_OUT" config/smoke/expected_checkpoints.tsv
