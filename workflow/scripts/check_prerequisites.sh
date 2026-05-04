#!/bin/bash
# =============================================================================
# check_prerequisites.sh
#
# Verifies that the base system-level tools required by the pipeline are
# available on the current host. Run this once before installing the
# workflow-specific conda environments.
#
# Checks:
#   1. conda            — needed to bootstrap snakemake_env and the auto-created
#                         vcf_annotation / sra_tools envs
#   2. git              — needed to clone the workflow repository
#   3. apptainer OR     — needed to run the containerized tools
#      singularity        (STAR, GATK, fastp, FastQC, MultiQC, bcftools, ...)
#   4. outbound HTTPS   — needed to download reference files and pull
#                         container images on first run
#
# Usage:
#   bash workflow/scripts/check_prerequisites.sh
#
# Exit code:
#   0 — all prerequisites satisfied
#   1 — one or more prerequisites missing
# =============================================================================

set -uo pipefail

missing=0

check_command() {
    local tool="$1"
    local version_flag="${2:---version}"
    if command -v "$tool" >/dev/null 2>&1; then
        local version
        version=$("$tool" "$version_flag" 2>&1 | head -1)
        printf "[OK]      %-12s %s\n" "$tool" "$version"
        return 0
    else
        printf "[MISSING] %-12s not found on PATH\n" "$tool"
        missing=$((missing + 1))
        return 1
    fi
}

echo "Checking base prerequisites for the GATK RNA-seq workflow..."
echo

# 1. Conda
check_command conda

# 2. Git
check_command git

# 3. Apptainer or Singularity (either is acceptable)
if command -v apptainer >/dev/null 2>&1; then
    printf "[OK]      %-12s %s\n" "apptainer" "$(apptainer --version 2>&1 | head -1)"
elif command -v singularity >/dev/null 2>&1; then
    printf "[OK]      %-12s %s\n" "singularity" "$(singularity --version 2>&1 | head -1)"
else
    printf "[MISSING] %-12s neither apptainer nor singularity found on PATH\n" "apptainer"
    missing=$((missing + 1))
fi

# 4. Outbound HTTPS (probe a small, stable target)
if command -v curl >/dev/null 2>&1; then
    if curl --silent --fail --max-time 10 --head https://github.com >/dev/null 2>&1; then
        printf "[OK]      %-12s reachable (https://github.com)\n" "outbound"
    else
        printf "[WARN]    %-12s could not reach https://github.com within 10s\n" "outbound"
        echo "          Reference downloads and container pulls may fail."
        echo "          On clusters that firewall compute nodes, run downloads"
        echo "          from a login or data-transfer node."
    fi
else
    printf "[WARN]    %-12s curl not available; skipping outbound HTTPS probe\n" "outbound"
fi

echo

if [ "$missing" -eq 0 ]; then
    echo "All required prerequisites are available."
    echo "Next step: create the Snakemake environment with"
    echo "    conda env create -f workflow/envs/snakemake.yaml"
    exit 0
else
    echo "Missing prerequisites: $missing"
    echo
    echo "These are system-level dependencies that cannot be installed by an"
    echo "unprivileged user. Two common causes:"
    echo "  - On clusters using environment modules (Lmod, Tcl), the tools are"
    echo "    installed but not visible until you run e.g. 'module load conda'."
    echo "    Check your cluster's documentation for the exact module names."
    echo "  - On clusters where the tool genuinely is not installed, contact"
    echo "    your administrator to request it."
    exit 1
fi
