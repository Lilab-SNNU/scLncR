#!/usr/bin/env bash
# scLncR raw FASTQ QC wrapper
#
# Purpose:
#   Run the optional scLncR raw FASTQ quality control module from the shell.
#
# Input:
#   -c, --config PATH   YAML config file, normally R/confings/config_QC.yaml
#
# Output:
#   The configured output_dir containing fastqc/, multiqc/, logs/, and
#   qc_summary.md.
#
# Important:
#   This wrapper only runs FastQC and MultiQC through the scLncR qc CLI.
#   It does not trim, filter, or modify original FASTQ files. With dry_run=true,
#   commands and FASTQ file lists are written but no QC tools are executed.

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: bash scripts/run_fastq_qc.sh -c <config.yaml>

Optional raw FASTQ QC for scLncR. Runs FastQC and MultiQC only.

Options:
  -c, --config FILE   Path to QC YAML config
  -h, --help          Show this help message
EOF
}

CONFIG=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -c|--config)
      CONFIG="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ -z "$CONFIG" ]]; then
  echo "Missing required -c/--config argument." >&2
  usage >&2
  exit 1
fi

SCRIPT_DIR="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd -P "$SCRIPT_DIR/.." && pwd)"

exec bash "$REPO_DIR/scLncR" qc -c "$CONFIG"
