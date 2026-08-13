#!/usr/bin/env bash
set -euo pipefail

: "${ESRI_IHS_MATRIX:?Set ESRI_IHS_MATRIX to the labeled EssMatIHS.csv file}"

script_dir=$(CDPATH= cd "$(dirname "$0")" && pwd)
repo_dir=$(dirname "$script_dir")
cd "$repo_dir"

run_id=$(date -u +%Y%m%dT%H%M%SZ)
output_root=${ESRI_OVERNIGHT_OUTPUT_DIR:-"results/overnight_powerlaw/$run_id"}
[ ! -e "$output_root" ] || { echo "Output already exists: $output_root" >&2; exit 1; }
mkdir -p "$output_root/networks"

cat > "$output_root/experiment_matrix.txt" <<'EOF'
network=directed truncated power law, mean out-degree 8, alpha 2.3, cap 128, seed 42
workers=1 4
10000_firms=ihs legacy linear
100000_firms=ihs
comparison=Julia versus FastCascade/R, one score per firm, tolerance 1e-2
EOF

for classification in ihs legacy linear; do
  cases="10000:1 10000:4"
  [ "$classification" != ihs ] || cases="$cases 100000:1 100000:4"
  echo "Starting $classification cases: $cases"
  ESRI_BENCHMARK_CASES="$cases" \
    ESRI_CLASSIFICATION="$classification" \
    ESRI_NETWORK_DIR="$output_root/networks" \
    ESRI_OUTPUT_DIR="$output_root/$classification" \
    benchmark/run_thread_memory_comparison.sh
done

results_tmp="$output_root/results.csv.tmp"
head -n 1 "$output_root/ihs/results.csv" > "$results_tmp"
for classification in ihs legacy linear; do
  tail -n +2 "$output_root/$classification/results.csv" >> "$results_tmp"
done
mv "$results_tmp" "$output_root/results.csv"
echo "Wrote $output_root/results.csv"
