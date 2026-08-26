#!/bin/sh
set -eu

: "${ESRI_IHS_MATRIX:?Set ESRI_IHS_MATRIX to the labeled EssMatIHS.csv file}"

workers="${ESRI_BENCHMARK_WORKERS:-8}"
sizes="${ESRI_FIRM_SIZES:-10000 20000 50000 100000}"
classifications="${ESRI_CLASSIFICATIONS:-ihs legacy linear}"
output_dir="${ESRI_OUTPUT_DIR:-results/full_esri_matrix_comparison}"
mkdir -p "$output_dir"

for firms in $sizes; do
  network_file="$output_dir/powerlaw_network_${firms}.csv"
  if [ ! -f "$network_file" ]; then
    julia --project=. benchmark/generate_powerlaw_network.jl "$firms" "$network_file"
  fi
  for classification in $classifications; do
    echo "Running $classification classification with $firms firms"
    ESRI_BENCHMARK_FIRMS="$firms" \
      ESRI_BENCHMARK_WORKERS="$workers" \
      ESRI_CLASSIFICATION="$classification" \
      ESRI_NETWORK_FILE="$network_file" \
      ESRI_OUTPUT_DIR="$output_dir" \
      OMP_NUM_THREADS=1 \
      OPENBLAS_NUM_THREADS=1 \
      MKL_NUM_THREADS=1 \
      VECLIB_MAXIMUM_THREADS=1 \
      Rscript benchmark/full_esri_fastcascade.R

    ESRI_BENCHMARK_FIRMS="$firms" \
      ESRI_BENCHMARK_WORKERS="$workers" \
      ESRI_CLASSIFICATION="$classification" \
      ESRI_NETWORK_FILE="$network_file" \
      ESRI_OUTPUT_DIR="$output_dir" \
      JULIA_NUM_THREADS="$workers" \
      OPENBLAS_NUM_THREADS=1 \
      OMP_NUM_THREADS=1 \
      julia --project=. benchmark/full_esri_julia.jl
  done
done

results_tmp="$output_dir/results.csv.tmp"
echo "classification,firms,workers,julia_total_s,cpp_total_s,speedup,julia_total_esri,cpp_total_esri,max_abs_score_difference,max_difference_firm,julia_score_at_max_difference,cpp_score_at_max_difference,mean_abs_score_difference,abs_total_esri_difference,comparison_tolerance,within_tolerance" > "$results_tmp"
for firms in $sizes; do
  for classification in $classifications; do
    tail -n 1 "$output_dir/${classification}_${firms}/comparison.csv" >> "$results_tmp"
  done
done
mv "$results_tmp" "$output_dir/results.csv"
echo "Wrote $output_dir/results.csv"
