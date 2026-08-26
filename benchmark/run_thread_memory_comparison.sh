#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  ESRI_IHS_MATRIX=/path/to/EssMatIHS.csv \
    benchmark/run_thread_memory_comparison.sh

Runs the full economy-wide ESRI workload on the same deterministic network:
10,000 firms at 2/4/8 Julia threads or R workers, and 100,000 firms at 8.
The shared network has a truncated power-law out-degree distribution.
Prepared-solve time comes from the existing runners. RAM is sampled peak RSS,
summed over an isolated process group so re-parented R workers are included.

Optional environment variables:
  ESRI_BENCHMARK_CASES      Space-separated firms:workers cases
                            (default: "10000:2 10000:4 10000:8 100000:8")
  ESRI_CLASSIFICATION       ihs, legacy, or linear (default: ihs)
  ESRI_CONVERGENCE_TOL      Shared convergence tolerance (default: 1e-2)
  ESRI_OUTPUT_DIR           Fresh output directory
  ESRI_NETWORK_DIR          Directory for shared generated networks
                            (default: output directory)
  ESRI_RSS_SAMPLE_SECONDS   RSS sampling interval (default: 0.1)
  ESRI_UNRESTRICT_BLAS      1 to leave BLAS/OpenMP thread counts unrestricted
                            (default: restricted to one)
EOF
}

process_group_rss_kib() {
  LC_ALL=C ps -axo pgid=,rss= | awk -v group="$1" '
    $1 == group { total += $2 }
    END {
      printf "%.0f\n", total
    }
  '
}

monitor_peak_rss() {
  group_id=$1
  root_pid=$2
  output_file=$3
  peak=0
  while kill -0 "$root_pid" 2>/dev/null; do
    if ! current=$(process_group_rss_kib "$group_id"); then
      echo "RSS sampling failed" >&2
      terminate_workload
      return 1
    fi
    if [ "$current" -gt "$peak" ]; then
      peak=$current
    fi
    sleep "$sample_seconds"
  done
  printf '%s\n' "$peak" > "$output_file"
}

measured_pid=
measured_pgid=
monitor_pid=
terminate_workload() {
  if [ -n "$measured_pgid" ]; then
    kill -TERM "-$measured_pgid" 2>/dev/null || :
    sleep 1
    kill -KILL "-$measured_pgid" 2>/dev/null || :
  elif [ -n "$measured_pid" ]; then
    kill "$measured_pid" 2>/dev/null || :
  fi
}

stop_measurement() {
  trap - HUP INT TERM
  [ -z "$monitor_pid" ] || kill "$monitor_pid" 2>/dev/null || :
  terminate_workload
  [ -z "$measured_pid" ] || wait "$measured_pid" 2>/dev/null || :
  [ -z "$monitor_pid" ] || wait "$monitor_pid" 2>/dev/null || :
  exit 130
}
trap stop_measurement HUP INT TERM

measure_peak_rss() {
  rss_file=$1
  shift
  set -m
  "$@" &
  measured_pid=$!
  set +m
  candidate_pgid=$(ps -o pgid= -p "$measured_pid" | awk 'NR == 1 { print $1 }') || candidate_pgid=
  driver_pgid=$(ps -o pgid= -p "$$" | awk 'NR == 1 { print $1 }') || driver_pgid=
  if [ -z "$candidate_pgid" ] || [ -z "$driver_pgid" ] || [ "$candidate_pgid" = "$driver_pgid" ]; then
    kill "$measured_pid" 2>/dev/null || :
    wait "$measured_pid" 2>/dev/null || :
    measured_pid=
    echo "Could not isolate workload process group" >&2
    return 1
  fi
  measured_pgid=$candidate_pgid
  monitor_peak_rss "$measured_pgid" "$measured_pid" "$rss_file" &
  monitor_pid=$!

  if wait "$measured_pid"; then
    status=0
  else
    status=$?
  fi
  if wait "$monitor_pid"; then
    monitor_status=0
  else
    monitor_status=$?
  fi
  if [ "$status" -ne 0 ] || [ "$monitor_status" -ne 0 ]; then
    terminate_workload
  fi
  [ "$status" -ne 0 ] || status=$monitor_status
  measured_pid=
  monitor_pid=
  measured_pgid=
  [ -s "$rss_file" ] || {
    echo "RSS monitor produced no result" >&2
    return 1
  }
  return "$status"
}

thread_env() {
  if [ "${ESRI_UNRESTRICT_BLAS:-0}" = 1 ]; then
    env -u OMP_NUM_THREADS -u OPENBLAS_NUM_THREADS -u MKL_NUM_THREADS \
      -u VECLIB_MAXIMUM_THREADS "$@"
  else
    env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
      VECLIB_MAXIMUM_THREADS=1 "$@"
  fi
}

sample_seconds=${ESRI_RSS_SAMPLE_SECONDS:-0.1}
awk -v value="$sample_seconds" 'BEGIN {
  if (value !~ /^[0-9]*[.]?[0-9]+$/ || value + 0 <= 0) exit 1
}' || { echo "ESRI_RSS_SAMPLE_SECONDS must be a positive number" >&2; exit 2; }
case "${1-}" in
  --help|-h)
    usage
    exit 0
    ;;
  --self-test)
    rss_file=$(mktemp "${TMPDIR:-/tmp}/esri-rss.XXXXXX")
    measure_peak_rss "$rss_file" sh -c 'sleep 0.5 & wait'
    rss=$(cat "$rss_file")
    rm -f "$rss_file"
    case "$rss" in
      ''|*[!0-9]*) echo "RSS self-test returned invalid output: $rss" >&2; exit 1 ;;
    esac
    [ "$rss" -gt 0 ] || { echo "RSS self-test measured zero KiB" >&2; exit 1; }
    [ "$rss" -lt 1048576 ] || { echo "RSS self-test included unrelated processes: $rss KiB" >&2; exit 1; }
    echo "RSS self-test passed: ${rss} KiB"
    exit 0
    ;;
  '') ;;
  *) usage >&2; exit 2 ;;
esac

: "${ESRI_IHS_MATRIX:?Set ESRI_IHS_MATRIX to the labeled EssMatIHS.csv file}"
[ -f "$ESRI_IHS_MATRIX" ] || { echo "File not found: $ESRI_IHS_MATRIX" >&2; exit 1; }
for required in awk julia ps Rscript; do
  command -v "$required" >/dev/null 2>&1 || { echo "Missing command: $required" >&2; exit 1; }
done

script_dir=$(CDPATH= cd "$(dirname "$0")" && pwd)
repo_dir=$(dirname "$script_dir")
cd "$repo_dir"

benchmark_cases=${ESRI_BENCHMARK_CASES:-"10000:2 10000:4 10000:8 100000:8"}
classification=${ESRI_CLASSIFICATION:-ihs}
run_id=$(date -u +%Y%m%dT%H%M%SZ)
output_dir=${ESRI_OUTPUT_DIR:-"results/thread_memory_comparison/$run_id"}
[ ! -e "$output_dir" ] || { echo "Output already exists: $output_dir" >&2; exit 1; }
mkdir -p "$output_dir"
network_dir=${ESRI_NETWORK_DIR:-$output_dir}
mkdir -p "$network_dir"

cat > "$output_dir/measurement.txt" <<EOF
timing=prepared solve reported by the existing Julia and FastCascade runners
memory=sampled peak aggregate RSS of the isolated workload process group
memory_scope=full invocation, including setup, warmup, solve, and validation
memory_caveat=includes re-parented R workers; aggregate RSS can count shared pages once per process
blas_openmp_restricted=$([ "${ESRI_UNRESTRICT_BLAS:-0}" = 1 ] && echo no || echo yes)
network=directed truncated power-law out-degree, mean 8, alpha 2.3, cap 128, seed 42, no self-links
network_dir=$network_dir
convergence_tolerance=${ESRI_CONVERGENCE_TOL:-1e-2}
rss_units=KiB from ps; results.csv converts to MiB
rss_sample_seconds=$sample_seconds
system=$(uname -srm)
benchmark_cases=$benchmark_cases
EOF

results_tmp="$output_dir/results.csv.tmp"
echo "classification,firms,workers,julia_total_s,julia_s_per_firm,cpp_total_s,cpp_s_per_firm,speedup,julia_peak_group_rss_mib,fastcascade_r_peak_group_rss_mib,fastcascade_r_over_julia_peak_rss,julia_total_esri,cpp_total_esri,max_abs_score_difference,max_difference_firm,julia_score_at_max_difference,cpp_score_at_max_difference,mean_abs_score_difference,abs_total_esri_difference,comparison_tolerance,within_tolerance" > "$results_tmp"

for benchmark_case in $benchmark_cases; do
  case "$benchmark_case" in
    *:*) ;;
    *) echo "Invalid ESRI_BENCHMARK_CASES entry: $benchmark_case" >&2; exit 2 ;;
  esac
  firms=${benchmark_case%:*}
  workers=${benchmark_case#*:}
  case "$firms,$workers" in
    *[!0-9,]*|,*|*,|0,*|*,0) echo "Invalid ESRI_BENCHMARK_CASES entry: $benchmark_case" >&2; exit 2 ;;
  esac
  network_file="$network_dir/powerlaw_network_${firms}.csv"
  if [ ! -f "$network_file" ]; then
    julia --project=. benchmark/generate_powerlaw_network.jl "$firms" "$network_file"
  fi
  worker_dir="$output_dir/workers_$workers/firms_$firms"
  mkdir -p "$worker_dir"
  echo "Running $firms firms with $workers workers: FastCascade/R"
  measure_peak_rss "$worker_dir/fastcascade_r_peak_group_rss_kib.txt" thread_env env \
    ESRI_BENCHMARK_FIRMS="$firms" \
    ESRI_BENCHMARK_WORKERS="$workers" \
    ESRI_CLASSIFICATION="$classification" \
    ESRI_NETWORK_FILE="$network_file" \
    ESRI_OUTPUT_DIR="$worker_dir" \
    ESRI_UNRESTRICT_BLAS="${ESRI_UNRESTRICT_BLAS:-0}" \
    Rscript benchmark/full_esri_fastcascade.R

  echo "Running $firms firms with $workers workers: Julia"
  measure_peak_rss "$worker_dir/julia_peak_group_rss_kib.txt" thread_env env \
    ESRI_BENCHMARK_FIRMS="$firms" \
    ESRI_BENCHMARK_WORKERS="$workers" \
    ESRI_CLASSIFICATION="$classification" \
    ESRI_NETWORK_FILE="$network_file" \
    ESRI_OUTPUT_DIR="$worker_dir" \
    JULIA_NUM_THREADS="$workers" \
    ESRI_UNRESTRICT_BLAS="${ESRI_UNRESTRICT_BLAS:-0}" \
    julia --project=. benchmark/full_esri_julia.jl

  comparison="$worker_dir/${classification}_${firms}/comparison.csv"
  [ -f "$comparison" ] || { echo "Missing comparison: $comparison" >&2; exit 1; }
  julia_rss=$(cat "$worker_dir/julia_peak_group_rss_kib.txt")
  fastcascade_r_rss=$(cat "$worker_dir/fastcascade_r_peak_group_rss_kib.txt")
  case "$julia_rss,$fastcascade_r_rss" in
    *[!0-9,]*|,*|*,|0,*|*,0) echo "Invalid RSS measurement for $workers workers" >&2; exit 1 ;;
  esac
  tail -n 1 "$comparison" | awk -F, \
    -v expected_classification="$classification" \
    -v expected_firms="$firms" \
    -v expected_workers="$workers" \
    -v julia_rss="$julia_rss" \
    -v fastcascade_r_rss="$fastcascade_r_rss" '
    {
      if (NF != 16 || $1 != expected_classification || $2 != expected_firms || \
          $3 != expected_workers) exit 1
      ratio = julia_rss == 0 ? 0 : fastcascade_r_rss / julia_rss
      printf "%s,%s,%s,%s,%.12f,%s,%.12f,%s,%.3f,%.3f,%.6f,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", \
        $1, $2, $3, $4, $4 / $2, $5, $5 / $2, $6, julia_rss / 1024, fastcascade_r_rss / 1024, ratio, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16
      valid = 1
    }
    END { if (!valid) exit 1 }
  ' >> "$results_tmp"
done

mv "$results_tmp" "$output_dir/results.csv"
echo "Wrote $output_dir/results.csv"
