args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop("usage: Rscript benchmark/esri_tutorial_runner.R <prepared_dir> <mode> <tol> <scenarios> <metrics_out> <scores_out> <gl_ncores|false> <gl_load_balance>")
}

prepared_dir <- normalizePath(args[1], mustWork = TRUE)
mode <- args[2]
tol <- as.numeric(args[3])
scenario_count <- as.integer(args[4])
metrics_out <- args[5]
scores_out <- args[6]
gl_ncores_raw <- tolower(args[7])
gl_load_balance_raw <- tolower(args[8])

if (!(mode %in% c("linear", "leontief"))) {
  stop("mode must be one of: linear, leontief")
}

if (gl_ncores_raw %in% c("false", "0")) {
  gl_ncores <- FALSE
} else {
  gl_ncores <- as.integer(args[7])
  if (is.na(gl_ncores) || gl_ncores <= 0) {
    stop("gl_ncores must be a positive integer or false")
  }
}

if (gl_load_balance_raw %in% c("true", "1", "yes", "on")) {
  gl_load_balance <- TRUE
} else if (gl_load_balance_raw %in% c("false", "0", "no", "off")) {
  gl_load_balance <- FALSE
} else {
  stop("gl_load_balance must be true or false")
}

suppressPackageStartupMessages(library(GLcascade))
suppressPackageStartupMessages(library(fastcascade))
suppressPackageStartupMessages(library(Matrix))

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

edges <- read.table(file.path(prepared_dir, "edges.tsv"), sep = "\t", col.names = c("i", "j", "x"))
industry_ids <- scan(file.path(prepared_dir, "industries.tsv"), what = integer(), quiet = TRUE)
industry_labels <- scan(file.path(prepared_dir, "industry_labels.tsv"), what = character(), quiet = TRUE)
row_weights <- scan(file.path(prepared_dir, "row_weights.tsv"), what = double(), quiet = TRUE)

n <- length(industry_ids)
if (n == 0) {
  stop("prepared industries.tsv is empty")
}
k <- min(scenario_count, n)

W <- Matrix::sparseMatrix(i = edges$i, j = edges$j, x = edges$x, dims = c(n, n))
p <- industry_labels[industry_ids]
ess_value <- if (mode == "linear") 1 else 2
ess_mat_sec <- matrix(ess_value, nrow = length(industry_labels), ncol = length(industry_labels))
rownames(ess_mat_sec) <- industry_labels
colnames(ess_mat_sec) <- industry_labels
psi_mat <- Matrix::Diagonal(n = n)[, seq_len(k), drop = FALSE]
h_weights <- cbind(row_weights)

start_elapsed <- proc.time()[["elapsed"]]
result <- GLcascade::GL_cascade(
  W = W,
  p = p,
  # Match the ESRI / legacy model's within-industry rationing term.
  p_market = p,
  p_sec_impacts = FALSE,
  ess_mat_sec = ess_mat_sec,
  h_weights = h_weights,
  sec_aggr_weights = FALSE,
  psi_mat = psi_mat,
  revenue = FALSE,
  costs = FALSE,
  track_h = FALSE,
  track_sector_impacts = FALSE,
  track_conv = TRUE,
  conv_type = 1,
  eps = tol,
  use_rcpp = TRUE,
  ncores = gl_ncores,
  run_id = paste0("esri_tutorial_", mode),
  load_balance = gl_load_balance
)
elapsed <- proc.time()[["elapsed"]] - start_elapsed

scores <- as.numeric(result$ESRI[, "ESRI_weight_1"])
iters <- as.numeric(result$ESRI_conv[, 3])
metrics <- data.frame(
  key = c("gl_total_s", "gl_score_min", "gl_score_max", "gl_max_iter", "gl_mean_iter", "scenario_count", "gl_ncores", "gl_load_balance"),
  value = c(
    elapsed,
    min(scores),
    max(scores),
    max(iters),
    mean(iters),
    length(scores),
    if (identical(gl_ncores, FALSE)) 0 else gl_ncores,
    as.integer(gl_load_balance)
  )
)

write.table(metrics, metrics_out, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(scores, scores_out, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
