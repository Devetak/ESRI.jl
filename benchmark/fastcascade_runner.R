args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop("usage: Rscript benchmark/fastcascade_runner.R <edges.tsv> <industries.tsv> <essential.tsv> <tol> <scenario_count|all> <metrics.tsv> <scores.tsv>")
}

edges_path <- args[1]
industries_path <- args[2]
essential_path <- args[3]
tol <- as.numeric(args[4])
scenario_count_raw <- args[5]
metrics_path <- args[6]
scores_path <- args[7]

if (!suppressWarnings(requireNamespace("GLcascade", quietly = TRUE))) {
  source_path <- Sys.getenv("GLCASCADE_SOURCE")
  if (!nzchar(source_path) || !file.exists(source_path)) {
    stop("GLcascade is not installed and GLCASCADE_SOURCE does not point to Diem_2022_GLCascade.R")
  }
  source(source_path)
} else {
  suppressPackageStartupMessages(library(GLcascade))
}
suppressPackageStartupMessages(library(fastcascade))
suppressPackageStartupMessages(library(Matrix))

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

edges <- read.table(edges_path, sep = "\t", col.names = c("i", "j", "x"))
industry_ids <- scan(industries_path, what = integer(), quiet = TRUE)
essential_flags <- scan(essential_path, what = integer(), quiet = TRUE)

n <- length(industry_ids)
if (length(essential_flags) == 0 || n == 0) {
  stop("industries and essential flags must be non-empty")
}

if (tolower(scenario_count_raw) %in% c("all", "full")) {
  k <- n
} else {
  k <- min(as.integer(scenario_count_raw), n)
  if (is.na(k) || k <= 0) {
    stop("scenario_count must be a positive integer or 'all'")
  }
}

W <- Matrix::sparseMatrix(i = edges$i, j = edges$j, x = edges$x, dims = c(n, n))
row_weights <- Matrix::rowSums(W)
nindustries <- length(essential_flags)
industry_labels <- as.character(seq_len(nindustries))
p <- industry_labels[industry_ids]

ess_values <- ifelse(essential_flags == 1L, 2L, 1L)
ess_mat_sec <- matrix(rep(ess_values, nindustries), nrow = nindustries, ncol = nindustries)
rownames(ess_mat_sec) <- industry_labels
colnames(ess_mat_sec) <- industry_labels

psi_mat <- Matrix::Diagonal(n = n)[, seq_len(k), drop = FALSE]
h_weights <- cbind(row_weights)

start_elapsed <- proc.time()[["elapsed"]]
result <- GL_cascade(
  W = W,
  p = p,
  p_market = FALSE,
  p_sec_impacts = FALSE,
  ess_mat_sec = ess_mat_sec,
  h_weights = h_weights,
  sec_aggr_weights = FALSE,
  psi_mat = psi_mat,
  revenue = FALSE,
  costs = FALSE,
  track_h = FALSE,
  track_sector_impacts = FALSE,
  track_conv = FALSE,
  conv_type = 1,
  eps = tol,
  use_rcpp = TRUE,
  ncores = FALSE,
  run_id = "fastcascade_subset",
  load_balance = FALSE
)
elapsed <- proc.time()[["elapsed"]] - start_elapsed

scores <- cbind(
  as.numeric(result$ESRI[, "ESRI_weight_1"]),
  as.numeric(result$ESRI[, "ESRI_d_weight_1"]),
  as.numeric(result$ESRI[, "ESRI_u_weight_1"])
)

metrics <- data.frame(
  key = c(
    "diem_total_s",
    "diem_score_min",
    "diem_score_max",
    "diem_downstream_score_min",
    "diem_downstream_score_max",
    "diem_upstream_score_min",
    "diem_upstream_score_max"
  ),
  value = c(
    elapsed,
    min(scores[, 1]),
    max(scores[, 1]),
    min(scores[, 2]),
    max(scores[, 2]),
    min(scores[, 3]),
    max(scores[, 3])
  )
)

write.table(metrics, metrics_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(scores, scores_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
