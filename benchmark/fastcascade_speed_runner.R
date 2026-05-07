args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9) {
  stop("usage: Rscript benchmark/fastcascade_speed_runner.R <prepared_dir> <tol> <scenario_count|all> <p_market_mode> <gl_ncores|false> <gl_load_balance> <scenario_batch_size|all> <metrics.tsv> <scores.tsv>")
}

prepared_dir <- normalizePath(args[1], mustWork = TRUE)
tol <- as.numeric(args[2])
scenario_count_raw <- tolower(args[3])
p_market_mode <- tolower(args[4])
gl_ncores_raw <- tolower(args[5])
gl_load_balance_raw <- tolower(args[6])
batch_size_raw <- tolower(args[7])
metrics_out <- args[8]
scores_out <- args[9]

suppressPackageStartupMessages(library(GLcascade))
suppressPackageStartupMessages(library(fastcascade))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(parallel))

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

if (gl_ncores_raw %in% c("false", "0")) {
  gl_ncores <- FALSE
} else {
  gl_ncores <- as.integer(gl_ncores_raw)
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

edges <- read.table(file.path(prepared_dir, "edges.tsv"), sep = "\t", col.names = c("i", "j", "x"))
industry_ids <- scan(file.path(prepared_dir, "industries.tsv"), what = integer(), quiet = TRUE)
industry_labels <- scan(file.path(prepared_dir, "industry_labels.tsv"), what = character(), quiet = TRUE)
row_weights <- scan(file.path(prepared_dir, "row_weights.tsv"), what = double(), quiet = TRUE)
essential_flags <- scan(file.path(prepared_dir, "essential.tsv"), what = integer(), quiet = TRUE)

n <- length(industry_ids)
if (n == 0 || length(essential_flags) == 0) {
  stop("prepared_dir is missing industries or essential flags")
}

if (scenario_count_raw %in% c("all", "full")) {
  scenario_count <- n
} else {
  scenario_count <- min(as.integer(scenario_count_raw), n)
  if (is.na(scenario_count) || scenario_count <= 0) {
    stop("scenario_count must be a positive integer or 'all'")
  }
}

if (batch_size_raw %in% c("all", "full")) {
  batch_size <- scenario_count
} else {
  batch_size <- as.integer(batch_size_raw)
  if (is.na(batch_size) || batch_size <= 0) {
    stop("scenario_batch_size must be a positive integer or 'all'")
  }
  batch_size <- min(batch_size, scenario_count)
}

W <- Matrix::sparseMatrix(i = edges$i, j = edges$j, x = edges$x, dims = c(n, n))
p <- industry_labels[industry_ids]
nindustries <- length(essential_flags)
ess_values <- ifelse(essential_flags == 1L, 2L, 1L)
ess_mat_sec <- matrix(rep(ess_values, nindustries), nrow = nindustries, ncol = nindustries)
rownames(ess_mat_sec) <- industry_labels
colnames(ess_mat_sec) <- industry_labels
h_weights <- cbind(row_weights)
sec_aggr_weights <- row_weights
sec_aggr_weights_up <- row_weights

if (p_market_mode %in% c("p", "industry", "sector")) {
  p_market <- p
} else if (p_market_mode %in% c("false", "none", "off")) {
  p_market <- FALSE
} else {
  stop("p_market_mode must be one of: false, p")
}

sectors <- sort(unique(p))
m <- length(sectors)
sectors_cons <- seq_along(sectors)
names(sectors_cons) <- sectors
p_cons <- unname(sectors_cons[as.character(p)])
psup <- Matrix::sparseMatrix(i = seq_len(n), j = p_cons, dims = c(n, m))
pi_abs <- Matrix::crossprod(W, psup)

if (length(p_market) <= 1) {
  substitution <- FALSE
  p_market_cons <- rep(1, n)
  psup_market <- Matrix::Matrix(0, ncol = 1, nrow = n)
} else {
  substitution <- TRUE
  sectors_market <- sort(unique(p_market))
  sectors_market_cons <- seq_along(sectors_market)
  names(sectors_market_cons) <- sectors_market
  p_market_cons <- unname(sectors_market_cons[as.character(p_market)])
  m_market <- length(unique(p_market))
  psup_market <- Matrix::sparseMatrix(i = seq_len(n), j = p_market_cons, dims = c(n, m_market))
}

uni_ess_mat_sec <- unique(as.matrix(ess_mat_sec), MARGIN = 2)
sector_cons_uni_ess_mat <- numeric(length(sectors_cons))
for (k in seq_len(dim(uni_ess_mat_sec)[2])) {
  ess_duplicates <- apply(ess_mat_sec, MARGIN = 2, function(x) identical(uni_ess_mat_sec[, k], x))
  sector_cons_uni_ess_mat[ess_duplicates] <- k
}
names(sector_cons_uni_ess_mat) <- sectors_cons
p_cons_uni_ess_mat <- sector_cons_uni_ess_mat[as.character(p_cons)]

out_str <- Matrix::rowSums(W)
in_str <- Matrix::colSums(W)
setup_elapsed <- proc.time()[["elapsed"]]
Lambda_d <- GLcascade::Lambda_d_calc_r(
  W = W,
  p = p,
  p_cons = p_cons,
  ess_mat_sec = ess_mat_sec,
  pi_abs = pi_abs,
  in_str = in_str,
  n = n
)
Lambda_u <- W / (((out_str == 0) * 1) + ((out_str > 0) * out_str))
setup_elapsed <- proc.time()[["elapsed"]] - setup_elapsed

if (substitution) {
  p_market_cons_cpp <- p_market_cons - 1
  psup_market_cpp <- Matrix::Matrix(psup_market * 1, sparse = TRUE)
} else {
  p_market_cons_cpp <- rep(0, n)
  psup_market_cpp <- Matrix::Matrix(psup_market * 1, sparse = TRUE)
}

context <- list(
  Lambda_d = Matrix::Matrix(Lambda_d, sparse = TRUE),
  Lambda_u = Matrix::Matrix(Lambda_u, sparse = TRUE),
  psup = Matrix::Matrix(psup * 1, sparse = TRUE),
  h_weights = h_weights,
  eps = tol,
  h_mat_round = FALSE,
  p_market_cons = p_market_cons_cpp,
  out_str = out_str,
  psup_market = psup_market_cpp,
  track_h = FALSE,
  track_sector_impacts = FALSE,
  track_conv = FALSE,
  conv_type = 1,
  psup_sec_impact = Matrix::Matrix(0, ncol = 1, nrow = n),
  sec_aggr_weights = sec_aggr_weights,
  sec_aggr_weights_up = sec_aggr_weights_up,
  sectors_cons = sectors_cons,
  p_cons_uni_ess_mat = p_cons_uni_ess_mat - 1,
  uni_ess_mat_sec = uni_ess_mat_sec,
  n = n,
  m = m,
  substitution = substitution
)

run_cpp_direct <- function(psi_mat_batch) {
  fastcascade::GL_cascade_dynamics_cpp(
    Lambda_d = context$Lambda_d,
    Lambda_u = context$Lambda_u,
    psi_mat = psi_mat_batch,
    psi_mat_up = psi_mat_batch,
    n = context$n,
    m = context$m,
    substitution = context$substitution,
    psup = context$psup,
    h_weights = context$h_weights,
    eps = context$eps,
    h_mat_round = context$h_mat_round,
    p_market_cons = context$p_market_cons,
    out_str = context$out_str,
    psup_market = context$psup_market,
    track_h = context$track_h,
    track_sector_impacts = context$track_sector_impacts,
    track_conv = context$track_conv,
    conv_type = context$conv_type,
    psup_sec_impact = context$psup_sec_impact,
    sec_aggr_weights = context$sec_aggr_weights,
    sec_aggr_weights_up = context$sec_aggr_weights_up,
    sectors_cons = context$sectors_cons,
    p_cons_uni_ess_mat = context$p_cons_uni_ess_mat,
    uni_ess_mat_sec = context$uni_ess_mat_sec
  )
}

cluster <- NULL
worker_names <- c(
  "Lambda_d", "Lambda_u", "n", "m", "substitution", "psup", "h_weights", "eps",
  "h_mat_round", "p_market_cons", "out_str", "psup_market", "track_h",
  "track_sector_impacts", "track_conv", "conv_type", "psup_sec_impact",
  "sec_aggr_weights", "sec_aggr_weights_up", "sectors_cons",
  "p_cons_uni_ess_mat", "uni_ess_mat_sec"
)

if (!identical(gl_ncores, FALSE)) {
  cluster <- parallel::makeCluster(gl_ncores, type = "PSOCK")
  parallel::clusterEvalQ(cluster, list(library(Matrix), library(fastcascade), library(Rcpp), library(RcppArmadillo)))
  parallel::clusterExport(
    cluster,
    worker_names,
    envir = list2env(context, parent = emptyenv())
  )
}

on.exit({
  if (!is.null(cluster)) {
    parallel::stopCluster(cluster)
  }
}, add = TRUE)

run_cpp_parallel <- function(psi_mat_batch) {
  psi_splits <- GLcascade::split_mat(psi_mat_batch, gl_ncores, gl_load_balance)
  psi_pairs <- lapply(psi_splits, function(split) list(split, split))

  worker_fn <- function(psi_pair) {
    fastcascade::GL_cascade_dynamics_cpp(
      Lambda_d = Lambda_d,
      Lambda_u = Lambda_u,
      psi_mat = psi_pair[[1]],
      psi_mat_up = psi_pair[[2]],
      n = n,
      m = m,
      substitution = substitution,
      psup = psup,
      h_weights = h_weights,
      eps = eps,
      h_mat_round = h_mat_round,
      p_market_cons = p_market_cons,
      out_str = out_str,
      psup_market = psup_market,
      track_h = track_h,
      track_sector_impacts = track_sector_impacts,
      track_conv = track_conv,
      conv_type = conv_type,
      psup_sec_impact = psup_sec_impact,
      sec_aggr_weights = sec_aggr_weights,
      sec_aggr_weights_up = sec_aggr_weights_up,
      sectors_cons = sectors_cons,
      p_cons_uni_ess_mat = p_cons_uni_ess_mat,
      uni_ess_mat_sec = uni_ess_mat_sec
    )
  }

  if (gl_load_balance) {
    parts <- parallel::parLapplyLB(cluster, psi_pairs, worker_fn)
  } else {
    parts <- parallel::parLapply(cluster, psi_pairs, worker_fn)
  }
  do.call(rbind, parts)
}

make_psi_batch <- function(indices) {
  Matrix::sparseMatrix(
    i = indices,
    j = seq_along(indices),
    x = 1,
    dims = c(n, length(indices))
  )
}

batch_starts <- seq.int(1, scenario_count, by = batch_size)
score_chunks <- vector("list", length(batch_starts))
batch_count <- length(batch_starts)

solve_elapsed <- proc.time()[["elapsed"]]
for (i in seq_along(batch_starts)) {
  start_idx <- batch_starts[i]
  stop_idx <- min(start_idx + batch_size - 1, scenario_count)
  indices <- start_idx:stop_idx
  psi_batch <- make_psi_batch(indices)
  output_mat <- if (identical(gl_ncores, FALSE)) run_cpp_direct(psi_batch) else run_cpp_parallel(psi_batch)
  score_chunks[[i]] <- as.numeric(output_mat[, 1])
}
solve_elapsed <- proc.time()[["elapsed"]] - solve_elapsed

scores <- unlist(score_chunks, use.names = FALSE)

metrics <- data.frame(
  key = c(
    "fast_setup_s",
    "fast_total_s",
    "fast_score_min",
    "fast_score_max",
    "scenario_count",
    "batch_count",
    "gl_ncores",
    "gl_load_balance",
    "p_market_used"
  ),
  value = c(
    setup_elapsed,
    solve_elapsed,
    min(scores),
    max(scores),
    length(scores),
    batch_count,
    if (identical(gl_ncores, FALSE)) 0 else gl_ncores,
    as.integer(gl_load_balance),
    as.integer(!identical(p_market, FALSE))
  )
)

write.table(metrics, metrics_out, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(scores, scores_out, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
