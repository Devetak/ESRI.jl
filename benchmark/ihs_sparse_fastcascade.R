library(GLcascade)
library(fastcascade)
library(Matrix)

matrix_path <- Sys.getenv("ESRI_IHS_MATRIX")
if (!file.exists(matrix_path)) stop("Set ESRI_IHS_MATRIX to the labeled EssMatIHS.csv file")

d <- read.csv(matrix_path, row.names = 1, check.names = FALSE, colClasses = "character")
E <- apply(d, 2L, as.integer)
rownames(E) <- rownames(d)
colnames(E) <- colnames(d)
stopifnot(nrow(E) == ncol(E), identical(rownames(E), colnames(E)), all(E %in% 0:2))

m <- nrow(E)
n <- as.integer(Sys.getenv("ESRI_BENCHMARK_FIRMS", as.character(2L * m)))
if (n < m) stop("ESRI_BENCHMARK_FIRMS must be at least the number of industries")
firms_per_industry <- n %/% m
extra_firms <- n %% m
p_cons <- c(rep(seq_len(m), each = firms_per_industry), seq_len(extra_firms))
p <- rownames(E)[p_cons]
closure_industries <- as.integer(round(seq(1L, m, length.out = 16L)))
closure_firms <- as.integer(
  (closure_industries - 1L) * firms_per_industry + 1L
)
stopifnot(length(unique(closure_firms)) == 16L, all(p_cons[closure_firms] == closure_industries))
offsets <- c(1L, 37L, 91L, 173L, 311L, 509L, 701L, 997L)
supplier <- rep(seq_len(n), each = length(offsets))
customer <- unlist(lapply(seq_len(n), function(i) ((i + offsets - 1L) %% n) + 1L), use.names = FALSE)
W <- sparseMatrix(
  i = supplier,
  j = customer,
  x = as.numeric(1 + ((13 * supplier + 7 * customer) %% 17)),
  dims = c(n, n)
)

# This is GLcascade's sparse preprocessing, performed once outside the timed kernel.
psup <- sparseMatrix(i = seq_len(n), j = p_cons, x = 1, dims = c(n, m))
in_str <- colSums(W)
out_str <- rowSums(W)
Lambda_d <- Matrix::drop0(Matrix(GLcascade::Lambda_d_calc_r(W, p, p_cons, E, crossprod(W, psup), in_str, n), sparse = TRUE))
Lambda_u <- Matrix(W / ifelse(out_str == 0, 1, out_str), sparse = TRUE)
unique_E <- unique(E, MARGIN = 2L)
profile <- vapply(seq_len(m), function(j) {
  which(vapply(seq_len(ncol(unique_E)), function(k) identical(E[, j], unique_E[, k]), logical(1)))[1]
}, integer(1))
psi <- sparseMatrix(i = closure_firms, j = seq_len(16L), x = 1, dims = c(n, 16L))
zero_sector <- sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = c(n, 1L))

run_cpp <- function() fastcascade::GL_cascade_dynamics_cpp(
  Lambda_d, Lambda_u, psi, psi, n, m, TRUE, psup, cbind(out_str), 1e-2, FALSE,
  as.integer(p_cons - 1L), out_str, psup, FALSE, FALSE, FALSE, 1,
  zero_sector, out_str, out_str, as.integer(seq_len(m)),
  as.integer(profile[p_cons] - 1L), unique_E
)

cpp_scores <- c(
  0.87837245554794330, 0.83928173885048851, 0.84026525189212131, 0.010491619358180508,
  0.83938875684430192, 0.83919913726570627, 0.83871134521356561, 0.83862604469123081,
  0.83842687824673123, 0.14782516685618610, 0.013224577635286240, 0.0064115243633385468,
  0.0046717737458951943, 0.83930463691132073, 0.011700530446422091, 0.0086842594773570356
)

invisible(run_cpp())
scores <- as.numeric(run_cpp()[, 1L])
reference_path <- Sys.getenv("ESRI_REFERENCE_SCORES", "")
if (nzchar(reference_path)) {
  write.table(scores, file = reference_path, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
if (n == 2L * m) {
  cat("cpp_scores=", paste(formatC(scores, digits = 17, format = "g"), collapse = ","), "\n", sep = "")
}
samples <- vapply(
  seq_len(as.integer(Sys.getenv("ESRI_BENCHMARK_SAMPLES", "5"))),
  function(i) system.time(invisible(run_cpp()))[["elapsed"]],
  numeric(1)
)
cat("firms=", n, " industries=", m, " edges=", length(W@x), " closures=16 closure_firms=", paste(closure_firms, collapse = ","), "\n", sep = "")
cat("downstream_nnz=", length(Lambda_d@x), "\n", sep = "")
cat("cpp_prepared_median_s=", median(samples), "\n", sep = "")
if (n == 2L * m) {
  cat("max_abs_error_to_reference=", max(abs(scores - cpp_scores)), "\n", sep = "")
  stopifnot(max(abs(scores - cpp_scores)) <= 1e-12)
}
