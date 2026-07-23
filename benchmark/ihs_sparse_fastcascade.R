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
n <- 2L * m
p <- rep(rownames(E), each = 2L)
p_cons <- match(p, rownames(E))
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
psi <- sparseMatrix(i = seq_len(16L), j = seq_len(16L), x = 1, dims = c(n, 16L))
zero_sector <- sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = c(n, 1L))

run_cpp <- function() fastcascade::GL_cascade_dynamics_cpp(
  Lambda_d, Lambda_u, psi, psi, n, m, TRUE, psup, cbind(out_str), 1e-2, FALSE,
  as.integer(p_cons - 1L), out_str, psup, FALSE, FALSE, FALSE, 1,
  zero_sector, out_str, out_str, as.integer(seq_len(m)),
  as.integer(profile[p_cons] - 1L), unique_E
)

cpp_scores <- c(
  0.87837245554794330, 0.87825862250992404, 0.87664923655232285, 0.87566688889824862,
  0.87524585194239457, 0.87433651919608424, 0.86060405584903987, 0.84701627242865141,
  0.85599031156433647, 0.85741739269669437, 0.85466742541113849, 0.85393912577670528,
  0.84935069451138689, 0.77626516574715354, 0.84447307983109421, 0.84456652158172307
)

invisible(run_cpp())
scores <- as.numeric(run_cpp()[, 1L])
samples <- vapply(
  seq_len(as.integer(Sys.getenv("ESRI_BENCHMARK_SAMPLES", "5"))),
  function(i) system.time(invisible(run_cpp()))[["elapsed"]],
  numeric(1)
)
cat("firms=", n, " industries=", m, " edges=", length(W@x), " closures=16\n", sep = "")
cat("downstream_nnz=", length(Lambda_d@x), "\n", sep = "")
cat("cpp_prepared_median_s=", median(samples), "\n", sep = "")
cat("max_abs_error_to_reference=", max(abs(scores - cpp_scores)), "\n", sep = "")
stopifnot(max(abs(scores - cpp_scores)) <= 1e-12)
