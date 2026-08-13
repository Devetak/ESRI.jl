library(GLcascade)
library(fastcascade)
library(Matrix)
library(parallel)

main <- function() {
  options(digits = 17)
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1"
  )

  matrix_path <- Sys.getenv("ESRI_IHS_MATRIX")
  if (!file.exists(matrix_path)) {
    stop("Set ESRI_IHS_MATRIX to the labeled EssMatIHS.csv file")
  }
  classification_name <- Sys.getenv("ESRI_CLASSIFICATION", "ihs")
  if (!classification_name %in% c("ihs", "legacy", "linear")) {
    stop("ESRI_CLASSIFICATION must be ihs, legacy, or linear")
  }
  n <- as.integer(Sys.getenv("ESRI_BENCHMARK_FIRMS", "10000"))
  workers <- as.integer(Sys.getenv("ESRI_BENCHMARK_WORKERS", "8"))
  convergence_tol <- as.numeric(Sys.getenv("ESRI_CONVERGENCE_TOL", "1e-2"))
  if (!is.finite(convergence_tol) || convergence_tol <= 0) {
    stop("ESRI_CONVERGENCE_TOL must be positive")
  }
  chunk_size <- as.integer(Sys.getenv("ESRI_CPP_CHUNK_SIZE", "256"))
  output_root <- Sys.getenv(
    "ESRI_OUTPUT_DIR",
    file.path("results", "full_esri_matrix_comparison")
  )
  network_path <- Sys.getenv("ESRI_NETWORK_FILE")
  if (!file.exists(network_path)) {
    stop("Set ESRI_NETWORK_FILE to the shared power-law edge CSV")
  }
  case_dir <- file.path(output_root, paste0(classification_name, "_", n))
  dir.create(case_dir, recursive = TRUE, showWarnings = FALSE)

  d <- read.csv(
    matrix_path,
    row.names = 1,
    check.names = FALSE,
    colClasses = "character"
  )
  ihs <- apply(d, 2L, as.integer)
  rownames(ihs) <- rownames(d)
  colnames(ihs) <- colnames(d)
  stopifnot(
    nrow(ihs) == ncol(ihs),
    identical(rownames(ihs), colnames(ihs)),
    all(ihs %in% 0:2)
  )

  m <- nrow(ihs)
  if (n < m) stop("ESRI_BENCHMARK_FIRMS must be at least ", m)
  if (classification_name == "ihs") {
    E <- ihs
  } else if (classification_name == "legacy") {
    essential_codes <- ifelse(seq_len(m) %% 2L == 0L, 2L, 1L)
    E <- matrix(rep(essential_codes, m), nrow = m, ncol = m)
  } else {
    E <- matrix(1L, nrow = m, ncol = m)
  }
  rownames(E) <- rownames(ihs)
  colnames(E) <- colnames(ihs)

  unique_E <- unique(E, MARGIN = 2L)
  expected_profiles <- if (classification_name == "ihs") 56L else 1L
  stopifnot(ncol(unique_E) == expected_profiles)
  profile <- vapply(seq_len(m), function(j) {
    which(vapply(
      seq_len(ncol(unique_E)),
      function(k) identical(E[, j], unique_E[, k]),
      logical(1)
    ))[1]
  }, integer(1))

  firms_per_industry <- n %/% m
  extra_firms <- n %% m
  p_cons <- c(rep(seq_len(m), each = firms_per_industry), seq_len(extra_firms))
  p <- rownames(E)[p_cons]
  network_values <- scan(network_path, what = double(), sep = ",", quiet = TRUE)
  stopifnot(length(network_values) == 24L * n)
  network_data <- matrix(network_values, ncol = 3L, byrow = TRUE)
  supplier <- as.integer(network_data[, 1L])
  customer <- as.integer(network_data[, 2L])
  weights <- network_data[, 3L]
  stopifnot(
    length(supplier) == 8L * n,
    all(supplier >= 1L & supplier <= n),
    all(customer >= 1L & customer <= n),
    all(supplier != customer),
    all(is.finite(weights) & weights > 0)
  )
  W <- sparseMatrix(
    i = supplier,
    j = customer,
    x = weights,
    dims = c(n, n)
  )
  stopifnot(length(W@x) == length(weights))
  rm(network_values, network_data, supplier, customer, weights)
  gc()

  psup <- sparseMatrix(i = seq_len(n), j = p_cons, x = 1, dims = c(n, m))
  in_str <- colSums(W)
  out_str <- rowSums(W)
  Lambda_d <- Matrix::drop0(Matrix(
    GLcascade::Lambda_d_calc_r(
      W,
      p,
      p_cons,
      E,
      crossprod(W, psup),
      in_str,
      n
    ),
    sparse = TRUE
  ))
  Lambda_u <- Matrix(W / ifelse(out_str == 0, 1, out_str), sparse = TRUE)
  zero_sector <- sparseMatrix(
    i = integer(),
    j = integer(),
    x = numeric(),
    dims = c(n, 1L)
  )

  run_native <- function(ids) {
    psi <- Matrix::sparseMatrix(
      i = ids,
      j = seq_along(ids),
      x = 1,
      dims = c(n, length(ids))
    )
    as.numeric(fastcascade::GL_cascade_dynamics_cpp(
      Lambda_d,
      Lambda_u,
      psi,
      psi,
      n,
      m,
      TRUE,
      psup,
      cbind(out_str),
      convergence_tol,
      FALSE,
      as.integer(p_cons - 1L),
      out_str,
      psup,
      FALSE,
      FALSE,
      FALSE,
      1,
      zero_sector,
      out_str,
      out_str,
      as.integer(seq_len(m)),
      as.integer(profile[p_cons] - 1L),
      unique_E
    )[, 1L])
  }
  run_batch <- function(task) {
    list(batch = task$batch, ids = task$ids, scores = run_native(task$ids))
  }

  starts <- seq.int(1L, n, by = chunk_size)
  tasks <- Map(
    function(batch, start) {
      list(batch = batch, ids = seq.int(start, min(start + chunk_size - 1L, n)))
    },
    seq_along(starts),
    starts
  )
  waves <- split(tasks, ceiling(seq_along(tasks) / (2L * workers)))
  config <- list(
    format_version = 1L,
    classification = classification_name,
    firms = n,
    workers = workers,
    chunk_size = chunk_size,
    convergence_tol = convergence_tol,
    network_md5 = unname(tools::md5sum(network_path))
  )
  config_path <- file.path(case_dir, "cpp_config.rds")
  if (file.exists(config_path)) {
    if (!identical(readRDS(config_path), config)) {
      stop("Existing C++ checkpoints use a different configuration")
    }
  } else {
    saveRDS(config, config_path)
  }

  cl <- parallel::makePSOCKcluster(workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, {
    library(fastcascade)
    library(Matrix)
    NULL
  })
  parallel::clusterExport(
    cl,
    c(
      "Lambda_d", "Lambda_u", "n", "m", "psup", "out_str", "p_cons",
      "profile", "unique_E", "zero_sector", "run_native", "run_batch"
    ),
    envir = environment()
  )
  invisible(parallel::clusterCall(cl, function() run_native(1L)))

  for (wave_index in seq_along(waves)) {
    checkpoint <- file.path(case_dir, sprintf("cpp_wave_%05d.rds", wave_index))
    if (file.exists(checkpoint)) next
    started <- proc.time()[["elapsed"]]
    parts <- parallel::parLapplyLB(cl, waves[[wave_index]], run_batch)
    elapsed <- unname(proc.time()[["elapsed"]] - started)
    tmp <- paste0(checkpoint, ".", Sys.getpid(), ".tmp")
    saveRDS(list(elapsed = elapsed, parts = parts), tmp)
    if (!file.rename(tmp, checkpoint)) stop("Could not save ", checkpoint)
    cat(
      "classification=", classification_name,
      " firms=", n,
      " wave=", wave_index, "/", length(waves),
      " elapsed_s=", elapsed, "\n",
      sep = ""
    )
    flush.console()
  }

  completed <- lapply(
    seq_along(waves),
    function(i) readRDS(file.path(case_dir, sprintf("cpp_wave_%05d.rds", i)))
  )
  results <- unlist(lapply(completed, `[[`, "parts"), recursive = FALSE)
  results <- results[order(vapply(results, `[[`, integer(1), "batch"))]
  ids <- unlist(lapply(results, `[[`, "ids"), use.names = FALSE)
  scores <- unlist(lapply(results, `[[`, "scores"), use.names = FALSE)
  stopifnot(identical(ids, seq_len(n)), length(scores) == n, all(is.finite(scores)))
  cpp_total_s <- sum(vapply(completed, `[[`, numeric(1), "elapsed"))

  write.table(
    cbind(ids, scores),
    file = file.path(case_dir, "cpp_scores.csv"),
    sep = ",",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  writeLines(
    c(
      paste0("classification=", classification_name),
      paste0("firms=", n),
      paste0("workers=", workers),
      paste0("cpp_total_s=", format(cpp_total_s, digits = 17)),
      paste0("cpp_total_esri=", format(sum(scores), digits = 17))
    ),
    file.path(case_dir, "cpp_summary.txt")
  )
  cat(
    "classification=", classification_name,
    " firms=", n,
    " workers=", workers,
    " cpp_total_s=", cpp_total_s,
    " total_esri=", sum(scores), "\n",
    sep = ""
  )
}

main()
