args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: Rscript benchmark/esri_tutorial_prepare.R <tutorial_root> <output_dir>")
}

tutorial_root <- normalizePath(args[1], mustWork = TRUE)
output_dir <- args[2]
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages(library(GLcascade))
suppressPackageStartupMessages(library(Matrix))

data_path <- file.path(tutorial_root, "data", "2020.SML.csv")
icio <- read.csv(data_path, check.names = FALSE, stringsAsFactors = FALSE)
rownames(icio) <- icio[[1]]
icio <- icio[, -1, drop = FALSE]

row_stop <- grep("^ROW_T$", rownames(icio))[1]
col_stop <- grep("^ROW_T$", colnames(icio))[1]
if (is.na(row_stop) || is.na(col_stop) || row_stop != col_stop) {
  stop("failed to align ROW_T boundary in 2020.SML.csv")
}

W_icio <- as.matrix(icio[seq_len(row_stop), seq_len(col_stop), drop = FALSE])
W_icio[W_icio <= 1] <- 0
W_icio <- Matrix::Matrix(W_icio, sparse = TRUE)

country_sector_pair <- rownames(W_icio)
sector_codes <- substr(country_sector_pair, 5, nchar(country_sector_pair))
sector_prefix <- substr(sector_codes, 1, 3)

nace_lookup <- as.character(GLcascade::nace_conv_mat[, "nace1_2"])
nace4_num <- as.character(GLcascade::nace_conv_mat[, "nace4_num"])
unique_prefix <- sort(unique(sector_prefix))
prefix_to_nace4 <- sapply(unique_prefix, function(prefix) {
  matches <- grep(prefix, nace_lookup, fixed = TRUE)
  if (length(matches) == 0) {
    stop(sprintf("failed to map ICIO sector prefix %s into nace_conv_mat", prefix))
  }
  min(matches)
})
names(prefix_to_nace4) <- unique_prefix

sector_nace2 <- substr(nace4_num[prefix_to_nace4[sector_prefix]], 1, 2)
sector_levels <- sort(unique(sector_nace2))
industry_ids <- match(sector_nace2, sector_levels)

row_weights <- Matrix::rowSums(W_icio)
col_weights <- Matrix::colSums(W_icio)
edge_summary <- Matrix::summary(W_icio)

write.table(
  cbind(edge_summary$i, edge_summary$j, edge_summary$x),
  file.path(output_dir, "edges.tsv"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
write.table(industry_ids, file.path(output_dir, "industries.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(sector_levels, file.path(output_dir, "industry_labels.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(country_sector_pair, file.path(output_dir, "firm_labels.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(row_weights, file.path(output_dir, "row_weights.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(col_weights, file.path(output_dir, "col_weights.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

metadata <- data.frame(
  key = c("n", "nnz", "nindustries"),
  value = c(nrow(W_icio), length(edge_summary$x), length(sector_levels))
)
write.table(metadata, file.path(output_dir, "metadata.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
