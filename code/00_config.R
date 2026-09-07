# Shared configuration for the SCZ gut mycobiome downstream analysis workflow.

get_code_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)))
  }
  if (file.exists("00_config.R")) {
    return(normalizePath(getwd(), mustWork = FALSE))
  }
  if (file.exists(file.path("code", "00_config.R"))) {
    return(normalizePath(file.path(getwd(), "code"), mustWork = FALSE))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

CODE_DIR <- get_code_dir()
REPO_DIR <- normalizePath(file.path(CODE_DIR, ".."), mustWork = FALSE)
DATA_DIR <- file.path(REPO_DIR, "data")
RAW_DIR <- file.path(DATA_DIR, "metadata")
PROCESSED_DIR <- file.path(DATA_DIR, "processed")
RF_DATA_DIR <- file.path(DATA_DIR, "random_forest")
RESULTS_DIR <- file.path(REPO_DIR, "results")

# Only these five files must be supplied by the user/repository.
RAW_INPUTS <- list(
  raw_sample_map = file.path(RAW_DIR, "sample_map"),
  raw_sample_info = file.path(RAW_DIR, "sample_info.xlsx"),
  raw_fungi_species = file.path(RAW_DIR, "fungi.profile"),
  raw_fungi_genus = file.path(RAW_DIR, "fungi.profile.genus"),
  raw_bacteria_species = file.path(RAW_DIR, "bacteria_species.profile")
)

# The files below are generated from RAW_INPUTS by code/01_prepare_inputs.R
# or by downstream analysis scripts. They do not need to be supplied manually.
GENERATED_PATHS <- list(
  metadata_baseline = file.path(DATA_DIR, "metadata", "sample_map_baseline.tsv"),
  metadata_all = file.path(DATA_DIR, "metadata", "sample_map_all.tsv"),
  clinical_metadata = file.path(DATA_DIR, "metadata", "clinical_metadata.tsv"),
  fungi_species = file.path(PROCESSED_DIR, "fungi.profile.species.tsv"),
  fungi_genus = file.path(PROCESSED_DIR, "fungi.profile.genus.tsv"),
  bacteria_species = file.path(PROCESSED_DIR, "bacteria.profile.species.tsv"),
  bacteria_genus = file.path(PROCESSED_DIR, "bacteria.profile.genus.tsv"),
  marker_species = file.path(PROCESSED_DIR, "marker_sp.tsv"),
  candidate_bacteria = file.path(PROCESSED_DIR, "candidate_bacteria.tsv"),
  rf_fungi = file.path(RF_DATA_DIR, "fungi_features.tsv"),
  rf_bacteria = file.path(RF_DATA_DIR, "bacteria_features.tsv"),
  rf_merged = file.path(RF_DATA_DIR, "merged_features.tsv")
)

PATHS <- c(RAW_INPUTS, GENERATED_PATHS)

GROUP_COLORS <- c("SCZ" = "#ea7827", "HC" = "#3288bd", "follow_up" = "#7f7f7f")

dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RF_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RESULTS_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RESULTS_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)

normalize_sample_id <- function(x) {
  x <- as.character(x)
  x <- gsub("\\.t2$", "-t2", x)
  x
}

normalize_group_labels <- function(sample_map, group_col = "Group") {
  sample_map[[group_col]] <- as.character(sample_map[[group_col]])
  sample_map[[group_col]][sample_map[[group_col]] == "flow_up"] <- "follow_up"
  sample_map[[group_col]] <- factor(sample_map[[group_col]], levels = c("HC", "SCZ", "follow_up"))
  sample_map
}

read_tsv <- function(path, row_names = NULL, check_exists = TRUE, ...) {
  if (check_exists && !file.exists(path)) {
    stop("Required file does not exist: ", path,
         "\nIf this is a generated file, run: Rscript code/01_prepare_inputs.R",
         call. = FALSE)
  }
  read.table(path, sep = "\t", header = TRUE, row.names = row_names,
             check.names = FALSE, comment.char = "", quote = "", ...)
}

write_tsv <- function(x, path, row_names = FALSE, ...) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = row_names, ...)
}

read_feature_table <- function(path) {
  tab <- read_tsv(path, row_names = 1)
  colnames(tab) <- normalize_sample_id(colnames(tab))
  tab[] <- lapply(tab, as.numeric)
  tab <- as.data.frame(tab, check.names = FALSE)
  tab[rowSums(tab, na.rm = TRUE) > 0, , drop = FALSE]
}

align_feature_table <- function(feature_table, sample_map, sample_col = "Sample") {
  sample_map[[sample_col]] <- normalize_sample_id(sample_map[[sample_col]])
  colnames(feature_table) <- normalize_sample_id(colnames(feature_table))
  samples <- intersect(colnames(feature_table), sample_map[[sample_col]])
  if (length(samples) == 0) {
    stop("No overlapping sample IDs between feature table and metadata.", call. = FALSE)
  }
  sample_map <- sample_map[match(samples, sample_map[[sample_col]]), , drop = FALSE]
  feature_table <- feature_table[, samples, drop = FALSE]
  feature_table <- feature_table[rowSums(feature_table, na.rm = TRUE) > 0, , drop = FALSE]
  list(feature_table = feature_table, sample_map = sample_map)
}

feature_to_rf_matrix <- function(feature_table, sample_map) {
  aligned <- align_feature_table(feature_table, sample_map)
  x <- as.data.frame(t(aligned$feature_table), check.names = FALSE)
  x <- cbind(Sample = rownames(x), x)
  rownames(x) <- NULL
  x
}

save_plot <- function(plot, path, width = 6, height = 5) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot = plot, width = width, height = height)
}

ensure_prepared_inputs <- function() {
  generated_files <- c(PATHS$metadata_baseline, PATHS$metadata_all, PATHS$clinical_metadata,
                       PATHS$fungi_species, PATHS$fungi_genus, PATHS$bacteria_species,
                       PATHS$rf_fungi, PATHS$rf_bacteria, PATHS$rf_merged)
  if (!all(file.exists(generated_files))) {
    source(file.path(CODE_DIR, "01_prepare_inputs.R"), chdir = TRUE)
  }
}
