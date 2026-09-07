library(dplyr)

code_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = FALSE))
source(file.path(code_dir, "00_config.R"), chdir = TRUE)
ensure_prepared_inputs()

if (!requireNamespace("Maaslin2", quietly = TRUE)) {
  stop("MaAsLin2 is not installed. Install MaAsLin2 or run the other downstream scripts without this optional sensitivity analysis.", call. = FALSE)
}

run_maaslin <- function(feature_table, sample_map, output_dir) {
  aligned <- align_feature_table(feature_table, sample_map)
  features <- as.data.frame(t(aligned$feature_table), check.names = FALSE)
  metadata <- data.frame(Group = aligned$sample_map$Group, row.names = aligned$sample_map$Sample)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  Maaslin2::Maaslin2(
    input_data = features,
    input_metadata = metadata,
    output = output_dir,
    fixed_effects = c("Group"),
    normalization = "NONE",
    transform = "LOG",
    analysis_method = "LM",
    correction = "BH",
    standardize = FALSE,
    plot_heatmap = FALSE,
    plot_scatter = FALSE
  )
}

sample_map <- normalize_group_labels(read_tsv(PATHS$metadata_baseline))
fungi <- read_feature_table(PATHS$fungi_species)
run_maaslin(fungi, sample_map, file.path(RESULTS_DIR, "maaslin2", "SCZ_vs_HC_species"))

all_meta <- normalize_group_labels(read_tsv(PATHS$metadata_all))
run_maaslin(fungi, all_meta, file.path(RESULTS_DIR, "maaslin2", "followup_species"))

