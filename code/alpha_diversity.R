library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)

code_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = FALSE))
source(file.path(code_dir, "00_config.R"), chdir = TRUE)
ensure_prepared_inputs()

calc_alpha <- function(feature_table, sample_map) {
  aligned <- align_feature_table(feature_table, sample_map)
  x <- t(aligned$feature_table)
  data.frame(
    Sample = rownames(x),
    Shannon = vegan::diversity(x, index = "shannon"),
    Simpson = vegan::diversity(x, index = "simpson"),
    Observed_species = rowSums(x > 0),
    stringsAsFactors = FALSE
  ) %>%
    left_join(aligned$sample_map, by = "Sample")
}

plot_alpha <- function(alpha_df, metric, output, groups = NULL) {
  df <- alpha_df
  if (!is.null(groups)) df <- df %>% filter(Group %in% groups)
  df$Group <- droplevels(df$Group)
  p <- ggplot(df, aes(x = Group, y = .data[[metric]], fill = Group)) +
    geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 1.2, alpha = 0.65) +
    scale_fill_manual(values = GROUP_COLORS, drop = FALSE) +
    theme_classic(base_size = 12) +
    labs(x = NULL, y = metric)
  save_plot(p, output, width = 4.2, height = 4)
}

fungi <- read_feature_table(PATHS$fungi_species)
baseline <- normalize_group_labels(read_tsv(PATHS$metadata_baseline))
all_meta <- normalize_group_labels(read_tsv(PATHS$metadata_all))

alpha_baseline <- calc_alpha(fungi, baseline)
alpha_all <- calc_alpha(fungi, all_meta)

write_tsv(alpha_baseline, file.path(RESULTS_DIR, "tables", "alpha_diversity_baseline.tsv"))
write_tsv(alpha_all, file.path(RESULTS_DIR, "tables", "alpha_diversity_all.tsv"))

plot_alpha(alpha_baseline, "Shannon", file.path(RESULTS_DIR, "figures", "figure1B_shannon.pdf"), c("HC", "SCZ"))
plot_alpha(alpha_baseline, "Simpson", file.path(RESULTS_DIR, "figures", "figure1C_simpson.pdf"), c("HC", "SCZ"))
plot_alpha(alpha_baseline, "Observed_species", file.path(RESULTS_DIR, "figures", "figure1D_observed_species.pdf"), c("HC", "SCZ"))
plot_alpha(alpha_all, "Shannon", file.path(RESULTS_DIR, "figures", "figure4A_shannon.pdf"))
plot_alpha(alpha_all, "Simpson", file.path(RESULTS_DIR, "figures", "figure4B_simpson.pdf"))
plot_alpha(alpha_all, "Observed_species", file.path(RESULTS_DIR, "figures", "figure4C_observed_species.pdf"))

