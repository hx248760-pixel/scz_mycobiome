library(dplyr)
library(ggplot2)
library(vegan)

code_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = FALSE))
source(file.path(code_dir, "00_config.R"), chdir = TRUE)
ensure_prepared_inputs()

run_pcoa <- function(feature_table, sample_map, output_prefix) {
  aligned <- align_feature_table(feature_table, sample_map)
  x <- t(aligned$feature_table)
  dist_mat <- vegan::vegdist(x, method = "bray")
  pcoa <- cmdscale(dist_mat, k = 2, eig = TRUE)
  scores <- data.frame(Sample = rownames(pcoa$points), PCoA1 = pcoa$points[, 1], PCoA2 = pcoa$points[, 2])
  scores <- left_join(scores, aligned$sample_map, by = "Sample")
  permanova <- vegan::adonis2(dist_mat ~ Group, data = aligned$sample_map, permutations = 999)

  write_tsv(scores, file.path(RESULTS_DIR, "tables", paste0(output_prefix, "_pcoa_scores.tsv")))
  write_tsv(as.data.frame(permanova), file.path(RESULTS_DIR, "tables", paste0(output_prefix, "_permanova.tsv")), row_names = TRUE)

  p <- ggplot(scores, aes(PCoA1, PCoA2, color = Group)) +
    geom_point(size = 2, alpha = 0.85) +
    stat_ellipse(level = 0.68, linewidth = 0.5, show.legend = FALSE) +
    scale_color_manual(values = GROUP_COLORS, drop = FALSE) +
    theme_classic(base_size = 12) +
    labs(x = "PCoA1", y = "PCoA2")
  save_plot(p, file.path(RESULTS_DIR, "figures", paste0(output_prefix, "_pcoa.pdf")), width = 5, height = 4.5)
}

fungi <- read_feature_table(PATHS$fungi_species)
run_pcoa(fungi, normalize_group_labels(read_tsv(PATHS$metadata_baseline)), "figure1")
run_pcoa(fungi, normalize_group_labels(read_tsv(PATHS$metadata_all)), "figure4")

