library(dplyr)
library(tidyr)
library(ggplot2)

code_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = FALSE))
source(file.path(code_dir, "00_config.R"), chdir = TRUE)
ensure_prepared_inputs()
if (!file.exists(PATHS$marker_species)) {
  source(file.path(CODE_DIR, "wilcoxon_daa.R"), chdir = TRUE)
}

fungi <- read_feature_table(PATHS$fungi_species)
marker <- read_tsv(PATHS$marker_species)
sample_map <- normalize_group_labels(read_tsv(PATHS$metadata_all))
candidate_taxa <- intersect(marker$name, rownames(fungi))
if (length(candidate_taxa) == 0) stop("No candidate fungal taxa are available for plotting.", call. = FALSE)

aligned <- align_feature_table(fungi[candidate_taxa, , drop = FALSE], sample_map)
plot_df <- as.data.frame(t(aligned$feature_table), check.names = FALSE)
plot_df$Sample <- rownames(plot_df)
plot_df <- plot_df %>%
  tidyr::pivot_longer(-Sample, names_to = "Taxon", values_to = "Abundance") %>%
  left_join(aligned$sample_map, by = "Sample")

p <- ggplot(plot_df, aes(Group, Abundance + 1e-9, fill = Group)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  facet_wrap(~ Taxon, scales = "free_y") +
  scale_y_log10() +
  scale_fill_manual(values = GROUP_COLORS, drop = FALSE) +
  theme_classic(base_size = 10) +
  labs(x = NULL, y = "Relative abundance + 1e-9")
save_plot(p, file.path(RESULTS_DIR, "figures", "figure4E_marker_boxplots.pdf"), width = 10, height = 7)

paired_tests <- plot_df %>%
  filter(Group %in% c("SCZ", "follow_up")) %>%
  mutate(pair_id = sub("-t2$", "", Sample)) %>%
  group_by(Taxon, pair_id) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  group_by(Taxon) %>%
  summarise(p.value = suppressWarnings(wilcox.test(Abundance ~ Group, paired = TRUE, exact = FALSE)$p.value),
            .groups = "drop") %>%
  mutate(q.value = p.adjust(p.value, method = "BH"))
write_tsv(paired_tests, file.path(RESULTS_DIR, "tables", "followup_candidate_taxa_paired_tests.tsv"))

