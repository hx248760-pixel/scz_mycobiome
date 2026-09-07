library(dplyr)

code_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = FALSE))
source(file.path(code_dir, "00_config.R"), chdir = TRUE)
ensure_prepared_inputs()

run_wilcoxon_daa <- function(feature_table, sample_map, case = "SCZ", control = "HC", pseudocount = 1e-9) {
  sample_map <- normalize_group_labels(sample_map)
  aligned <- align_feature_table(feature_table, sample_map)
  dt <- aligned$feature_table
  meta <- aligned$sample_map
  case_samples <- meta$Sample[meta$Group == case]
  control_samples <- meta$Sample[meta$Group == control]

  res <- lapply(rownames(dt), function(feature) {
    x_case <- as.numeric(dt[feature, case_samples])
    x_control <- as.numeric(dt[feature, control_samples])
    avg_case <- mean(x_case, na.rm = TRUE)
    avg_control <- mean(x_control, na.rm = TRUE)
    test <- suppressWarnings(wilcox.test(x_case, x_control, exact = FALSE))
    fold_change <- ifelse(avg_case >= avg_control,
                          (avg_case + pseudocount) / (avg_control + pseudocount),
                          (avg_control + pseudocount) / (avg_case + pseudocount))
    data.frame(
      name = feature,
      Avg.HC = avg_control,
      Avg.SCZ = avg_case,
      log2FC = log2((avg_case + pseudocount) / (avg_control + pseudocount)),
      fold_change = fold_change,
      enriched = ifelse(avg_case >= avg_control, case, control),
      pvalue = test$p.value,
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows() %>%
    mutate(qvalue = p.adjust(pvalue, method = "BH")) %>%
    arrange(pvalue)
  res
}

sample_map <- normalize_group_labels(read_tsv(PATHS$metadata_baseline))

fungi_species_res <- run_wilcoxon_daa(read_feature_table(PATHS$fungi_species), sample_map)
fungi_species_candidates <- fungi_species_res %>% filter(pvalue < 0.05, fold_change > 1.2)
write_tsv(fungi_species_res, file.path(RESULTS_DIR, "tables", "wilcoxon_fungal_species_all.tsv"))
write_tsv(fungi_species_candidates, file.path(RESULTS_DIR, "tables", "wilcoxon_fungal_species_candidates.tsv"))
write_tsv(fungi_species_candidates["name"], PATHS$marker_species)

fungi_genus_res <- run_wilcoxon_daa(read_feature_table(PATHS$fungi_genus), sample_map)
write_tsv(fungi_genus_res, file.path(RESULTS_DIR, "tables", "wilcoxon_fungal_genus_all.tsv"))
write_tsv(fungi_genus_res %>% filter(pvalue < 0.05, fold_change > 1.2),
          file.path(RESULTS_DIR, "tables", "wilcoxon_fungal_genus_candidates.tsv"))

bacteria_species_res <- run_wilcoxon_daa(read_feature_table(PATHS$bacteria_species), sample_map)
bacteria_candidates <- bacteria_species_res %>% filter(pvalue < 0.05, fold_change > 1.2)
write_tsv(bacteria_species_res, file.path(RESULTS_DIR, "tables", "wilcoxon_bacterial_species_all.tsv"))
write_tsv(bacteria_candidates, file.path(RESULTS_DIR, "tables", "wilcoxon_bacterial_species_candidates.tsv"))
write_tsv(bacteria_candidates["name"], PATHS$candidate_bacteria)

