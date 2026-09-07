library(dplyr)
library(ggplot2)

code_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = FALSE))
source(file.path(code_dir, "00_config.R"), chdir = TRUE)

auc_files <- list.files(file.path(RESULTS_DIR, "random_forest"), pattern = "\\.auc\\.tsv$", full.names = TRUE)
if (length(auc_files) == 0) stop("No random-forest AUC files found. Run code/randomforest.txt first.", call. = FALSE)
auc <- bind_rows(lapply(auc_files, read_tsv))
write_tsv(auc, file.path(RESULTS_DIR, "tables", "random_forest_auc_summary.tsv"))

p <- ggplot(auc, aes(x = model, y = auc, fill = model)) +
  geom_col(width = 0.65) +
  ylim(0, 1) +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Cross-validated AUC")
save_plot(p, file.path(RESULTS_DIR, "figures", "figure3_random_forest_auc.pdf"), width = 4.5, height = 4)

