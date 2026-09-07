# Compatibility wrapper for random-forest AUC plotting.
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = FALSE)), "plot_randomforest_results.R"), chdir = TRUE)

