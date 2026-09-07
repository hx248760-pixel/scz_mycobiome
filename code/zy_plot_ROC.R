# Compatibility wrapper for the original ROC plotting script.
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = FALSE)), "plot_randomforest_results.R"), chdir = TRUE)

