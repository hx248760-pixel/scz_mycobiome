# Compatibility wrapper for the original alpha-diversity helper.
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = FALSE)), "alpha_diversity.R"), chdir = TRUE)

