# Compatibility wrapper for the original all-sample PCoA script.
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = FALSE)), "PCoA.R"), chdir = TRUE)

