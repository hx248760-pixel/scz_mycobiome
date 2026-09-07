library(readxl)
library(dplyr)

code_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = FALSE))
source(file.path(code_dir, "00_config.R"), chdir = TRUE)

clean_names <- function(x) {
  x <- gsub("[\r\n]+", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

sample_map <- read_tsv(PATHS$raw_sample_map)
sample_map$Sample <- normalize_sample_id(sample_map$Sample)
sample_map <- normalize_group_labels(sample_map)
sample_map$Subject_ID <- sub("-t2$", "", sample_map$Sample)
sample_map$Timepoint <- ifelse(sample_map$Group == "follow_up", "follow_up", "baseline")
sample_map <- sample_map[, c("Sample", "Group", "Subject_ID", "Timepoint")]

write_tsv(sample_map, PATHS$metadata_all)
write_tsv(sample_map %>% filter(Group %in% c("HC", "SCZ")), PATHS$metadata_baseline)

sample_info <- readxl::read_excel(PATHS$raw_sample_info, sheet = 1, skip = 1)
colnames(sample_info) <- clean_names(colnames(sample_info))
sample_info$Sample <- normalize_sample_id(sample_info$Sample)

clinical <- sample_info %>%
  mutate(
    Subject_ID = Sample,
    Timepoint = "baseline"
  ) %>%
  left_join(sample_map[, c("Sample", "Group")], by = "Sample", suffix = c("", ".map"))

if ("Group.map" %in% colnames(clinical)) {
  clinical$Group <- ifelse(is.na(clinical$Group), as.character(clinical$Group.map), as.character(clinical$Group))
  clinical$Group.map <- NULL
}

write_tsv(as.data.frame(clinical), PATHS$clinical_metadata)

fungi_species <- read_feature_table(PATHS$raw_fungi_species)
fungi_genus <- read_feature_table(PATHS$raw_fungi_genus)
bacteria_species <- read_feature_table(PATHS$raw_bacteria_species)

write_tsv(fungi_species, PATHS$fungi_species, row_names = TRUE, col.names = NA)
write_tsv(fungi_genus, PATHS$fungi_genus, row_names = TRUE, col.names = NA)
write_tsv(bacteria_species, PATHS$bacteria_species, row_names = TRUE, col.names = NA)

fungi_rf <- feature_to_rf_matrix(fungi_species, sample_map %>% filter(Group %in% c("HC", "SCZ")))
bacteria_rf <- feature_to_rf_matrix(bacteria_species, sample_map %>% filter(Group %in% c("HC", "SCZ")))
merged_rf <- dplyr::left_join(fungi_rf, bacteria_rf, by = "Sample", suffix = c("_fungi", "_bacteria"))

write_tsv(fungi_rf, PATHS$rf_fungi)
write_tsv(bacteria_rf, PATHS$rf_bacteria)
write_tsv(merged_rf, PATHS$rf_merged)

message("Prepared standardized downstream input files from the supplied raw metadata and abundance tables.")

