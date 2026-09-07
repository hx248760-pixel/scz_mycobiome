library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(ggplot2)

code_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = FALSE))
source(file.path(code_dir, "00_config.R"), chdir = TRUE)
ensure_prepared_inputs()
if (!file.exists(PATHS$marker_species) || !file.exists(PATHS$candidate_bacteria)) {
  source(file.path(CODE_DIR, "wilcoxon_daa.R"), chdir = TRUE)
}

read_candidate_names <- function(path) {
  if (!file.exists(path)) return(character(0))
  tab <- read_tsv(path)
  if ("name" %in% colnames(tab)) tab$name else tab[[1]]
}

matrix_by_features <- function(abundance, features, sample_map) {
  if (length(features) > 0) abundance <- abundance[intersect(features, rownames(abundance)), , drop = FALSE]
  aligned <- align_feature_table(abundance, sample_map)
  out <- as.data.frame(t(aligned$feature_table), check.names = FALSE)
  out$Sample <- rownames(out)
  out
}

pairwise_spearman <- function(x, y, source_label) {
  out <- list()
  n <- 1
  for (a in colnames(x)) {
    for (b in colnames(y)) {
      ok <- complete.cases(x[[a]], y[[b]])
      if (sum(ok) >= 5) {
        test <- suppressWarnings(cor.test(x[[a]][ok], y[[b]][ok], method = "spearman", exact = FALSE))
        out[[n]] <- data.frame(A = a, B = b, corr = unname(test$estimate), pvalue = test$p.value,
                               source = source_label, stringsAsFactors = FALSE)
        n <- n + 1
      }
    }
  }
  bind_rows(out)
}

feature_filter <- function(tab, min_frequency = 0.10, min_mean = 1e-5) {
  keep <- colMeans(tab > 0, na.rm = TRUE) > min_frequency & colMeans(tab, na.rm = TRUE) > min_mean
  tab[, keep, drop = FALSE]
}

build_group_network <- function(group_name, fungi_tab, bacteria_tab, clinical_tab) {
  samples <- Reduce(intersect, list(fungi_tab$Sample, bacteria_tab$Sample, clinical_tab$Sample))
  fungi <- fungi_tab %>% filter(Sample %in% samples) %>% arrange(match(Sample, samples)) %>% select(-Sample)
  bacteria <- bacteria_tab %>% filter(Sample %in% samples) %>% arrange(match(Sample, samples)) %>% select(-Sample)
  clinical <- clinical_tab %>% filter(Sample %in% samples) %>% arrange(match(Sample, samples)) %>% select(-Sample, -Group, -Subject_ID, -Timepoint, any_of("Diagnosis"))
  clinical <- clinical[, vapply(clinical, is.numeric, logical(1)), drop = FALSE]
  fungi <- feature_filter(fungi)
  bacteria <- feature_filter(bacteria)

  edges_all <- bind_rows(
    pairwise_spearman(bacteria, fungi, "bacteria_fungi"),
    pairwise_spearman(clinical, fungi, "clinical_fungi"),
    pairwise_spearman(bacteria, clinical, "bacteria_clinical")
  ) %>%
    mutate(qvalue = p.adjust(pvalue, method = "BH"), Group = group_name)

  edges_filtered <- edges_all %>% filter(abs(corr) >= 0.3, qvalue < 0.05)
  graph <- graph_from_data_frame(edges_filtered[, c("A", "B")], directed = FALSE)
  nodes <- data.frame(name = V(graph)$name, degree = degree(graph), Group = group_name) %>%
    mutate(type = case_when(
      name %in% colnames(fungi) ~ "fungi",
      name %in% colnames(bacteria) ~ "bacteria",
      TRUE ~ "clinical"
    ))
  list(all_edges = edges_all, filtered_edges = edges_filtered, nodes = nodes)
}

sample_map <- normalize_group_labels(read_tsv(PATHS$metadata_baseline))
clinical <- read_tsv(PATHS$clinical_metadata) %>% mutate(Sample = normalize_sample_id(Sample))
clinical <- inner_join(sample_map, clinical, by = c("Sample", "Group"), suffix = c("", ".clinical"))

fungi <- read_feature_table(PATHS$fungi_species)
bacteria <- read_feature_table(PATHS$bacteria_species)
fungi_tab <- matrix_by_features(fungi, read_candidate_names(PATHS$marker_species), sample_map)
bacteria_tab <- matrix_by_features(bacteria, read_candidate_names(PATHS$candidate_bacteria), sample_map)

all_edges <- list()
filtered_edges <- list()
node_stats <- list()

for (grp in c("HC", "SCZ")) {
  smp <- sample_map$Sample[sample_map$Group == grp]
  res <- build_group_network(
    grp,
    fungi_tab %>% filter(Sample %in% smp),
    bacteria_tab %>% filter(Sample %in% smp),
    clinical %>% filter(Sample %in% smp)
  )
  all_edges[[grp]] <- res$all_edges
  filtered_edges[[grp]] <- res$filtered_edges
  node_stats[[grp]] <- res$nodes

  if (nrow(res$filtered_edges) > 0) {
    graph <- graph_from_data_frame(res$filtered_edges[, c("A", "B", "corr", "source")],
                                   directed = FALSE, vertices = res$nodes)
    p <- ggraph(graph, layout = "fr") +
      geom_edge_link(aes(color = corr > 0, width = abs(corr)), alpha = 0.45) +
      geom_node_point(aes(fill = type), shape = 21, size = 2.5) +
      theme_void() +
      ggtitle(paste0(grp, " exploratory correlation network"))
    save_plot(p, file.path(RESULTS_DIR, "figures", paste0("figure2_", grp, "_network.pdf")), width = 7, height = 6)
  }
}

write_tsv(bind_rows(all_edges), file.path(RESULTS_DIR, "tables", "network_all_tested_pairs.tsv"))
write_tsv(bind_rows(filtered_edges), file.path(RESULTS_DIR, "tables", "network_filtered_edges.tsv"))
write_tsv(bind_rows(node_stats), file.path(RESULTS_DIR, "tables", "network_node_degrees.tsv"))

