# =============================================================================
# microbiome_quick_compare.R
# Author: Catarina A Bettencourt
# Date: 02/27/2026
#
# Description:
#   Interactive pipeline for comparing gut microbiome composition between a
#   disease group and healthy controls using curatedMetagenomicData.
#   Supports relative abundance, pathway abundance, and pathway coverage.
#   Outputs plots (.png) and statistical results (.csv) to a timestamped
#   project folder for downstream reporting.
#
# Usage:
#   Run interactively in RStudio, or from the terminal:
#     Rscript microbiome_quick_compare.R
#
# Outputs (written to ./<PROJECT_ID>/):
#   data/     — cached phyloseq objects (.rds)
#   plots/    — diversity and top-feature plots (.png)
#   results/  — statistical tables (.csv)
#
# Notes:
#   Creating the phyloseq objects is extremely taxing on a laptop. Sample data
#   available on github. Formatting suggestions provided by Claude.
# =============================================================================


# -- 0. Dependencies ----------------------------------------------------------

packages <- c(
  "curatedMetagenomicData", "ggplot2", "tidyverse", "vegan",
  "phyloseq", "ggsignif", "ggpubr", "cowplot"
)

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

citation("curatedMetagenomicData")

# -- 1. User Inputs -----------------------------------------------------------

resume_choice <- menu(
  c("Start new analysis", "Resume existing project"),
  title = "New or existing project?"
)

if (resume_choice == 2) {
  FOLDER_PATH   <- readline(prompt = "Paste project folder name: ")
  config        <- readRDS(file.path(FOLDER_PATH, "data", "config.rds"))
  DIAGNOSIS         <- config$DIAGNOSIS
  SAMPLE_LOCATION   <- config$SAMPLE_LOCATION
  data_to_pull      <- config$data_to_pull
  acceptable_abx    <- config$acceptable_abx
  acceptable_gender <- config$acceptable_gender
  MIN_READS         <- config$MIN_READS
  MIN_BMI           <- config$MIN_BMI
  MAX_BMI           <- config$MAX_BMI
  meta_combined     <- readRDS(file.path(FOLDER_PATH, "data", "meta_combined.rds"))
  cat("Resumed project:", FOLDER_PATH, "\n")
  
} else {

  meta             <- sampleMetadata
  diag_options     <- unique(meta$disease)
  loc_options      <- unique(meta$body_site)
  
  date             <- readline(prompt = "Enter date (YYYYMMDD): ")
  DIAGNOSIS        <- diag_options[menu(diag_options, title = "Select Diagnosis:")]
  SAMPLE_LOCATION  <- loc_options[menu(loc_options,   title = "Select Sample Location:")]
  
  # Data type selection
  # NOTE: will add "gene_families" once server access is available
  data_options <- list(
    "relative_abundance",
    "pathway_abundance",
    "pathway_coverage",
    list("pathway_abundance", "pathway_coverage"),
    list("relative_abundance", "pathway_abundance", "pathway_coverage")
  )
  
  data_option_labels <- c(
    "Relative Abundance only",
    "Pathway Abundance only",
    "Pathway Coverage only",
    "Pathway Abundance + Pathway Coverage",
    "All: Relative Abundance + Pathway Abundance + Pathway Coverage"
  )
  
  data_choice  <- menu(data_option_labels, title = "Select Types of Data to Pull:")
  data_to_pull <- data_options[[data_choice]]
  if (is.character(data_to_pull)) data_to_pull <- list(data_to_pull)
  
  # Antibiotic filter
  abx_choice     <- menu(
    choices = c("Only Antibiotic Use", "Exclude Antibiotic Use", "Include All"),
    title   = "Antibiotic Filtering:"
  )
  acceptable_abx <- switch(abx_choice, c("yes"), c("no", NA), c("yes", "no", NA))
  abx_label      <- switch(abx_choice, "abx-only", "no-abx", "all-pts")
  
  # Gender filter
  gen_choice        <- menu(c("Male", "Female", "Include All"), title = "Gender Filter:")
  acceptable_gender <- switch(gen_choice, "male", "female", unique(meta$gender))
  gen_selected      <- switch(gen_choice, "sex-M", "sex-F", "sex-All")
  
  # Numeric filters
  input_reads   <- readline(prompt = "Min Read Length (press Enter for 0): ")
  MIN_READS     <- if (input_reads   == "") 0   else as.integer(input_reads)
  
  input_min_bmi <- readline(prompt = "Min BMI (press Enter for 0): ")
  MIN_BMI       <- if (input_min_bmi == "") 0   else as.integer(input_min_bmi)
  
  input_max_bmi <- readline(prompt = "Max BMI (press Enter for no max): ")
  MAX_BMI       <- if (input_max_bmi == "") Inf else as.integer(input_max_bmi)
  
  
  # -- 2. Cohort Filtering ------------------------------------------------------
  
  disease_data <- meta %>% filter(
    disease                 == DIAGNOSIS,
    body_site               == SAMPLE_LOCATION,
    antibiotics_current_use %in% acceptable_abx,
    gender                  %in% acceptable_gender,
    BMI                     >= MIN_BMI,
    BMI                     <= MAX_BMI,
    minimum_read_length     >= MIN_READS
  )
  
  control_data <- meta %>% filter(
    disease                 == "healthy",
    body_site               == SAMPLE_LOCATION,
    antibiotics_current_use %in% acceptable_abx,
    gender                  %in% acceptable_gender,
    BMI                     >= MIN_BMI,
    BMI                     <= MAX_BMI,
    minimum_read_length     >= MIN_READS
  )
  
  cat("\nCohort found:", nrow(disease_data), "with", DIAGNOSIS,
      "and", nrow(control_data), "healthy controls.\n")
  
  if (menu(c("Yes", "No"), title = "Proceed?") == 2) stop("Analysis halted by user.")
  
  
  # -- 3. Project Folder Setup --------------------------------------------------
  
  PROJECT_ID  <- paste0(date, "_", DIAGNOSIS, "_", SAMPLE_LOCATION, "_",
                        abx_label, "_", gen_selected,
                        "_BMIrange-", MIN_BMI, "-", MAX_BMI,
                        "_minreads-", MIN_READS)
  FOLDER_PATH <- PROJECT_ID
  
  subdirs <- c("data", "plots", "results")
  for (d in subdirs) {
    path <- file.path(FOLDER_PATH, d)
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  }
  message("Project folder ready: ", FOLDER_PATH)
  
  
  # -- 4. Age Matching -----------------------------------------------------------
  
  disease_age_min  <- min(disease_data$age, na.rm = TRUE)
  disease_age_max  <- max(disease_data$age, na.rm = TRUE)
  n_diseased       <- nrow(disease_data)
  
  set.seed(123)
  age_matched_controls <- control_data %>%
    filter(age >= disease_age_min, age <= disease_age_max) %>%
    slice_sample(n = n_diseased)
  
  meta_combined         <- rbind(disease_data, age_matched_controls)
  meta_combined$disease <- as.factor(meta_combined$disease)
  
  saveRDS(meta_combined, file.path(FOLDER_PATH, "data", "meta_combined.rds"))
  saveRDS(list(
    DIAGNOSIS         = DIAGNOSIS,
    SAMPLE_LOCATION   = SAMPLE_LOCATION,
    data_to_pull      = data_to_pull,
    acceptable_abx    = acceptable_abx,
    acceptable_gender = acceptable_gender,
    MIN_READS         = MIN_READS,
    MIN_BMI           = MIN_BMI,
    MAX_BMI           = MAX_BMI
  ), file.path(FOLDER_PATH, "data", "config.rds"))
  
  cat("Age-matched cohort:", nrow(meta_combined), "total samples",
      "(age range", disease_age_min, "–", disease_age_max, ")\n")
}

# =============================================================================
# FUNCTIONS
# =============================================================================

# -- pull_sample_data ---------------------------------------------------------
# Fetches a curatedMetagenomicData data type, builds a phyloseq object, and
# caches it as an .rds so subsequent runs load instantly.

pull_sample_data <- function(my_metadata, data_type) {
  gc()

  data_dir  <- file.path(FOLDER_PATH, "data")
  file_name <- file.path(data_dir, paste0(data_type, "_ps.rds"))

  if (file.exists(file_name)) {
    message("Loading from cache: ", file_name)
    ps_data <- readRDS(file_name)

  } else {
    message("Pulling ", data_type, ": this may take a few minutes...")
    se_object <- suppressMessages(returnSamples(my_metadata, dataType = data_type))

    # Filter rare gene families to save RAM
    if (data_type == "gene_families") {
      message("Filtering rare genes...")
      present_in_samples <- rowSums(assay(se_object) > 0)
      se_object          <- se_object[present_in_samples >= ncol(se_object) * 0.05, ]
    }

    counts_matrix <- assay(se_object)
    metadata_df   <- as.data.frame(colData(se_object))
    tax_raw       <- rownames(se_object)

    # Parse taxonomy if hierarchical
    if (any(grepl("k__", tax_raw))) {
      taxonomy_mat <- tax_raw %>%
        as.data.frame() %>%
        tidyr::separate(
          col  = 1,
          into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
          sep  = "\\|", fill = "right"
        ) %>%
        dplyr::mutate(across(everything(), ~ stringr::str_remove(.x, "^[kpcofgs]__"))) %>%
        as.matrix()
      rownames(taxonomy_mat) <- tax_raw
    } else {
      taxonomy_mat           <- as.matrix(tax_raw)
      colnames(taxonomy_mat) <- "Feature"
      rownames(taxonomy_mat) <- tax_raw
    }

    ps_data <- phyloseq(
      otu_table(as.matrix(counts_matrix), taxa_are_rows = TRUE),
      sample_data(metadata_df),
      tax_table(taxonomy_mat)
    )

    saveRDS(ps_data, file_name)
    message("Cached to: ", file_name)
  }

  return(list(
    ps         = ps_data,
    data_name  = data_type,
    clean_name = tools::toTitleCase(gsub("_", " ", data_type))
  ))
}


# -- save_csv()  --------------------------------------------------------------
# Thin wrapper so every CSV save is one line and self-documents its path.

save_csv <- function(df, file_stem) {
  path <- file.path(FOLDER_PATH, "results", paste0(file_stem, ".csv"))
  write.csv(df, path, row.names = FALSE)
  message("Saved: ", path)
}


# -- create_diversity_plots()  ------------------------------------------------
# Alpha (Shannon, Wilcoxon-annotated) + Beta (PCoA Bray-Curtis) side by side.
# Skips alpha for pathway_coverage (bounded metric, Shannon not meaningful).

create_diversity_plots <- function(ps_obj, file_type = ".png") {

  save_dir <- file.path(FOLDER_PATH, "plots")
  groups   <- as.character(unique(sample_data(ps_obj$ps)$disease))

  # Beta diversity — always produced
  ord    <- ordinate(ps_obj$ps, method = "PCoA", distance = "bray")
  p_beta <- plot_ordination(ps_obj$ps, ord, color = "disease") +
    stat_ellipse() +
    theme_minimal() +
    labs(title = paste("Beta Diversity (PCoA, Bray-Curtis) —", ps_obj$clean_name)) +
    theme(legend.title = element_blank())

  ggsave(file.path(save_dir, paste0(ps_obj$data_name, "_beta-diversity", file_type)),
         plot = p_beta, width = 12, height = 8, dpi = 300)

  # Alpha diversity - skipped for pathway_coverage
  if (ps_obj$data_name == "pathway_coverage") {
    message("Skipping alpha diversity for pathway_coverage (bounded metric).")
    return(list(beta_plot = p_beta, alpha_plot = NULL, combined_plot = p_beta))
  }

  richness_vals <- estimate_richness(ps_obj$ps, measures = "Shannon")
  max_height    <- max(richness_vals$Shannon, na.rm = TRUE)

  p_alpha <- plot_richness(ps_obj$ps, x = "disease", measures = "Shannon") +
    geom_boxplot(aes(fill = disease), alpha = 0.7, outlier.shape = NA) +
    theme_bw() +
    labs(title = "Alpha Diversity (Shannon)", x = "Diagnosis", y = "Shannon Index") +
    theme(legend.position = "none") +
    geom_signif(
      comparisons      = list(groups),
      map_signif_level = TRUE,
      test             = "wilcox.test",
      y_position       = max_height + 0.5
    )

  ggsave(file.path(save_dir, paste0(ps_obj$data_name, "_alpha-diversity", file_type)),
         plot = p_alpha, width = 12, height = 8, dpi = 300)

  # Combined panel
  title_grob <- ggdraw() +
    draw_label(
      paste0("Within and Between Group Differences — ", ps_obj$clean_name),
      fontface = "bold", size = 14, x = 0.5, hjust = 0.5
    )

  combined <- plot_grid(
    title_grob,
    plot_grid(p_alpha, p_beta, labels = c("A", "B"), ncol = 2, align = "h", axis = "bt"),
    ncol = 1, rel_heights = c(0.08, 1)
  )

  ggsave(file.path(save_dir, paste0(ps_obj$data_name, "_combined-diversity", file_type)),
         plot = combined, width = 14, height = 8, dpi = 300)

  message("Diversity plots saved for: ", ps_obj$clean_name)
  return(list(alpha_plot = p_alpha, 
              beta_plot = p_beta, 
              combined_plot = combined))
}


# -- run_permanova()  ---------------------------------------------------------
# PERMANOVA (adonis2) + betadisper dispersion test.

run_permanova <- function(ps_obj) {

  dist_matrix <- phyloseq::distance(ps_obj$ps, method = "bray")
  sample_df   <- as(sample_data(ps_obj$ps), "data.frame")

  permanova_result <- adonis2(dist_matrix ~ disease, data = sample_df)
  bd               <- betadisper(dist_matrix, sample_df$disease)
  bd_result        <- permutest(bd)

  # Flatten to data frames and save
  permanova_df  <- cbind(Term = rownames(as.data.frame(permanova_result)),
                         as.data.frame(permanova_result))
  betadisper_df <- cbind(Term = rownames(as.data.frame(bd_result$tab)),
                         as.data.frame(bd_result$tab))

  save_csv(permanova_df,  paste0(ps_obj$data_name, "_permanova"))
  save_csv(betadisper_df, paste0(ps_obj$data_name, "_betadisper"))

  return(list(
    permanova    = permanova_result,
    betadisper   = bd_result,
    ps_data_name = ps_obj$data_name
  ))
}


# -- create_top_group_plots()  ------------------------------------------------
# Faceted boxplots of the top N (5) features per group, with Wilcoxon significance bars.

create_top_group_plots <- function(ps_obj, group_var, rank_name, top_n = 5, file_type = ".png") {

  group_names <- unique(as.character(sample_data(ps_obj$ps)[[group_var]]))

  # Union of top N features from each group to capture group-specific signatures
  top_features <- c()
  for (g in group_names) {
    ps_sub       <- prune_samples(sample_data(ps_obj$ps)[[group_var]] == g, ps_obj$ps)
    top_features <- c(top_features, names(sort(taxa_sums(ps_sub), decreasing = TRUE))[1:top_n])
  }

  ps_top  <- prune_taxa(unique(top_features), ps_obj$ps)
  df_plot <- psmelt(ps_top)
  
  y_label <- switch(ps_obj$data_name,
                    "pathway_coverage"    = "Coverage (0–1)",
                    "pathway_abundance"   = "Pathway Abundance",
                    "relative_abundance"  = "Relative Abundance",
                    "gene_families"       = "Gene Families"
  )

  p <- ggplot(df_plot, aes(x = !!sym(group_var), y = Abundance, fill = !!sym(group_var))) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    facet_wrap(as.formula(paste("~", if (rank_name == "Feature") "OTU" else rank_name)), 
               scales = "free_y") +
    theme_bw() +
    labs(
      title = paste0("Top ", top_n, " ", rank_name, " — ", DIAGNOSIS, " vs Healthy Controls"),
      y     = y_label,
      x     = tools::toTitleCase(group_var)
    ) +
    theme(
      legend.position = "none",
      strip.text      = element_text(face = "bold.italic")
    ) +
    geom_signif(
      comparisons      = list(group_names),
      map_signif_level = TRUE,
      test             = "wilcox.test"
    ) + 
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

  save_path <- file.path(FOLDER_PATH, "plots",
                         paste0(ps_obj$data_name, "_top-", top_n, "-", rank_name, file_type))
  ggsave(save_path, plot = p, width = 12, height = 10, dpi = 300)
  message("Saved: ", save_path)

  return(list(top_plot = p, 
              ps_data_name = ps_obj$data_name, 
              top_df = df_plot))
}


# -- run_wilcox_taxa() --------------------------------------------------------
# Wilcoxon rank-sum test on the top N features, with FDR correction.

run_wilcox_taxa <- function(ps_obj, rank_level, top_n = 5) {
  feature_col <- if (rank_level == "Feature") "OTU" else rank_level
  
  # Build melted data frame — prune before melt for pathways to save memory
  if (rank_level == "Feature") {
    ps_glom  <- ps_obj$ps
    # Get union of top N per group BEFORE melting to keep memory manageable
    group_names    <- unique(as.character(sample_data(ps_glom)$disease))
    top_features   <- c()
    for (g in group_names) {
      ps_sub       <- prune_samples(sample_data(ps_glom)$disease == g, ps_glom)
      top_features <- c(top_features,
                        names(sort(taxa_sums(ps_sub), decreasing = TRUE))[1:top_n])
    }
    top_features <- unique(top_features)
    ps_glom      <- prune_taxa(top_features, ps_glom)
    df_melt      <- psmelt(ps_glom)
    top_labels   <- top_features
    
  } else {
    ps_glom      <- tax_glom(ps_obj$ps, rank_level)
    df_melt      <- psmelt(ps_glom)
    group_names  <- unique(as.character(sample_data(ps_glom)$disease))
    top_features <- c()
    for (g in group_names) {
      ps_sub       <- prune_samples(sample_data(ps_glom)$disease == g, ps_glom)
      top_features <- c(top_features,
                        names(sort(taxa_sums(ps_sub), decreasing = TRUE))[1:top_n])
    }
    top_labels <- unique(top_features)
    
    # Map OTU names to rank labels (e.g. OTU ID -> Genus name)
    tax_df     <- as.data.frame(tax_table(ps_glom))
    top_labels <- tax_df[top_labels, rank_level]
  }
  
  results <- lapply(top_labels, function(feat) {
    test_data <- df_melt[df_melt[[feature_col]] == feat, ]
    res       <- wilcox.test(Abundance ~ disease, data = test_data, exact = FALSE)
    data.frame(Feature = feat, p_value = res$p.value)
  })
  
  sig_table            <- do.call(rbind, results) %>% arrange(p_value)
  sig_table$p_adjusted <- p.adjust(sig_table$p_value, method = "BH")
  save_csv(sig_table, paste0(ps_obj$data_name, "_wilcox_", tolower(rank_level)))
  
  return(list(
    results      = sig_table,
    df_melt      = df_melt,
    ps_data_name = ps_obj$data_name,
    rank_level   = rank_level
  ))
}


# -- build_direction_table() --------------------------------------------------
# Merges Wilcoxon results with per-group mean abundances and a direction label
# (increased / decreased in disease group vs healthy).

build_direction_table <- function(wilcox_results, group_col = "disease",
                                  disease_label = DIAGNOSIS) {
  sig_table  <- wilcox_results$results
  df_melt    <- wilcox_results$df_melt
  rank_level <- wilcox_results$rank_level

  
  feature_col <- if (rank_level == "Feature") "OTU" else rank_level
  
  summary_table <- df_melt %>%
    filter(.data[[feature_col]] %in% sig_table$Feature) %>%
    group_by(.data[[feature_col]], .data[[group_col]]) %>%
    summarize(mean_abundance = mean(Abundance), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = all_of(group_col), values_from = mean_abundance)
  
  colnames(summary_table)[1] <- "Feature"
  summary_table <- as.data.frame(summary_table)
  
  summary_table <- summary_table %>%
    mutate(direction = ifelse(
      .data[[disease_label]] > .data[["healthy"]],
      paste0("Increased in ", disease_label),
      paste0("Decreased in ", disease_label)
    ))

  final_results <- merge(sig_table, summary_table, by = "Feature")
  save_csv(final_results, paste0(wilcox_results$ps_data_name, "_direction_table"))

  return(final_results)
}


# =============================================================================
# MAIN — Data Pull & Analysis Loop
# =============================================================================

# Pull all selected data types (uses cache if already downloaded)
ps_list        <- lapply(data_to_pull, pull_sample_data, my_metadata = meta_combined)
names(ps_list) <- unlist(data_to_pull)

all_results <- list()

for (dtype in names(ps_list)) {
  # Format of messages provided by Claude
  ps_obj <- ps_list[[dtype]]
  cat("\n\n══════════════════════════════════════════════════\n")
  cat("  Analysing:", ps_obj$clean_name, "\n")
  cat("══════════════════════════════════════════════════\n")

  results <- list(data_type = dtype)

  # 1. Diversity plots
  results$diversity <- create_diversity_plots(ps_obj)

  # 2. PERMANOVA + betadisper
  results$permanova <- run_permanova(ps_obj)

  # 3. Top features + Wilcoxon + FDR + direction table
  if (dtype == "relative_abundance") {
    # Agglomerate to genus level for relative abundance
    ps_genus            <- ps_obj
    ps_genus$ps         <- tax_glom(ps_obj$ps, "Genus")
    ps_genus$data_name  <- paste0(dtype, "_genus")
    ps_genus$clean_name <- "Relative Abundance (Genus)"

    results$top_plot  <- create_top_group_plots(ps_genus, "disease", "Genus", top_n = 5)
    results$wilcox    <- run_wilcox_taxa(ps_genus, "Genus", top_n = 5)
    results$direction <- build_direction_table(results$wilcox)

  } else {
    # Pathway abundance & coverage use "Feature" as the rank label
    results$top_plot  <- create_top_group_plots(ps_obj, "disease", "Feature", top_n = 5)
    results$wilcox    <- run_wilcox_taxa(ps_obj, "Feature", top_n = 5)
    results$direction <- build_direction_table(results$wilcox)
  }

  all_results[[dtype]] <- results
}

cat("\n\n✓ Analysis complete. All outputs written to:", FOLDER_PATH, "\n")
