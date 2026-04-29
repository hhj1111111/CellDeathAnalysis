# =============================================================================
# TCGA Pan-Cancer Analysis v2 - With Crosstalk Scoring & Subtype Classification
# =============================================================================
#
# Enhanced analysis using CellDeathAnalysis v0.4.0 new features:
# 1. Crosstalk-aware scoring (novel algorithm)
# 2. Cell death subtype classification (novel algorithm)
# 3. Comparison between z-score and crosstalk methods
#
# Author: Keran Sun
# Email: s1214844197@163.com
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Environment Setup
# -----------------------------------------------------------------------------

packages_cran <- c("tidyverse", "survival", "survminer", "ggpubr",
                   "pheatmap", "RColorBrewer", "ggrepel", "patchwork",
                   "corrplot", "reshape2")

packages_bioc <- c("TCGAbiolinks", "SummarizedExperiment", "DESeq2",
                   "ComplexHeatmap", "circlize")

for (pkg in packages_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in packages_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
}

library(tidyverse)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(patchwork)
library(CellDeathAnalysis)

dir.create("results", showWarnings = FALSE)
dir.create("results/figures", showWarnings = FALSE)
dir.create("results/tables", showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Download TCGA Pan-Cancer Data (same as v1)
# -----------------------------------------------------------------------------

cancer_types <- c(
  "TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL",
  "TCGA-COAD", "TCGA-DLBC", "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC",
  "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LAML", "TCGA-LGG",
  "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", "TCGA-OV",
  "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC",
  "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM",
  "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM"
)

download_tcga_data <- function(project, save_dir = "data/tcga") {
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  save_file <- file.path(save_dir, paste0(project, "_data.rds"))

  if (file.exists(save_file)) {
    message("Loading cached data for ", project)
    return(readRDS(save_file))
  }

  message("Downloading ", project, "...")

  tryCatch({
    query <- GDCquery(
      project = project,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts"
    )

    GDCdownload(query, method = "api", files.per.chunk = 50)
    data <- GDCprepare(query, summarizedExperiment = TRUE)

    if ("tpm_unstrand" %in% assayNames(data)) {
      expr <- assay(data, "tpm_unstrand")
    } else if ("unstranded" %in% assayNames(data)) {
      expr <- assay(data, "unstranded")
    } else {
      expr <- assay(data, 1)
    }

    gene_info <- rowData(data)
    if ("gene_name" %in% colnames(gene_info)) {
      rownames(expr) <- gene_info$gene_name
    }

    dup_genes <- duplicated(rownames(expr))
    expr <- expr[!dup_genes, ]

    sample_info <- as.data.frame(colData(data))

    if ("sample_type" %in% colnames(sample_info)) {
      sample_info$sample_type_simple <- ifelse(
        grepl("Normal", sample_info$sample_type, ignore.case = TRUE),
        "Normal", "Tumor"
      )
    } else if ("shortLetterCode" %in% colnames(sample_info)) {
      sample_info$sample_type_simple <- ifelse(
        sample_info$shortLetterCode == "NT", "Normal", "Tumor"
      )
    } else {
      sample_code <- substr(colnames(expr), 14, 15)
      sample_info$sample_type_simple <- ifelse(
        as.numeric(sample_code) < 10, "Tumor", "Normal"
      )
    }

    clinical <- tryCatch({
      GDCquery_clinic(project, type = "clinical")
    }, error = function(e) {
      message("  Warning: Could not fetch clinical data: ", e$message)
      NULL
    })

    result <- list(expr = expr, sample_info = sample_info,
                   clinical = clinical, project = project)
    saveRDS(result, save_file)
    return(result)

  }, error = function(e) {
    message("Error downloading ", project, ": ", e$message)
    return(NULL)
  })
}

# Download all cancer types
tcga_data_list <- list()
for (cancer in cancer_types) {
  result <- download_tcga_data(cancer)
  if (!is.null(result)) tcga_data_list[[cancer]] <- result
}

tcga_data_list <- tcga_data_list[!sapply(tcga_data_list, is.null)]
saveRDS(tcga_data_list, "results/tcga_data_list.rds")
message("\nSuccessfully downloaded ", length(tcga_data_list), " cancer types")

# -----------------------------------------------------------------------------
# 2. Calculate Cell Death Scores - BOTH Methods
# -----------------------------------------------------------------------------

calculate_cancer_scores_both <- function(tcga_data) {
  if (is.null(tcga_data)) return(NULL)

  expr <- log2(tcga_data$expr + 1)

  # Method 1: z-score (original)
  scores_zscore <- calculate_death_score(
    expr, pathways = "all", method = "zscore",
    min_genes = 5, verbose = FALSE
  )

  # Method 2: crosstalk-aware (novel)
  scores_crosstalk <- calculate_death_score(
    expr, pathways = "all", method = "crosstalk",
    min_genes = 5, verbose = FALSE
  )

  # Add sample info to both
  add_info <- function(scores) {
    scores$sample_id <- rownames(scores)
    if ("barcode" %in% colnames(tcga_data$sample_info)) {
      match_idx <- match(scores$sample_id, tcga_data$sample_info$barcode)
    } else {
      match_idx <- match(scores$sample_id, rownames(tcga_data$sample_info))
    }
    scores$sample_type <- tcga_data$sample_info$sample_type_simple[match_idx]
    scores$project <- tcga_data$project
    return(scores)
  }

  list(
    zscore = add_info(scores_zscore),
    crosstalk = add_info(scores_crosstalk)
  )
}

# Calculate scores for all cancers
all_scores_zscore <- list()
all_scores_crosstalk <- list()

for (cancer in names(tcga_data_list)) {
  message("Calculating scores for ", cancer)
  result <- calculate_cancer_scores_both(tcga_data_list[[cancer]])
  if (!is.null(result)) {
    all_scores_zscore[[cancer]] <- result$zscore
    all_scores_crosstalk[[cancer]] <- result$crosstalk
  }
}

# Combine scores
combine_scores <- function(score_list) {
  combined <- do.call(rbind, score_list)
  rownames(combined) <- NULL
  combined[!is.na(combined$sample_type), ]
}

pan_cancer_zscore <- combine_scores(all_scores_zscore)
pan_cancer_crosstalk <- combine_scores(all_scores_crosstalk)

saveRDS(pan_cancer_zscore, "results/pan_cancer_zscore.rds")
saveRDS(pan_cancer_crosstalk, "results/pan_cancer_crosstalk.rds")

message("Total samples: ", nrow(pan_cancer_zscore))

# -----------------------------------------------------------------------------
# 3. Compare Z-score vs Crosstalk Methods (NEW - Figure 3D)
# -----------------------------------------------------------------------------

pathway_cols <- names(get_death_geneset("all", type = "simple"))

# Calculate correlation between methods per cancer
method_comparison <- list()

for (cancer in unique(pan_cancer_zscore$project)) {
  zscore_sub <- pan_cancer_zscore[pan_cancer_zscore$project == cancer, ]
  crosstalk_sub <- pan_cancer_crosstalk[pan_cancer_crosstalk$project == cancer, ]

  # Match samples
  common_samples <- intersect(zscore_sub$sample_id, crosstalk_sub$sample_id)
  if (length(common_samples) < 10) next

  zscore_matched <- zscore_sub[match(common_samples, zscore_sub$sample_id), ]
  crosstalk_matched <- crosstalk_sub[match(common_samples, crosstalk_sub$sample_id), ]

  for (pw in pathway_cols) {
    if (!pw %in% colnames(zscore_matched) || !pw %in% colnames(crosstalk_matched)) next

    r_val <- cor(zscore_matched[[pw]], crosstalk_matched[[pw]],
                 use = "complete.obs", method = "spearman")

    method_comparison[[paste(cancer, pw, sep = "_")]] <- data.frame(
      cancer = cancer,
      pathway = pw,
      spearman_r = r_val,
      stringsAsFactors = FALSE
    )
  }
}

method_comp_df <- do.call(rbind, method_comparison)

# Summary: mean correlation per pathway
method_comp_summary <- method_comp_df %>%
  group_by(pathway) %>%
  summarise(
    mean_r = mean(spearman_r, na.rm = TRUE),
    sd_r = sd(spearman_r, na.rm = TRUE),
    min_r = min(spearman_r, na.rm = TRUE),
    max_r = max(spearman_r, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_r)

write.csv(method_comp_summary, "results/tables/method_comparison_summary.csv", row.names = FALSE)

# Heatmap of method comparison
comp_matrix <- method_comp_df %>%
  select(cancer, pathway, spearman_r) %>%
  pivot_wider(names_from = pathway, values_from = spearman_r) %>%
  column_to_rownames("cancer") %>%
  as.matrix()

pdf("results/figures/Figure3D_method_comparison.pdf", width = 12, height = 10)

col_fun_comp <- colorRamp2(c(0.8, 0.95, 1.0), c("white", "#FDBF6F", "#E31A1C"))

Heatmap(
  comp_matrix,
  name = "Spearman r",
  col = col_fun_comp,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_labels = gsub("TCGA-", "", rownames(comp_matrix)),
  column_names_rot = 45,
  column_title = "Z-score vs Crosstalk Method Comparison",
  row_title = "Cancer Type",
  heatmap_legend_param = list(title = "Spearman\nCorrelation")
)

dev.off()

# -----------------------------------------------------------------------------
# 4. Pan-Cancer Heatmap with Crosstalk Scores (Figure 3A updated)
# -----------------------------------------------------------------------------

mean_scores_crosstalk <- pan_cancer_crosstalk %>%
  group_by(project, sample_type) %>%
  summarise(across(all_of(pathway_cols), mean, na.rm = TRUE), .groups = "drop")

tumor_matrix_crosstalk <- mean_scores_crosstalk %>%
  filter(sample_type == "Tumor") %>%
  select(-sample_type) %>%
  column_to_rownames("project") %>%
  as.matrix()

tumor_matrix_scaled_crosstalk <- scale(tumor_matrix_crosstalk)

col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

pdf("results/figures/Figure3A_pancancer_heatmap_crosstalk.pdf", width = 12, height = 10)

ht <- Heatmap(
  tumor_matrix_scaled_crosstalk,
  name = "Z-score",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  column_names_rot = 45,
  row_labels = gsub("TCGA-", "", rownames(tumor_matrix_scaled_crosstalk)),
  column_title = "Cell Death Pathways (Crosstalk-Aware)",
  row_title = "TCGA Cancer Types"
)

draw(ht)
dev.off()

# -----------------------------------------------------------------------------
# 5. Tumor vs Normal Comparison - Both Methods (Figure 3B)
# -----------------------------------------------------------------------------

compare_tumor_normal <- function(scores_df, pathway) {
  tumor <- scores_df[[pathway]][scores_df$sample_type == "Tumor"]
  normal <- scores_df[[pathway]][scores_df$sample_type == "Normal"]

  if (length(normal) < 3) {
    return(data.frame(
      pathway = pathway,
      mean_tumor = mean(tumor, na.rm = TRUE),
      mean_normal = NA,
      log2FC = NA,
      p_value = NA,
      stringsAsFactors = FALSE
    ))
  }

  test <- wilcox.test(tumor, normal)

  data.frame(
    pathway = pathway,
    mean_tumor = mean(tumor, na.rm = TRUE),
    mean_normal = mean(normal, na.rm = TRUE),
    log2FC = log2(mean(tumor, na.rm = TRUE) / mean(normal, na.rm = TRUE)),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}

# Compare for both methods
run_comparison <- function(all_scores, method_name) {
  comparison_results <- list()
  for (cancer in names(tcga_data_list)) {
    cancer_scores <- all_scores[[cancer]]
    if (sum(cancer_scores$sample_type == "Normal") >= 3) {
      results <- lapply(pathway_cols, function(p) compare_tumor_normal(cancer_scores, p))
      results_df <- do.call(rbind, results)
      results_df$cancer <- cancer
      comparison_results[[cancer]] <- results_df
    }
  }
  comparison_df <- do.call(rbind, comparison_results)
  comparison_df$p_adjust <- p.adjust(comparison_df$p_value, method = "BH")
  comparison_df$method <- method_name
  return(comparison_df)
}

comp_zscore <- run_comparison(all_scores_zscore, "zscore")
comp_crosstalk <- run_comparison(all_scores_crosstalk, "crosstalk")

# Combine and save
comp_combined <- rbind(comp_zscore, comp_crosstalk)
write.csv(comp_combined, "results/tables/tumor_vs_normal_both_methods.csv", row.names = FALSE)

# Heatmap for crosstalk method
fc_matrix_crosstalk <- comp_crosstalk %>%
  select(cancer, pathway, log2FC) %>%
  pivot_wider(names_from = pathway, values_from = log2FC) %>%
  column_to_rownames("cancer") %>%
  as.matrix()

fc_matrix_crosstalk[is.na(fc_matrix_crosstalk)] <- 0
fc_matrix_crosstalk[is.infinite(fc_matrix_crosstalk)] <- 0

pdf("results/figures/Figure3B_tumor_normal_FC_crosstalk.pdf", width = 12, height = 8)

col_fun_fc <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

Heatmap(
  fc_matrix_crosstalk,
  name = "log2FC",
  col = col_fun_fc,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_labels = gsub("TCGA-", "", rownames(fc_matrix_crosstalk)),
  column_names_rot = 45,
  column_title = "Tumor vs Normal (Crosstalk-Aware, log2FC)",
  row_title = "Cancer Type"
)

dev.off()

# -----------------------------------------------------------------------------
# 6. Cell Death Subtype Classification (NEW - Figure 6)
# -----------------------------------------------------------------------------

message("\n========================================")
message("  Cell Death Subtype Classification")
message("========================================\n")

# For each cancer type with enough samples, classify subtypes
subtype_results <- list()

for (cancer in unique(pan_cancer_crosstalk$project)) {
  cancer_data <- pan_cancer_crosstalk[pan_cancer_crosstalk$project == cancer, ]
  tumor_data <- cancer_data[cancer_data$sample_type == "Tumor", ]

  if (nrow(tumor_data) < 30) {
    message("Skipping ", cancer, " - too few tumor samples (", nrow(tumor_data), ")")
    next
  }

  message("Classifying subtypes for ", cancer, " (n=", nrow(tumor_data), ")")

  # Get scores matrix
  scores_mat <- tumor_data[, pathway_cols, drop = FALSE]
  rownames(scores_mat) <- tumor_data$sample_id

  # Get survival data
  clinical <- tcga_data_list[[cancer]]$clinical
  time_vec <- NULL
  status_vec <- NULL

  if (!is.null(clinical) && "submitter_id" %in% colnames(clinical)) {
    patient_ids <- substr(tumor_data$sample_id, 1, 12)
    clinical_match <- clinical[match(patient_ids, clinical$submitter_id), ]

    if ("days_to_death" %in% colnames(clinical)) {
      time_vec <- ifelse(
        !is.na(clinical_match$days_to_death),
        as.numeric(clinical_match$days_to_death),
        as.numeric(clinical_match$days_to_last_follow_up)
      )
    }

    if ("vital_status" %in% colnames(clinical)) {
      status_vec <- ifelse(
        tolower(clinical_match$vital_status) %in% c("dead", "deceased"), 1, 0
      )
    }
  }

  # Run subtype classification
  tryCatch({
    subtype_result <- classify_death_subtypes(
      scores_mat,
      time = time_vec,
      status = status_vec,
      max_k = 6,
      n_reps = 500,
      seed = 42,
      verbose = FALSE
    )

    subtype_results[[cancer]] <- subtype_result

    # Save individual result
    saveRDS(subtype_result, paste0("results/tables/", cancer, "_subtypes.rds"))

    message("  Optimal k = ", subtype_result$optimal_k)
    message("  Subtypes: ", paste(subtype_result$subtype_names, collapse = ", "))

  }, error = function(e) {
    message("  Error: ", e$message)
  })
}

saveRDS(subtype_results, "results/subtype_results_all.rds")

# -----------------------------------------------------------------------------
# 7. Subtype Summary Statistics (NEW - Table 3)
# -----------------------------------------------------------------------------

subtype_summary <- list()

for (cancer in names(subtype_results)) {
  res <- subtype_results[[cancer]]

  for (sname in names(res$subtype_names)) {
    subtype_summary[[paste(cancer, sname)]] <- data.frame(
      cancer = cancer,
      cluster = sname,
      subtype_name = res$subtype_names[sname],
      n_samples = sum(res$subtypes == sname),
      total_samples = length(res$subtypes),
      pct = 100 * sum(res$subtypes == sname) / length(res$subtypes),
      stringsAsFactors = FALSE
    )
  }
}

subtype_summary_df <- do.call(rbind, subtype_summary)
write.csv(subtype_summary_df, "results/tables/subtype_summary.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 8. Survival Analysis with Crosstalk Scores (Figure 4 updated)
# -----------------------------------------------------------------------------

perform_survival_analysis <- function(tcga_data, scores) {
  if (is.null(tcga_data) || is.null(tcga_data$clinical)) return(NULL)

  clinical <- tcga_data$clinical
  if (!"submitter_id" %in% colnames(clinical)) return(NULL)

  scores_tumor <- scores[scores$sample_type == "Tumor", ]
  if (nrow(scores_tumor) < 20) return(NULL)

  scores_tumor$patient_id <- substr(scores_tumor$sample_id, 1, 12)

  clinical$OS_time <- NA
  clinical$OS_status <- NA

  if ("days_to_death" %in% colnames(clinical)) {
    clinical$OS_time <- ifelse(
      !is.na(clinical$days_to_death),
      as.numeric(clinical$days_to_death),
      as.numeric(clinical$days_to_last_follow_up)
    )
  } else if ("days_to_last_follow_up" %in% colnames(clinical)) {
    clinical$OS_time <- as.numeric(clinical$days_to_last_follow_up)
  }

  if ("vital_status" %in% colnames(clinical)) {
    clinical$OS_status <- ifelse(
      tolower(clinical$vital_status) %in% c("dead", "deceased"), 1, 0
    )
  }

  matched <- merge(
    scores_tumor,
    clinical[, c("submitter_id", "OS_time", "OS_status")],
    by.x = "patient_id", by.y = "submitter_id", all.x = FALSE
  )

  matched <- matched %>%
    filter(!is.na(OS_time) & !is.na(OS_status) & OS_time > 0)

  if (nrow(matched) < 30) return(NULL)
  return(matched)
}

# Survival analysis for crosstalk scores
survival_results_crosstalk <- list()

for (cancer in names(tcga_data_list)) {
  message("Survival analysis (crosstalk) for ", cancer)

  if (!cancer %in% names(all_scores_crosstalk)) next

  matched <- perform_survival_analysis(
    tcga_data_list[[cancer]], all_scores_crosstalk[[cancer]]
  )

  if (is.null(matched)) next

  for (pathway in pathway_cols) {
    if (!pathway %in% colnames(matched)) next
    if (all(is.na(matched[[pathway]]))) next

    tryCatch({
      pathway_values <- matched[[pathway]]
      med_val <- median(pathway_values, na.rm = TRUE)
      matched$group <- ifelse(pathway_values > med_val, "High", "Low")
      matched$group <- factor(matched$group, levels = c("Low", "High"))

      group_table <- table(matched$group)
      if (any(group_table < 10)) next

      cox_data <- data.frame(
        time = matched$OS_time,
        status = matched$OS_status,
        score = pathway_values
      )
      cox_data <- cox_data[complete.cases(cox_data), ]

      cox_fit <- coxph(Surv(time, status) ~ score, data = cox_data)
      cox_sum <- summary(cox_fit)

      surv_data <- data.frame(
        time = matched$OS_time,
        status = matched$OS_status,
        group = matched$group
      )
      surv_data <- surv_data[complete.cases(surv_data), ]

      surv_diff <- survdiff(Surv(time, status) ~ group, data = surv_data)
      log_rank_p <- 1 - pchisq(surv_diff$chisq, df = 1)

      survival_results_crosstalk[[paste(cancer, pathway, sep = "_")]] <- data.frame(
        cancer = cancer,
        pathway = pathway,
        method = "crosstalk",
        n = nrow(cox_data),
        hr = cox_sum$conf.int[1, 1],
        hr_lower = cox_sum$conf.int[1, 3],
        hr_upper = cox_sum$conf.int[1, 4],
        cox_p = cox_sum$coefficients[1, 5],
        log_rank_p = log_rank_p,
        stringsAsFactors = FALSE
      )

    }, error = function(e) {})
  }
}

# Also run for zscore for comparison
survival_results_zscore <- list()

for (cancer in names(tcga_data_list)) {
  if (!cancer %in% names(all_scores_zscore)) next

  matched <- perform_survival_analysis(
    tcga_data_list[[cancer]], all_scores_zscore[[cancer]]
  )

  if (is.null(matched)) next

  for (pathway in pathway_cols) {
    if (!pathway %in% colnames(matched)) next
    if (all(is.na(matched[[pathway]]))) next

    tryCatch({
      pathway_values <- matched[[pathway]]
      med_val <- median(pathway_values, na.rm = TRUE)
      matched$group <- ifelse(pathway_values > med_val, "High", "Low")
      matched$group <- factor(matched$group, levels = c("Low", "High"))

      group_table <- table(matched$group)
      if (any(group_table < 10)) next

      cox_data <- data.frame(
        time = matched$OS_time,
        status = matched$OS_status,
        score = pathway_values
      )
      cox_data <- cox_data[complete.cases(cox_data), ]

      cox_fit <- coxph(Surv(time, status) ~ score, data = cox_data)
      cox_sum <- summary(cox_fit)

      surv_data <- data.frame(
        time = matched$OS_time,
        status = matched$OS_status,
        group = matched$group
      )
      surv_data <- surv_data[complete.cases(surv_data), ]

      surv_diff <- survdiff(Surv(time, status) ~ group, data = surv_data)
      log_rank_p <- 1 - pchisq(surv_diff$chisq, df = 1)

      survival_results_zscore[[paste(cancer, pathway, sep = "_")]] <- data.frame(
        cancer = cancer,
        pathway = pathway,
        method = "zscore",
        n = nrow(cox_data),
        hr = cox_sum$conf.int[1, 1],
        hr_lower = cox_sum$conf.int[1, 3],
        hr_upper = cox_sum$conf.int[1, 4],
        cox_p = cox_sum$coefficients[1, 5],
        log_rank_p = log_rank_p,
        stringsAsFactors = FALSE
      )

    }, error = function(e) {})
  }
}

# Combine survival results
survival_df_crosstalk <- if (length(survival_results_crosstalk) > 0) {
  df <- do.call(rbind, survival_results_crosstalk)
  df$cox_p_adj <- p.adjust(df$cox_p, method = "BH")
  df$log_rank_p_adj <- p.adjust(df$log_rank_p, method = "BH")
  df
} else data.frame()

survival_df_zscore <- if (length(survival_results_zscore) > 0) {
  df <- do.call(rbind, survival_results_zscore)
  df$cox_p_adj <- p.adjust(df$cox_p, method = "BH")
  df$log_rank_p_adj <- p.adjust(df$log_rank_p, method = "BH")
  df
} else data.frame()

survival_combined <- rbind(survival_df_crosstalk, survival_df_zscore)
write.csv(survival_combined, "results/tables/survival_analysis_both_methods.csv", row.names = FALSE)

# Forest plot for ferroptosis (crosstalk)
if (nrow(survival_df_crosstalk) > 0) {
  ferro_survival <- survival_df_crosstalk %>%
    filter(pathway == "ferroptosis") %>%
    arrange(hr)

  if (nrow(ferro_survival) > 0) {
    ferro_survival$cancer_abbr <- gsub("TCGA-", "", ferro_survival$cancer)
    ferro_survival$cancer_abbr <- factor(ferro_survival$cancer_abbr,
                                          levels = ferro_survival$cancer_abbr)

    pdf("results/figures/Figure4A_ferroptosis_forest_crosstalk.pdf", width = 10, height = 12)

    p_forest <- ggplot(ferro_survival, aes(x = hr, y = cancer_abbr)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
      geom_errorbarh(aes(xmin = hr_lower, xmax = hr_upper), height = 0.2) +
      geom_point(aes(color = cox_p_adj < 0.05), size = 3) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                         labels = c("TRUE" = "p < 0.05", "FALSE" = "NS"),
                         name = "Significance") +
      scale_x_log10() +
      theme_bw() +
      theme(legend.position = "bottom", axis.text.y = element_text(size = 10)) +
      labs(x = "Hazard Ratio (95% CI)", y = NULL,
           title = "Ferroptosis Score and Overall Survival (Crosstalk-Aware)")

    print(p_forest)
    dev.off()
  }
}

# -----------------------------------------------------------------------------
# 9. Subtype Survival Analysis (NEW - Figure 6B)
# -----------------------------------------------------------------------------

for (cancer in names(subtype_results)) {
  res <- subtype_results[[cancer]]

  if (!is.null(res$survival_result) && !is.null(res$survival_result$fit)) {
    pdf(paste0("results/figures/Figure6B_subtype_survival_", gsub("TCGA-", "", cancer), ".pdf"),
        width = 8, height = 6)

    p <- plot_subtype_survival(res)
    print(p)

    dev.off()
    message("Subtype survival plot saved for ", cancer)
  }
}

# -----------------------------------------------------------------------------
# 10. Subtype Heatmap (NEW - Figure 6A)
# -----------------------------------------------------------------------------

for (cancer in names(subtype_results)) {
  res <- subtype_results[[cancer]]
  cancer_data <- pan_cancer_crosstalk[pan_cancer_crosstalk$project == cancer, ]
  tumor_data <- cancer_data[cancer_data$sample_type == "Tumor", ]
  scores_mat <- tumor_data[, pathway_cols, drop = FALSE]
  rownames(scores_mat) <- tumor_data$sample_id

  pdf(paste0("results/figures/Figure6A_subtype_heatmap_", gsub("TCGA-", "", cancer), ".pdf"),
      width = 10, height = 8)

  tryCatch({
    p <- plot_subtype_heatmap(res, scores_mat)
    if (!is.null(p)) print(p)
  }, error = function(e) {
    message("  Heatmap error for ", cancer, ": ", e$message)
  })

  dev.off()
}

# -----------------------------------------------------------------------------
# 11. Pathway Correlation Analysis (Figure 5)
# -----------------------------------------------------------------------------

cor_matrix_crosstalk <- cor(
  pan_cancer_crosstalk[, pathway_cols],
  use = "pairwise.complete.obs",
  method = "spearman"
)

pdf("results/figures/Figure5_pathway_correlation_crosstalk.pdf", width = 10, height = 10)

col_fun_cor <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

Heatmap(
  cor_matrix_crosstalk,
  name = "Correlation",
  col = col_fun_cor,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_names_rot = 45,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", cor_matrix_crosstalk[i, j]), x, y, gp = gpar(fontsize = 8))
  },
  column_title = "Cell Death Pathway Correlation (Crosstalk-Aware, Spearman)"
)

dev.off()

# Compare correlation matrices
cor_matrix_zscore <- cor(
  pan_cancer_zscore[, pathway_cols],
  use = "pairwise.complete.obs",
  method = "spearman"
)

# Correlation difference (should show reduced cross-pathway correlation)
cor_diff <- cor_matrix_zscore - cor_matrix_crosstalk

pdf("results/figures/Figure5C_correlation_difference.pdf", width = 10, height = 10)

col_fun_diff <- colorRamp2(c(0, 0.1, 0.3), c("white", "#FDBF6F", "#E31A1C"))

Heatmap(
  cor_diff,
  name = "r difference",
  col = col_fun_diff,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_names_rot = 45,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.3f", cor_diff[i, j]), x, y, gp = gpar(fontsize = 8))
  },
  column_title = "Correlation Reduction: Z-score minus Crosstalk"
)

dev.off()

# -----------------------------------------------------------------------------
# 12. Summary Report
# -----------------------------------------------------------------------------

summary_stats <- list(
  n_cancers = length(unique(pan_cancer_zscore$project)),
  n_samples_total = nrow(pan_cancer_zscore),
  n_tumor = sum(pan_cancer_zscore$sample_type == "Tumor", na.rm = TRUE),
  n_normal = sum(pan_cancer_zscore$sample_type == "Normal", na.rm = TRUE),
  n_pathways = length(pathway_cols),
  n_genes_total = length(get_all_death_genes()),
  n_subtypes_classified = length(subtype_results),
  mean_method_correlation = mean(method_comp_summary$mean_r),
  n_significant_survival_crosstalk = if(nrow(survival_df_crosstalk) > 0)
    sum(survival_df_crosstalk$cox_p_adj < 0.05, na.rm = TRUE) else 0,
  n_significant_survival_zscore = if(nrow(survival_df_zscore) > 0)
    sum(survival_df_zscore$cox_p_adj < 0.05, na.rm = TRUE) else 0
)

saveRDS(summary_stats, "results/summary_statistics_v2.rds")

cat("\n")
cat("=============================================================\n")
cat("    TCGA Pan-Cancer Analysis v2 Summary                      \n")
cat("=============================================================\n")
cat("\n")
cat("Data Summary:\n")
cat("  Cancer types analyzed:", summary_stats$n_cancers, "\n")
cat("  Total samples:", summary_stats$n_samples_total, "\n")
cat("  Tumor samples:", summary_stats$n_tumor, "\n")
cat("  Normal samples:", summary_stats$n_normal, "\n")
cat("\n")
cat("Gene Sets:\n")
cat("  Cell death pathways:", summary_stats$n_pathways, "\n")
cat("  Total genes:", summary_stats$n_genes_total, "\n")
cat("\n")
cat("Novel Algorithms:\n")
cat("  Subtypes classified:", summary_stats$n_subtypes_classified, "cancer types\n")
cat("  Mean method correlation:", round(summary_stats$mean_method_correlation, 3), "\n")
cat("\n")
cat("Survival Analysis:\n")
cat("  Significant (crosstalk):", summary_stats$n_significant_survival_crosstalk, "\n")
cat("  Significant (zscore):", summary_stats$n_significant_survival_zscore, "\n")
cat("\n")
cat("Output Files:\n")
cat("  - results/pan_cancer_zscore.rds\n")
cat("  - results/pan_cancer_crosstalk.rds\n")
cat("  - results/subtype_results_all.rds\n")
cat("  - results/tables/method_comparison_summary.csv\n")
cat("  - results/tables/tumor_vs_normal_both_methods.csv\n")
cat("  - results/tables/survival_analysis_both_methods.csv\n")
cat("  - results/tables/subtype_summary.csv\n")
cat("  - results/figures/Figure3A_pancancer_heatmap_crosstalk.pdf\n")
cat("  - results/figures/Figure3B_tumor_normal_FC_crosstalk.pdf\n")
cat("  - results/figures/Figure3D_method_comparison.pdf\n")
cat("  - results/figures/Figure4A_ferroptosis_forest_crosstalk.pdf\n")
cat("  - results/figures/Figure5_pathway_correlation_crosstalk.pdf\n")
cat("  - results/figures/Figure5C_correlation_difference.pdf\n")
cat("  - results/figures/Figure6A_subtype_heatmap_*.pdf\n")
cat("  - results/figures/Figure6B_subtype_survival_*.pdf\n")
cat("\n")
cat("=============================================================\n")

cat("\nAnalysis v2 completed successfully!\n")
