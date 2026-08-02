#' =============================================================================
#' Single-Cell Analysis Module for CellDeathAnalysis
#' =============================================================================

#' Calculate Death Scores for Single-Cell Data
#'
#' Calculate cell death pathway scores for single-cell RNA-seq data.
#' Supports Seurat and SingleCellExperiment objects.
#'
#' @param object A Seurat object, SingleCellExperiment object, or expression matrix.
#' @param pathways Character vector of pathways to analyze, or "all".
#' @param method Scoring method: "aucell" (recommended for scRNA-seq), "seurat",
#'   "mean", or "zscore".
#' @param assay For Seurat objects, which assay to use. Default is "RNA".
#' @param slot For Seurat objects, which slot to use. Default is "data".
#' @param ncores Number of cores for parallel processing. Default is 1.
#' @param add_to_object Logical. Whether to add scores to the object metadata.
#'   Default is TRUE.
#'
#' @return If add_to_object is TRUE, returns the modified object with scores
#'   added to metadata. Otherwise, returns a data frame of scores.
#'
#' @export
#'
#' @examples
#' \dontrun
#' # With Seurat object
#' seurat_obj <- sc_death_score(seurat_obj, method = "aucell")
#'
#' # View scores
#' head(seurat_obj@meta.data[, grep("death_", names(seurat_obj@meta.data))])
#'
#' # Visualize
#' FeaturePlot(seurat_obj, features = "death_ferroptosis")
#' }
#'
sc_death_score <- function(object,species = 'human',
                            pathways = "all",
                            method = c("aucell", "seurat", "mean", "zscore"),
                            assay = "RNA",
                            slot = "data",
                            ncores = 1,
                            add_to_object = TRUE) {

  method <- match.arg(method)
 
  # Get expression matrix based on object type
  if (inherits(object, "Seurat")) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Package 'Seurat' is required")
    }
    expr <- Seurat::GetAssayData(object, assay = assay, slot = slot)
    expr <- as.matrix(expr)
   
  } else if (inherits(object, "SingleCellExperiment")) {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("Package 'SingleCellExperiment' is required")
    }
    expr <- SingleCellExperiment::logcounts(object)
    expr <- as.matrix(expr)
   
  } else if (is.matrix(object) || inherits(object, "dgCMatrix")) {
    expr <- as.matrix(object)
    add_to_object <- FALSE
   
  } else {
    stop("object must be a Seurat object, SingleCellExperiment, or matrix")
  }
 
  # Get gene sets
  genesets <- get_death_geneset(pathways, type = "all")


  library(SeuratExtend)
  # 根据物种处理
  if(species == "mouse"){
    genesets_use <- lapply(genesets, function(gs){
      unique(na.omit(
        SeuratExtend::HumanToMouseGenesymbol(gs)$MGI.symbol
      ))
    })
  } else {
    genesets_use <- genesets
  }


  
  if (!is.list(genesets)) {
    genesets <- list(genesets)
    names(genesets) <- pathways
  }
 
  # Calculate scores based on method
  if (method == "aucell") {
    scores <- .sc_score_aucell(expr, genesets, ncores)
   
  } else if (method == "seurat") {
    scores <- .sc_score_seurat(expr, genesets)
   
  } else if (method == "mean") {
    scores <- .sc_score_mean(expr, genesets)
   
  } else if (method == "zscore") {
    scores <- .sc_score_zscore(expr, genesets)
  }
 
  # Add prefix to column names
  colnames(scores) <- paste0("death_", colnames(scores))
 
  # Add to object or return
if (add_to_object) {
    if (inherits(object, "Seurat")) {
      for (col in colnames(scores)) {
        object[[col]] <- scores[, col]
      }
      return(object)
     
    } else if (inherits(object, "SingleCellExperiment")) {
      for (col in colnames(scores)) {
        SingleCellExperiment::colData(object)[[col]] <- scores[, col]
      }
      return(object)
    }
  }
 
  return(scores)
}


#' Internal: AUCell scoring for single-cell
#' @keywords internal
.sc_score_aucell <- function(expr, genesets, ncores = 1) {
 
  if (!requireNamespace("AUCell", quietly = TRUE)) {
    stop("Package 'AUCell' is required. Install with: BiocManager::install('AUCell')")
  }
 
  # Build rankings
  rankings <- AUCell::AUCell_buildRankings(expr, nCores = ncores, plotStats = FALSE)
 
  # Calculate AUC
  auc <- AUCell::AUCell_calcAUC(genesets, rankings, nCores = ncores)
 
  # Extract scores
  scores <- t(AUCell::getAUC(auc))
  scores <- as.data.frame(scores)
 
  return(scores)
}


#' Internal: Seurat AddModuleScore-like scoring
#' @keywords internal
.sc_score_seurat <- function(expr, genesets) {
 
  scores <- matrix(NA, nrow = ncol(expr), ncol = length(genesets))
  rownames(scores) <- colnames(expr)
  colnames(scores) <- names(genesets)
 
  # Calculate background
  gene_means <- rowMeans(expr)
  gene_bins <- cut(gene_means, breaks = 25, labels = FALSE)
 
  for (i in seq_along(genesets)) {
    genes <- intersect(genesets[[i]], rownames(expr))
   
    if (length(genes) < 3) {
      scores[, i] <- NA
      next
    }
   
    # Get control genes (same expression level distribution)
    pathway_bins <- gene_bins[genes]
    control_genes <- lapply(unique(pathway_bins), function(b) {
      pool <- names(gene_bins)[gene_bins == b]
      pool <- setdiff(pool, genes)
      if (length(pool) > 100) pool <- sample(pool, 100)
      pool
    })
    control_genes <- unique(unlist(control_genes))
   
    # Calculate scores
    pathway_mean <- colMeans(expr[genes, , drop = FALSE])
    control_mean <- colMeans(expr[control_genes, , drop = FALSE])
   
    scores[, i] <- pathway_mean - control_mean
  }
 
  return(as.data.frame(scores))
}


#' Internal: Mean expression scoring for single-cell
#' @keywords internal
.sc_score_mean <- function(expr, genesets) {
 
  scores <- matrix(NA, nrow = ncol(expr), ncol = length(genesets))
  rownames(scores) <- colnames(expr)
  colnames(scores) <- names(genesets)
 
  for (i in seq_along(genesets)) {
    genes <- intersect(genesets[[i]], rownames(expr))
    if (length(genes) >= 3) {
      scores[, i] <- colMeans(expr[genes, , drop = FALSE])
    }
  }
 
  return(as.data.frame(scores))
}


#' Internal: Z-score for single-cell
#' @keywords internal
.sc_score_zscore <- function(expr, genesets) {
 
  # Z-score normalize expression
  expr_scaled <- t(scale(t(expr)))
 
  scores <- matrix(NA, nrow = ncol(expr), ncol = length(genesets))
  rownames(scores) <- colnames(expr)
  colnames(scores) <- names(genesets)
 
  for (i in seq_along(genesets)) {
    genes <- intersect(genesets[[i]], rownames(expr_scaled))
    if (length(genes) >= 3) {
      scores[, i] <- colMeans(expr_scaled[genes, , drop = FALSE], na.rm = TRUE)
    }
  }
 
  return(as.data.frame(scores))
}


#' Visualize Death Scores on UMAP/tSNE
#'
#' Create feature plots for death pathway scores on dimensionality reduction.
#'
#' @param object A Seurat object with death scores calculated.
#' @param pathways Character vector of pathways to plot.
#' @param reduction Which dimensionality reduction to use. Default is "umap".
#' @param ncol Number of columns in the plot grid.
#' @param colors Color gradient for scores.
#' @param pt_size Point size. Default is 0.5.
#'
#' @return A ggplot object or list of plots.
#'
#' @import ggplot2
#' @export
#'
sc_plot_death_feature <- function(object,
                                   pathways = c("ferroptosis", "pyroptosis"),
                                   reduction = "umap",
                                   ncol = 2,
                                   colors = c("lightgrey", "blue", "red"),
                                   pt_size = 0.5) {
 
  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }
 
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required")
  }
 
  # Get coordinates
  coords <- as.data.frame(Seurat::Embeddings(object, reduction = reduction))
  colnames(coords) <- c("Dim1", "Dim2")
 
  # Create plots
  plots <- list()
 
  for (pathway in pathways) {
    score_col <- paste0("death_", pathway)
   
    if (!score_col %in% colnames(object@meta.data)) {
      warning("Score '", score_col, "' not found. Run sc_death_score() first.")
      next
    }
   
    plot_data <- coords
    plot_data$score <- object@meta.data[[score_col]]
   
    p <- ggplot(plot_data, aes(x = Dim1, y = Dim2, color = score)) +
      geom_point(size = pt_size) +
      scale_color_gradientn(colors = colors, name = "Score") +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      ) +
      labs(
        x = paste0(toupper(reduction), "_1"),
        y = paste0(toupper(reduction), "_2"),
        title = tools::toTitleCase(pathway)
      )
   
    plots[[pathway]] <- p
  }
 
  # Combine plots
  if (length(plots) == 1) {
    return(plots[[1]])
  }
 
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- patchwork::wrap_plots(plots, ncol = ncol)
    return(combined)
  } else {
    return(plots)
  }
}


#' Compare Death Scores Across Cell Types/Clusters
#'
#' Compare cell death pathway scores across different cell populations.
#'
#' @param object A Seurat object with death scores.
#' @param group_by Column in metadata to group cells by (e.g., "seurat_clusters").
#' @param pathways Character vector of pathways to compare.
#' @param test_method Statistical test: "wilcox" or "t.test".
#'
#' @return A data frame with comparison results.
#'
#' @export
#'
sc_compare_clusters <- function(object,
                                 group_by = "seurat_clusters",
                                 pathways = "all",
                                 test_method = c("wilcox", "t.test")) {
 
  test_method <- match.arg(test_method)
 
  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }
 
  if (!group_by %in% colnames(object@meta.data)) {
    stop("Column '", group_by, "' not found in metadata")
  }
 
  # Get pathways
  if (pathways[1] == "all") {
    score_cols <- grep("^death_", colnames(object@meta.data), value = TRUE)
    pathways <- gsub("^death_", "", score_cols)
  } else {
    score_cols <- paste0("death_", pathways)
  }
 
  groups <- unique(object@meta.data[[group_by]])
  results <- list()
 
  for (pathway in pathways) {
    score_col <- paste0("death_", pathway)
   
    if (!score_col %in% colnames(object@meta.data)) {
      next
    }
   
    scores <- object@meta.data[[score_col]]
    group_labels <- object@meta.data[[group_by]]
   
    # Calculate statistics for each group
    for (g in groups) {
      group_scores <- scores[group_labels == g]
      other_scores <- scores[group_labels != g]
     
      # Statistical test
      if (test_method == "wilcox") {
        test <- wilcox.test(group_scores, other_scores)
      } else {
        test <- t.test(group_scores, other_scores)
      }
     
      results[[paste(pathway, g, sep = "_")]] <- data.frame(
        pathway = pathway,
        cluster = g,
        n_cells = length(group_scores),
        mean_score = mean(group_scores, na.rm = TRUE),
        median_score = median(group_scores, na.rm = TRUE),
        mean_other = mean(other_scores, na.rm = TRUE),
        log2FC = log2(mean(group_scores, na.rm = TRUE) / mean(other_scores, na.rm = TRUE)),
        p_value = test$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
 
  result_df <- do.call(rbind, results)
  rownames(result_df) <- NULL
 
  # Adjust p-values
  result_df$p_adjust <- p.adjust(result_df$p_value, method = "BH")
 
  # Sort
  result_df <- result_df[order(result_df$pathway, result_df$p_value), ]
 
  return(result_df)
}


#' Plot Death Scores by Cluster - Heatmap
#'
#' Create a heatmap showing average death scores for each cluster.
#'
#' @param object A Seurat object with death scores.
#' @param group_by Column in metadata to group cells by.
#' @param pathways Pathways to include.
#' @param scale Whether to scale scores. Default is TRUE.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @export
#'
sc_plot_cluster_heatmap <- function(object,
                                     group_by = "seurat_clusters",
                                     pathways = "all",
                                     scale = TRUE) {
 
  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }
 
  # Get pathways
  if (pathways[1] == "all") {
    score_cols <- grep("^death_", colnames(object@meta.data), value = TRUE)
  } else {
    score_cols <- paste0("death_", pathways)
    score_cols <- score_cols[score_cols %in% colnames(object@meta.data)]
  }
 
  if (length(score_cols) == 0) {
    stop("No death scores found. Run sc_death_score() first.")
  }
 
  # Calculate mean scores per cluster
  clusters <- unique(object@meta.data[[group_by]])
  mean_scores <- matrix(NA, nrow = length(clusters), ncol = length(score_cols))
  rownames(mean_scores) <- clusters
  colnames(mean_scores) <- gsub("^death_", "", score_cols)
 
  for (i in seq_along(clusters)) {
    idx <- object@meta.data[[group_by]] == clusters[i]
    mean_scores[i, ] <- colMeans(object@meta.data[idx, score_cols, drop = FALSE],
                                  na.rm = TRUE)
  }
 
  # Scale if requested
  if (scale) {
    mean_scores <- scale(mean_scores)
  }
 
  # Convert to long format
  plot_data <- as.data.frame(mean_scores)
  plot_data$cluster <- rownames(plot_data)
  plot_data <- tidyr::pivot_longer(
    plot_data,
    cols = -cluster,
    names_to = "pathway",
    values_to = "score"
  )
 
  # Create heatmap
  p <- ggplot(plot_data, aes(x = pathway, y = cluster, fill = score)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, name = ifelse(scale, "Z-score", "Score")
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    ) +
    labs(
      x = "Death Pathway",
      y = group_by,
      title = "Cell Death Scores by Cluster"
    )
 
  return(p)
}


#' Plot Death Scores by Cluster - Violin
#'
#' Create violin plots showing death score distributions for each cluster.
#'
#' @param object A Seurat object with death scores.
#' @param pathways Pathways to plot.
#' @param group_by Column in metadata to group cells by.
#' @param ncol Number of columns in plot grid.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @export
#'
sc_plot_cluster_violin <- function(object,
                                    pathways = c("ferroptosis", "pyroptosis"),
                                    group_by = "seurat_clusters",
                                    ncol = 2) {
 
  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }
 
  plots <- list()
 
  for (pathway in pathways) {
    score_col <- paste0("death_", pathway)
   
    if (!score_col %in% colnames(object@meta.data)) {
      warning("Score '", score_col, "' not found")
      next
    }
   
    plot_data <- data.frame(
      cluster = object@meta.data[[group_by]],
      score = object@meta.data[[score_col]]
    )
   
    p <- ggplot(plot_data, aes(x = cluster, y = score, fill = cluster)) +
      geom_violin(scale = "width") +
      geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
      theme_bw() +
      theme(legend.position = "none") +
      labs(
        x = group_by,
        y = "Score",
        title = tools::toTitleCase(pathway)
      )
   
    plots[[pathway]] <- p
  }
 
  if (length(plots) == 1) {
    return(plots[[1]])
  }
 
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- patchwork::wrap_plots(plots, ncol = ncol)
    return(combined)
  }
 
  return(plots)
}


#' Identify Cells with High Death Pathway Activity
#'
#' Identify cells with significantly high death pathway scores.
#'
#' @param object A Seurat object with death scores.
#' @param pathway Which pathway to use.
#' @param method Method for identifying high cells: "quantile", "mad", or "kmeans".
#' @param threshold For quantile method, the quantile threshold. Default is 0.9.
#' @param n_mad For MAD method, number of MADs above median. Default is 2.
#'
#' @return A character vector of cell names with high scores.
#'
#' @export
#'
sc_identify_high_death <- function(object,
                                    pathway = "ferroptosis",
                                    method = c("quantile", "mad", "kmeans"),
                                    threshold = 0.9,
                                    n_mad = 2) {
 
  method <- match.arg(method)
 
  score_col <- paste0("death_", pathway)
 
  if (!score_col %in% colnames(object@meta.data)) {
    stop("Score '", score_col, "' not found. Run sc_death_score() first.")
  }
 
  scores <- object@meta.data[[score_col]]
  names(scores) <- colnames(object)
 
  if (method == "quantile") {
    cutoff <- quantile(scores, probs = threshold, na.rm = TRUE)
    high_cells <- names(scores)[scores >= cutoff]
   
  } else if (method == "mad") {
    med <- median(scores, na.rm = TRUE)
    mad_val <- mad(scores, na.rm = TRUE)
    cutoff <- med + n_mad * mad_val
    high_cells <- names(scores)[scores >= cutoff]
   
  } else if (method == "kmeans") {
    km <- kmeans(scores[!is.na(scores)], centers = 2)
    high_cluster <- which.max(km$centers)
    high_cells <- names(scores)[!is.na(scores)][km$cluster == high_cluster]
  }
 
  return(high_cells)
}


#' Differential Death Analysis Between Conditions
#'
#' Compare death pathway scores between conditions within each cluster.
#'
#' @param object A Seurat object with death scores.
#' @param condition_col Column containing condition information.
#' @param cluster_col Column containing cluster information.
#' @param pathways Pathways to analyze.
#'
#' @return A data frame with differential analysis results.
#'
#' @export
#'
sc_diff_death <- function(object,
                           condition_col,
                           cluster_col = "seurat_clusters",
                           pathways = "all") {
 
  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }
 
  if (!condition_col %in% colnames(object@meta.data)) {
    stop("Column '", condition_col, "' not found")
  }
 
  conditions <- unique(object@meta.data[[condition_col]])
  if (length(conditions) != 2) {
    stop("Exactly 2 conditions required")
  }
 
  clusters <- unique(object@meta.data[[cluster_col]])
 
  # Get pathways
  if (pathways[1] == "all") {
    score_cols <- grep("^death_", colnames(object@meta.data), value = TRUE)
    pathways <- gsub("^death_", "", score_cols)
  }
 
  results <- list()
 
  for (cluster in clusters) {
    cluster_idx <- object@meta.data[[cluster_col]] == cluster
   
    for (pathway in pathways) {
      score_col <- paste0("death_", pathway)
     
      if (!score_col %in% colnames(object@meta.data)) next
     
      scores <- object@meta.data[[score_col]]
      conds <- object@meta.data[[condition_col]]
     
      scores1 <- scores[cluster_idx & conds == conditions[1]]
      scores2 <- scores[cluster_idx & conds == conditions[2]]
     
      if (length(scores1) < 3 || length(scores2) < 3) next
     
      test <- wilcox.test(scores2, scores1)
     
      results[[paste(cluster, pathway, sep = "_")]] <- data.frame(
        cluster = cluster,
        pathway = pathway,
        n_cond1 = length(scores1),
        n_cond2 = length(scores2),
        mean_cond1 = mean(scores1, na.rm = TRUE),
        mean_cond2 = mean(scores2, na.rm = TRUE),
        log2FC = log2(mean(scores2, na.rm = TRUE) / mean(scores1, na.rm = TRUE)),
        p_value = test$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
 
  result_df <- do.call(rbind, results)
  rownames(result_df) <- NULL
 
  result_df$p_adjust <- p.adjust(result_df$p_value, method = "BH")
  result_df <- result_df[order(result_df$p_value), ]
 
  # Rename columns
  names(result_df)[5:6] <- paste0("mean_", conditions)
 
  return(result_df)
}


#' Plot Trajectory Death Scores
#'
#' Visualize death scores along a pseudotime trajectory.
#'
#' @param object A Seurat object with pseudotime and death scores.
#' @param pseudotime_col Column containing pseudotime values.
#' @param pathways Pathways to plot.
#' @param smooth Whether to add smoothed trend line. Default is TRUE.
#' @param ncol Number of columns in plot grid.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @export
#'
sc_plot_trajectory <- function(object,
                                pseudotime_col = "pseudotime",
                                pathways = c("ferroptosis", "apoptosis"),
                                smooth = TRUE,
                                ncol = 2) {
 
  if (!pseudotime_col %in% colnames(object@meta.data)) {
    stop("Pseudotime column '", pseudotime_col, "' not found")
  }
 
  plots <- list()
 
  for (pathway in pathways) {
    score_col <- paste0("death_", pathway)
   
    if (!score_col %in% colnames(object@meta.data)) {
      warning("Score '", score_col, "' not found")
      next
    }
   
    plot_data <- data.frame(
      pseudotime = object@meta.data[[pseudotime_col]],
      score = object@meta.data[[score_col]]
    )
   
    p <- ggplot(plot_data, aes(x = pseudotime, y = score)) +
      geom_point(alpha = 0.3, size = 0.5) +
      theme_bw() +
      labs(
        x = "Pseudotime",
        y = "Score",
        title = tools::toTitleCase(pathway)
      )
   
    if (smooth) {
      p <- p + geom_smooth(method = "loess", color = "red", se = TRUE)
    }
   
    plots[[pathway]] <- p
  }
 
  if (length(plots) == 1) {
    return(plots[[1]])
  }
 
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- patchwork::wrap_plots(plots, ncol = ncol)
    return(combined)
  }
 
  return(plots)
}
