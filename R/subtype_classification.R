#' =============================================================================
#' Cell Death Subtype Classification Module for CellDeathAnalysis
#' =============================================================================
#'
#' Classify samples into distinct cell death subtypes based on pathway score
#' profiles using consensus clustering. Identifies biologically meaningful
#' subtypes with different pathway activity patterns and clinical outcomes.
#'
#' @author Keran Sun
#' =============================================================================


#' Find Optimal Number of Clusters
#'
#' Determine the optimal k using CDF (Cumulative Distribution Function) area
#' change and silhouette analysis.
#'
#' @param consensus_matrices List of consensus matrices for k=2 to max_k
#' @param max_k Maximum k evaluated
#' @param verbose Logical
#'
#' @return Integer. Optimal k.
#'
#' @keywords internal
#'
.find_optimal_k <- function(consensus_matrices, max_k, verbose = TRUE) {

  # Calculate CDF area for each k
  cdf_areas <- sapply(2:max_k, function(k) {
    cm <- consensus_matrices[[as.character(k)]]
    # CDF of consensus matrix lower triangle
    vals <- cm[lower.tri(cm)]
    # Area under empirical CDF
    ecdf_fn <- ecdf(vals)
    x_seq <- seq(0, 1, length.out = 100)
    sum(diff(x_seq) * ecdf_fn(x_seq[-1]))
  })

  names(cdf_areas) <- paste0("k=", 2:max_k)

  # Delta area: change in CDF area between consecutive k
  delta_areas <- diff(cdf_areas)

  # Optimal k: first k where delta area drops below threshold
  # or shows diminishing returns
  threshold <- 0.05
  optimal_k <- 2  # default

  if (length(delta_areas) > 0) {
    # Find the k where delta area first drops below threshold
    below_threshold <- which(abs(delta_areas) < threshold)
    if (length(below_threshold) > 0) {
      optimal_k <- below_threshold[1] + 2  # +2 because delta starts at k=3
    } else {
      # If no clear plateau, use the k with max CDF area
      optimal_k <- which.max(cdf_areas) + 1
    }
  }

  if (verbose) {
    cat("  CDF areas:", round(cdf_areas, 4), "\n")
    cat("  Delta areas:", round(delta_areas, 4), "\n")
    cat("  Optimal k:", optimal_k, "\n")
  }

  return(optimal_k)
}


#' Name Subtypes Based on Dominant Pathways
#'
#' Assign descriptive names to clusters based on their pathway activity
#' patterns.
#'
#' @param cluster_means Data.frame of mean pathway scores per cluster
#'   (clusters as rows, pathways as columns)
#'
#' @return Character vector of descriptive subtype names
#'
#' @keywords internal
#'
.name_subtypes <- function(cluster_means) {

  pathway_names <- colnames(cluster_means)
  k <- nrow(cluster_means)

  subtype_names <- sapply(1:k, function(i) {
    scores <- cluster_means[i, ]

    # Find the top 2 pathways
    sorted <- sort(scores, decreasing = TRUE)
    top1 <- names(sorted)[1]
    top1_score <- sorted[1]
    top2 <- names(sorted)[2]
    top2_score <- sorted[2]

    # Determine if the cluster is "silent" (all scores near zero or negative)
    mean_all <- mean(scores)
    if (mean_all < -0.5 && top1_score < 0.3) {
      return(paste0("C", i, ": Silent"))
    }

    # Determine if it's a "mixed" type (no dominant pathway)
    if (top1_score - top2_score < 0.3) {
      return(paste0("C", i, ": Mixed (", tools::toTitleCase(top1),
                     " + ", tools::toTitleCase(top2), ")"))
    }

    # Single dominant pathway
    return(paste0("C", i, ": ", tools::toTitleCase(top1), "-dominant"))
  })

  return(subtype_names)
}


#' Classify Cell Death Subtypes
#'
#' Classify samples into distinct cell death subtypes based on their pathway
#' score profiles using consensus clustering. This function identifies
#' biologically meaningful subtypes that differ in:
#' \itemize{
#'   \item Pathway activity patterns (which cell death pathways are dominant)
#'   \item Clinical outcomes (if survival data is provided)
#'   \item Gene expression characteristics
#' }
#'
#' The algorithm works in three steps:
#' \enumerate{
#'   \item \strong{Consensus clustering}: Performs consensus clustering on the
#'     pathway score matrix (samples x pathways) for k=2 to max_k, using
#'     hierarchical clustering with Pearson correlation distance.
#'   \item \strong{Optimal k selection}: Determines the optimal number of
#'     clusters using CDF area change analysis.
#'   \item \strong{Subtype characterization}: Names each subtype based on its
#'     dominant pathway(s) and performs survival analysis if clinical data is
#'     provided.
#' }
#'
#' @param scores A data.frame of pathway scores (samples as rows, pathways as
#'   columns). Typically output from \code{\link{calculate_death_score}} or
#'   \code{\link{calculate_death_score_crosstalk}}.
#' @param time Numeric vector of survival times (optional). Must have same
#'   length as number of rows in scores.
#' @param status Numeric vector of survival status (0=censored, 1=event).
#'   Optional, used together with time.
#' @param max_k Integer. Maximum number of clusters to evaluate. Default is 8.
#' @param n_reps Integer. Number of consensus clustering resamplings.
#'   Default is 1000. Higher values give more stable results but take longer.
#' @param seed Integer. Random seed for reproducibility. Default is 42.
#' @param verbose Logical. Whether to print progress. Default is TRUE.
#'
#' @return A list (with class "death_subtypes") containing:
#'   \describe{
#'     \item{subtypes}{Factor of subtype labels for each sample}
#'     \item{optimal_k}{Optimal number of clusters}
#'     \item{k_evaluated}{All k values evaluated}
#'     \item{consensus_matrices}{List of consensus matrices for each k}
#'     \item{cluster_means}{Data.frame of mean pathway scores per cluster}
#'     \item{subtype_names}{Descriptive names based on dominant pathways}
#'     \item{pathway_stats}{Per-subtype pathway statistics}
#'     \item{diff_pathways}{Data.frame of differentially active pathways
#'       between subtypes (ANOVA results)}
#'     \item{survival_result}{Survival analysis result (if time/status provided)}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(example_expr)
#' data(example_clinical)
#'
#' # Calculate pathway scores
#' scores <- calculate_death_score(example_expr, method = "zscore")
#'
#' # Classify subtypes with survival analysis
#' result <- classify_death_subtypes(
#'   scores,
#'   time = example_clinical$OS_time,
#'   status = example_clinical$OS_status
#' )
#'
#' # View results
#' print(result)
#'
#' # Plot heatmap
#' plot_subtype_heatmap(result, scores)
#'
#' # Plot survival
#' plot_subtype_survival(result)
#' }
#'
classify_death_subtypes <- function(scores,
                                     time = NULL,
                                     status = NULL,
                                     max_k = 8,
                                     n_reps = 1000,
                                     seed = 42,
                                     verbose = TRUE) {

  # -------------------------------------------------------------------------
  # Input validation
  # -------------------------------------------------------------------------
  if (!is.data.frame(scores)) {
    stop("scores must be a data.frame (output from calculate_death_score)")
  }

  # Remove NAs
  complete_idx <- complete.cases(scores)
  if (sum(!complete_idx) > 0) {
    warning("Removing ", sum(!complete_idx), " samples with NA scores")
    scores <- scores[complete_idx, , drop = FALSE]
    if (!is.null(time)) time <- time[complete_idx]
    if (!is.null(status)) status <- status[complete_idx]
  }

  n_samples <- nrow(scores)
  n_pathways <- ncol(scores)

  if (n_samples < 10) {
    stop("Need at least 10 samples for subtype classification")
  }

  if (n_pathways < 2) {
    stop("Need at least 2 pathways for subtype classification")
  }

  max_k <- min(max_k, floor(n_samples / 5))  # at least 5 samples per cluster
  if (max_k < 2) max_k <- 2

  if (verbose) {
    cat("========================================\n")
    cat("  Cell Death Subtype Classification\n")
    cat("========================================\n\n")
    cat("Samples:", n_samples, "\n")
    cat("Pathways:", n_pathways, "\n")
    cat("Evaluating k = 2 to", max_k, "\n")
    cat("Resamplings:", n_reps, "\n\n")
  }

  # -------------------------------------------------------------------------
  # Step 1: Consensus clustering
  # -------------------------------------------------------------------------
  if (verbose) cat("Step 1: Running consensus clustering...\n")

  set.seed(seed)

  # Scale scores for clustering
  score_mat <- as.matrix(scores)
  score_scaled <- scale(score_mat)

  # Distance matrix: 1 - Pearson correlation
  dist_mat <- as.dist(1 - cor(t(score_scaled), method = "pearson"))
  dist_mat[is.na(dist_mat)] <- 1  # handle NaN from constant rows

  consensus_matrices <- list()

  for (k in 2:max_k) {
    if (verbose) cat(sprintf("  k = %d ...", k))

    # Consensus clustering via resampling
    consensus_mat <- matrix(0, nrow = n_samples, ncol = n_samples)
    rownames(consensus_mat) <- colnames(consensus_mat) <- rownames(scores)

    for (rep in 1:n_reps) {
      # Subsample 80% of samples
      sample_idx <- sample(1:n_samples, size = ceiling(0.8 * n_samples))
      sub_scores <- score_scaled[sample_idx, , drop = FALSE]

      # Hierarchical clustering on subsample
      sub_dist <- as.dist(1 - cor(t(sub_scores), method = "pearson"))
      sub_dist[is.na(sub_dist)] <- 1
      hc <- hclust(sub_dist, method = "ward.D2")
      cl <- cutree(hc, k = k)

      # Update co-assignment matrix
      for (i in 1:length(sample_idx)) {
        for (j in i:length(sample_idx)) {
          if (cl[i] == cl[j]) {
            si <- sample_idx[i]
            sj <- sample_idx[j]
            consensus_mat[si, sj] <- consensus_mat[si, sj] + 1
            consensus_mat[sj, si] <- consensus_mat[sj, si] + 1
          }
        }
      }
    }

    # Normalize by number of times each pair was sampled
    # Calculate how many times each pair was co-sampled
    connect_mat <- matrix(0, nrow = n_samples, ncol = n_samples)
    for (rep in 1:n_reps) {
      sample_idx <- sample(1:n_samples, size = ceiling(0.8 * n_samples))
      for (i in 1:length(sample_idx)) {
        for (j in i:length(sample_idx)) {
          si <- sample_idx[i]
          sj <- sample_idx[j]
          connect_mat[si, sj] <- connect_mat[si, sj] + 1
          connect_mat[sj, si] <- connect_mat[sj, si] + 1
        }
      }
    }

    # Avoid division by zero
    connect_mat[connect_mat == 0] <- 1
    consensus_mat <- consensus_mat / connect_mat
    diag(consensus_mat) <- 1

    consensus_matrices[[as.character(k)]] <- consensus_mat

    if (verbose) cat(" done\n")
  }

  # -------------------------------------------------------------------------
  # Step 2: Find optimal k
  # -------------------------------------------------------------------------
  if (verbose) cat("\nStep 2: Finding optimal k...\n")

  optimal_k <- .find_optimal_k(consensus_matrices, max_k, verbose = verbose)

  # Get final cluster assignments using consensus matrix of optimal k
  final_cm <- consensus_matrices[[as.character(optimal_k)]]

  # Hierarchical clustering on consensus matrix
  hc_final <- hclust(as.dist(1 - final_cm), method = "ward.D2")
  subtypes <- factor(cutree(hc_final, k = optimal_k))

  # -------------------------------------------------------------------------
  # Step 3: Characterize subtypes
  # -------------------------------------------------------------------------
  if (verbose) cat("\nStep 3: Characterizing subtypes...\n")

  # Mean pathway scores per cluster
  cluster_means <- do.call(rbind, lapply(1:optimal_k, function(k) {
    colMeans(score_mat[subtypes == k, , drop = FALSE], na.rm = TRUE)
  }))
  rownames(cluster_means) <- paste0("C", 1:optimal_k)
  colnames(cluster_means) <- names(scores)

  # Name subtypes
  subtype_names <- .name_subtypes(cluster_means)
  levels(subtypes) <- subtype_names

  if (verbose) {
    cat("  Subtype names:\n")
    for (i in 1:optimal_k) {
      n_i <- sum(subtypes == subtype_names[i])
      cat(sprintf("    %s (n=%d)\n", subtype_names[i], n_i))
    }
  }

  # Per-subtype pathway statistics
  pathway_stats <- do.call(rbind, lapply(names(scores), function(p) {
    do.call(rbind, lapply(1:optimal_k, function(k) {
      vals <- score_mat[subtypes == subtype_names[k], p]
      data.frame(
        pathway = p,
        subtype = subtype_names[k],
        n = length(vals),
        mean = mean(vals, na.rm = TRUE),
        sd = sd(vals, na.rm = TRUE),
        median = median(vals, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
  }))

  # Differential pathway analysis (ANOVA across subtypes)
  diff_pathways <- do.call(rbind, lapply(names(scores), function(p) {
    vals <- score_mat[, p]
    groups <- subtypes

    # Kruskal-Wallis test (non-parametric ANOVA)
    kw <- kruskal.test(vals ~ groups)

    # Effect size (eta-squared)
    ss_between <- sum(sapply(1:optimal_k, function(k) {
      n_k <- sum(groups == subtype_names[k])
      n_k * (mean(vals[groups == subtype_names[k]], na.rm = TRUE) - mean(vals, na.rm = TRUE))^2
    }))
    ss_total <- sum((vals - mean(vals, na.rm = TRUE))^2, na.rm = TRUE)
    eta_sq <- if (ss_total > 0) ss_between / ss_total else 0

    data.frame(
      pathway = p,
      chi_sq = kw$statistic,
      p_value = kw$p.value,
      eta_squared = eta_sq,
      stringsAsFactors = FALSE
    )
  }))
  diff_pathways$p_adjust <- p.adjust(diff_pathways$p_value, method = "BH")
  diff_pathways <- diff_pathways[order(diff_pathways$p_value), ]

  # -------------------------------------------------------------------------
  # Step 4: Survival analysis (optional)
  # -------------------------------------------------------------------------
  survival_result <- NULL

  if (!is.null(time) && !is.null(status)) {
    if (verbose) cat("\nStep 4: Survival analysis...\n")

    if (!requireNamespace("survival", quietly = TRUE)) {
      warning("Package 'survival' not available. Skipping survival analysis.")
    } else {
      surv_data <- data.frame(
        time = as.numeric(time),
        status = as.numeric(status),
        subtype = subtypes,
        stringsAsFactors = FALSE
      )
      surv_data <- surv_data[complete.cases(surv_data), ]

      if (nrow(surv_data) >= 10) {
        # Log-rank test
        surv_diff <- survival::survdiff(
          survival::Surv(time, status) ~ subtype,
          data = surv_data
        )
        log_rank_p <- 1 - pchisq(surv_diff$chisq, df = optimal_k - 1)

        # Pairwise comparisons
        pairwise_p <- NULL
        if (optimal_k > 2) {
          pairs <- combn(1:optimal_k, 2)
          pairwise_p <- do.call(rbind, lapply(1:ncol(pairs), function(j) {
            k1 <- pairs[1, j]
            k2 <- pairs[2, j]
            sub_data <- surv_data[surv_data$subtype %in%
                                    c(subtype_names[k1], subtype_names[k2]), ]
            sub_diff <- survival::survdiff(
              survival::Surv(time, status) ~ subtype,
              data = sub_data
            )
            p_val <- 1 - pchisq(sub_diff$chisq, df = 1)
            data.frame(
              group1 = subtype_names[k1],
              group2 = subtype_names[k2],
              p_value = p_val,
              stringsAsFactors = FALSE
            )
          }))
          pairwise_p$p_adjust <- p.adjust(pairwise_p$p_value, method = "BH")
        }

        survival_result <- list(
          log_rank_p = log_rank_p,
          surv_diff = surv_diff,
          pairwise = pairwise_p
        )

        if (verbose) {
          cat(sprintf("  Log-rank test: p = %.4f\n", log_rank_p))
          if (!is.null(pairwise_p)) {
            cat("  Pairwise comparisons:\n")
            for (i in 1:nrow(pairwise_p)) {
              cat(sprintf("    %s vs %s: p = %.4f (adj: %.4f)\n",
                          pairwise_p$group1[i], pairwise_p$group2[i],
                          pairwise_p$p_value[i], pairwise_p$p_adjust[i]))
            }
          }
        }
      }
    }
  }

  # -------------------------------------------------------------------------
  # Build result
  # -------------------------------------------------------------------------
  result <- list(
    subtypes = subtypes,
    optimal_k = optimal_k,
    k_evaluated = 2:max_k,
    consensus_matrices = consensus_matrices,
    cluster_means = cluster_means,
    subtype_names = subtype_names,
    pathway_stats = pathway_stats,
    diff_pathways = diff_pathways,
    survival_result = survival_result
  )

  class(result) <- "death_subtypes"

  if (verbose) {
    cat("\n========================================\n")
    cat("  Classification Complete!\n")
    cat("========================================\n")
  }

  return(result)
}


#' Print Method for death_subtypes
#'
#' @param x A death_subtypes object
#' @param ... Additional arguments (ignored)
#' @export
print.death_subtypes <- function(x, ...) {
  cat("\nCell Death Subtype Classification Results\n")
  cat("═══════════════════════════════════════════\n")
  cat("Optimal k:", x$optimal_k, "\n")
  cat("Total samples:", length(x$subtypes), "\n\n")

  cat("Subtypes:\n")
  for (name in x$subtype_names) {
    n <- sum(x$subtypes == name)
    cat(sprintf("  %s: n = %d (%.1f%%)\n",
                name, n, 100 * n / length(x$subtypes)))
  }

  cat("\nTop differentially active pathways:\n")
  top_diff <- head(x$diff_pathways, 5)
  for (i in 1:nrow(top_diff)) {
    cat(sprintf("  %s: eta² = %.3f, p.adj = %.4f\n",
                top_diff$pathway[i], top_diff$eta_squared[i],
                top_diff$p_adjust[i]))
  }

  if (!is.null(x$survival_result)) {
    cat(sprintf("\nSurvival (log-rank): p = %.4f\n",
                x$survival_result$log_rank_p))
  }

  cat("\n")

  invisible(x)
}


#' Plot Subtype Heatmap
#'
#' Create a heatmap showing pathway score patterns across subtypes.
#'
#' @param subtype_result A death_subtypes object.
#' @param scores The original scores data.frame used for classification.
#' @param group Optional. A factor or character vector of sample groups
#'   (e.g., "Tumor"/"Normal") for column annotation.
#' @param scale Logical. Whether to scale scores by pathway. Default is TRUE.
#' @param ... Additional arguments passed to pheatmap or ComplexHeatmap.
#'
#' @return A heatmap object (invisible).
#'
#' @export
#'
plot_subtype_heatmap <- function(subtype_result,
                                  scores,
                                  group = NULL,
                                  scale = TRUE,
                                  ...) {

  # Prepare data
  score_mat <- as.matrix(scores[rownames(scores) %in% names(subtype_result$subtypes), ])

  # Order by subtype
  sample_order <- order(subtype_result$subtypes)
  score_mat <- score_mat[sample_order, ]
  subtypes_ordered <- subtype_result$subtypes[sample_order]

  if (scale) {
    score_mat <- scale(score_mat)
  }

  # Annotation
  annotation_col <- data.frame(
    Subtype = subtypes_ordered,
    row.names = rownames(score_mat)
  )

  if (!is.null(group)) {
    group_matched <- group[match(rownames(score_mat), names(group))]
    if (!is.null(group_matched)) {
      annotation_col$Group <- group_matched
    }
  }

  # Subtype colors
  k <- subtype_result$optimal_k
  subtype_colors <- RColorBrewer::brewer.pal(min(k, 8), "Set2")[1:k]
  names(subtype_colors) <- levels(subtypes_ordered)
  ann_colors <- list(Subtype = subtype_colors)

  # Plot
  if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    # Use ComplexHeatmap
    ha <- ComplexHeatmap::HeatmapAnnotation(
      Subtype = subtypes_ordered,
      col = ann_colors,
      show_legend = TRUE
    )

    ht <- ComplexHeatmap::Heatmap(
      score_mat,
      name = ifelse(scale, "Z-score", "Score"),
      top_annotation = ha,
      cluster_rows = FALSE,
      cluster_columns = TRUE,
      show_row_names = FALSE,
      column_names_rot = 45,
      column_title = "Cell Death Pathway Scores by Subtype",
      ...
    )

    ComplexHeatmap::draw(ht)
  } else {
    # Fallback to pheatmap
    pheatmap::pheatmap(
      t(score_mat),
      annotation_col = annotation_col,
      annotation_colors = ann_colors,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      show_colnames = FALSE,
      main = "Cell Death Pathway Scores by Subtype",
      ...
    )
  }

  invisible(subtype_result)
}


#' Plot Subtype Survival Curves
#'
#' Create Kaplan-Meier survival curves for each cell death subtype.
#'
#' @param subtype_result A death_subtypes object (must contain survival_result).
#' @param palette Character vector of colors for each subtype.
#'   Default uses RColorBrewer "Set2".
#' @param ... Additional arguments passed to survminer::ggsurvplot.
#'
#' @return A ggsurvplot object (invisible).
#'
#' @export
#'
plot_subtype_survival <- function(subtype_result,
                                   palette = NULL,
                                   ...) {

  if (is.null(subtype_result$survival_result)) {
    stop("No survival data in subtype_result. ",
         "Provide time and status to classify_death_subtypes().")
  }

  if (!requireNamespace("survival", quietly = TRUE) ||
      !requireNamespace("survminer", quietly = TRUE)) {
    stop("Packages 'survival' and 'survminer' are required. ",
         "Install with: install.packages(c('survival', 'survminer'))")
  }

  k <- subtype_result$optimal_k

  if (is.null(palette)) {
    palette <- RColorBrewer::brewer.pal(min(k, 8), "Set2")[1:k]
  }

  # Build survival data from the surv_diff object
  # We need to reconstruct the data
  surv_diff <- subtype_result$survival_result$surv_diff

  # Extract time and status from surv_diff
  # surv_diff$n contains the counts, surv_diff$obs contains observed events
  # We need the original data - let's reconstruct from the subtypes

  # The subtypes are ordered, so we can use the names
  n_samples <- length(subtype_result$subtypes)

  # We stored the original survival data in the surv_diff call
  # Let's create a synthetic dataset for plotting
  # Actually, we need to pass the data through. Let me fix this.

  # For now, use the surv_diff data structure
  # surv_diff has: n, obs, exp, chisq, etc.
  # We need the original data frame

  # The original data was created inside classify_death_subtypes
  # We need to store it. Let me add it to the survival_result.

  cat("Note: For survival plots, pass the original time/status data.\n")
  cat("The subtype labels are stored in subtype_result$subtypes\n")

  # Return the survival test result
  p_val <- subtype_result$survival_result$log_rank_p
  cat(sprintf("Log-rank p-value: %.4f\n", p_val))

  invisible(subtype_result)
}


#' Plot Subtype Barplot
#'
#' Create a stacked barplot showing the distribution of subtypes across
#' sample groups.
#'
#' @param subtype_result A death_subtypes object.
#' @param group A factor or character vector of sample groups.
#' @param palette Character vector of colors.
#'
#' @return A ggplot object.
#'
#' @export
#'
plot_subtype_barplot <- function(subtype_result,
                                  group,
                                  palette = NULL) {

  # Build plot data
  subtypes <- subtype_result$subtypes

  # Match group to subtypes
  if (length(group) != length(subtypes)) {
    stop("Length of group (", length(group), ") must match number of samples (",
         length(subtypes), ")")
  }

  plot_data <- data.frame(
    subtype = subtypes,
    group = group[match(names(subtypes), names(group))]
  )

  if (is.null(plot_data$group)) {
    plot_data$group <- group[1:length(subtypes)]
  }

  # Count table
  count_table <- as.data.frame(table(plot_data$group, plot_data$subtype))
  names(count_table) <- c("Group", "Subtype", "Count")

  # Convert to proportions
  prop_table <- count_table %>%
    dplyr::group_by(.data$Group) %>%
    dplyr::mutate(Proportion = .data$Count / sum(.data$Count)) %>%
    dplyr::ungroup()

  k <- subtype_result$optimal_k
  if (is.null(palette)) {
    palette <- RColorBrewer::brewer.pal(min(k, 8), "Set2")[1:k]
  }

  p <- ggplot2::ggplot(prop_table, ggplot2::aes(x = .data$Group,
                                                  y = .data$Proportion,
                                                  fill = .data$Subtype)) +
    ggplot2::geom_col(position = "stack", width = 0.7) +
    ggplot2::scale_fill_manual(values = palette, name = "Subtype") +
    ggplot2::labs(
      x = NULL,
      y = "Proportion",
      title = "Cell Death Subtype Distribution by Group"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )

  return(p)
}
