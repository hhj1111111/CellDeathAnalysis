#' =============================================================================
#' Enrichment Analysis Module for CellDeathAnalysis
#' =============================================================================

#' Over-Representation Analysis for Death Pathways
#'
#' Perform over-representation analysis (ORA) to test if death pathway genes
#' are enriched in a gene list.
#'
#' @param gene_list Character vector of gene symbols (e.g., DEGs).
#' @param pathways Character vector of pathways to test, or "all".
#' @param background Character vector of background genes. If NULL, uses all
#'   genes in death_genesets.
#' @param min_genes Minimum number of overlapping genes required. Default is 3.
#' @param p_adjust Method for p-value adjustment. Default is "BH".
#'
#' @return A data frame with enrichment results.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example with differentially expressed genes
#' degs <- c("GPX4", "SLC7A11", "ACSL4", "TFRC", "NLRP3", "CASP1", "GSDMD")
#' result <- death_enrich_ora(degs)
#' print(result)
#' }
#'
death_enrich_ora <- function(gene_list,species = 'human',
                              pathways = "all",
                              background = NULL,
                              min_genes = 3,
                              p_adjust = "BH") {
  
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
  
  # Set background
  if (is.null(background)) {
    background <- unique(unlist(genesets))
  }
  
  # Filter gene list to background
  gene_list <- intersect(gene_list, background)
  
  if (length(gene_list) == 0) {
    stop("No genes from gene_list found in background")
  }
  
  # Perform Fisher's exact test for each pathway
  results <- lapply(names(genesets), function(pathway) {
    pathway_genes <- genesets[[pathway]]
    pathway_genes <- intersect(pathway_genes, background)
    
    # Build contingency table
    # a: in list AND in pathway
    # b: in list but NOT in pathway
    # c: NOT in list but in pathway
    # d: NOT in list AND NOT in pathway
    
    a <- length(intersect(gene_list, pathway_genes))
    b <- length(setdiff(gene_list, pathway_genes))
    c <- length(setdiff(pathway_genes, gene_list))
    d <- length(background) - a - b - c
    
    # Skip if too few overlapping genes
    if (a < min_genes) {
      return(NULL)
    }
    
    # Fisher's exact test
    mat <- matrix(c(a, b, c, d), nrow = 2)
    test <- fisher.test(mat, alternative = "greater")
    
    # Calculate fold enrichment
    expected <- length(gene_list) * length(pathway_genes) / length(background)
    fold_enrichment <- a / expected
    
    data.frame(
      pathway = pathway,
      n_overlap = a,
      n_pathway = length(pathway_genes),
      n_query = length(gene_list),
      n_background = length(background),
      fold_enrichment = fold_enrichment,
      p_value = test$p.value,
      genes = paste(intersect(gene_list, pathway_genes), collapse = "/"),
      stringsAsFactors = FALSE
    )
  })
  
  # Combine results
  results <- do.call(rbind, results[!sapply(results, is.null)])
  
  if (is.null(results) || nrow(results) == 0) {
    warning("No pathways with >= ", min_genes, " overlapping genes")
    return(NULL)
  }
  
  # Adjust p-values
  results$p_adjust <- p.adjust(results$p_value, method = p_adjust)
  
  # Sort by p-value
  results <- results[order(results$p_value), ]
  rownames(results) <- NULL
  
  # Add significance
  results$significance <- ifelse(results$p_adjust < 0.001, "***",
                                  ifelse(results$p_adjust < 0.01, "**",
                                         ifelse(results$p_adjust < 0.05, "*", "")))
  
  class(results) <- c("death_ora", "data.frame")
  
  return(results)
}


#' Print Method for death_ora
#' @param x A death_ora object
#' @param ... Additional arguments
#' @export
print.death_ora <- function(x, ...) {
  cat("\n")
  cat("========================================\n")
  cat("  Over-Representation Analysis Results\n")
  cat("========================================\n\n")
  
  cat("Query genes:", x$n_query[1], "\n")
  cat("Background genes:", x$n_background[1], "\n")
  cat("Pathways tested:", nrow(x), "\n")
  cat("Significant (adj.p < 0.05):", sum(x$p_adjust < 0.05), "\n\n")
  
  # Print results
  print_df <- x[, c("pathway", "n_overlap", "fold_enrichment", "p_value", "p_adjust", "significance")]
  print_df$fold_enrichment <- round(print_df$fold_enrichment, 2)
  print_df$p_value <- format.pval(print_df$p_value, digits = 2)
  print_df$p_adjust <- format.pval(print_df$p_adjust, digits = 2)
  
  print(print_df, row.names = FALSE)
  cat("\n")
  
  invisible(x)
}


#' Gene Set Enrichment Analysis for Death Pathways
#'
#' Perform GSEA using ranked gene list against death pathways.
#'
#' @param ranked_genes A named numeric vector of gene scores (e.g., log2FC or
#'   -log10(p) * sign(FC)). Names should be gene symbols.
#' @param pathways Character vector of pathways to test, or "all".
#' @param nperm Number of permutations. Default is 1000.
#' @param min_size Minimum gene set size. Default is 10.
#' @param max_size Maximum gene set size. Default is 500.
#'
#' @return A data frame with GSEA results.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create ranked gene list (e.g., from DESeq2 results)
#' ranked <- setNames(de_results$log2FoldChange, de_results$gene)
#' ranked <- sort(ranked, decreasing = TRUE)
#'
#' result <- death_enrich_gsea(ranked)
#' }
#'
death_enrich_gsea <- function(ranked_genes,
                               pathways = "all",
                               nperm = 1000,
                               min_size = 10,
                               max_size = 500) {
  
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required. Install with: BiocManager::install('fgsea')")
  }
  
  # Get gene sets
  genesets <- get_death_geneset(pathways, type = "all")
  
  if (!is.list(genesets)) {
    genesets <- list(genesets)
    names(genesets) <- pathways
  }
  
  # Ensure ranked_genes is sorted
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  # Run fgsea
  fgsea_res <- fgsea::fgsea(
    pathways = genesets,
    stats = ranked_genes,
    minSize = min_size,
    maxSize = max_size,
    nperm = nperm
  )
  
  # Convert to data frame
  results <- as.data.frame(fgsea_res)
  results <- results[order(results$pval), ]
  
  # Add significance
  results$significance <- ifelse(results$padj < 0.001, "***",
                                  ifelse(results$padj < 0.01, "**",
                                         ifelse(results$padj < 0.05, "*", "")))
  
  # Direction
  results$direction <- ifelse(results$NES > 0, "Up", "Down")
  
  class(results) <- c("death_gsea", "data.frame")
  
  return(results)
}


#' Print Method for death_gsea
#' @param x A death_gsea object
#' @param ... Additional arguments
#' @export
print.death_gsea <- function(x, ...) {
  cat("\n")
  cat("========================================\n")
  cat("  Gene Set Enrichment Analysis Results\n")
  cat("========================================\n\n")
  
  cat("Pathways tested:", nrow(x), "\n")
  cat("Significant (adj.p < 0.05):", sum(x$padj < 0.05, na.rm = TRUE), "\n")
  cat("  Upregulated:", sum(x$padj < 0.05 & x$NES > 0, na.rm = TRUE), "\n")
  cat("  Downregulated:", sum(x$padj < 0.05 & x$NES < 0, na.rm = TRUE), "\n\n")
  
  # Print results
  print_df <- x[, c("pathway", "size", "NES", "pval", "padj", "direction", "significance")]
  print_df$NES <- round(print_df$NES, 3)
  print_df$pval <- format.pval(print_df$pval, digits = 2)
  print_df$padj <- format.pval(print_df$padj, digits = 2)
  
  print(print_df, row.names = FALSE)
  cat("\n")
  
  invisible(x)
}


#' Plot Enrichment Results - Bar Plot
#'
#' Create a bar plot showing enrichment results.
#'
#' @param enrich_result A death_ora or death_gsea object.
#' @param top_n Number of top pathways to show.
#' @param show_all Logical. Whether to show all pathways or only significant.
#' @param colors Colors for bars.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @export
#'
plot_enrichment_bar <- function(enrich_result,
                                 top_n = 10,
                                 show_all = FALSE,
                                 colors = c("Up" = "#E41A1C", "Down" = "#377EB8")) {
  
  if (inherits(enrich_result, "death_ora")) {
    # ORA result
    plot_data <- enrich_result
    if (!show_all) {
      plot_data <- plot_data[plot_data$p_adjust < 0.05, ]
    }
    plot_data <- head(plot_data, top_n)
    
    plot_data$pathway <- factor(plot_data$pathway, 
                                 levels = rev(plot_data$pathway))
    
    p <- ggplot(plot_data, aes(x = -log10(p_adjust), y = pathway)) +
      geom_col(aes(fill = fold_enrichment), width = 0.7) +
      geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
      scale_fill_gradient(low = "lightblue", high = "darkblue", 
                          name = "Fold\nEnrichment") +
      theme_bw() +
      labs(
        x = "-log10(Adjusted p-value)",
        y = NULL,
        title = "Death Pathway Enrichment (ORA)"
      )
    
  } else if (inherits(enrich_result, "death_gsea")) {
    # GSEA result
    plot_data <- enrich_result
    if (!show_all) {
      plot_data <- plot_data[plot_data$padj < 0.05, ]
    }
    
    if (nrow(plot_data) == 0) {
      warning("No significant pathways to plot")
      return(NULL)
    }
    
    plot_data <- head(plot_data[order(abs(plot_data$NES), decreasing = TRUE), ], top_n)
    plot_data$pathway <- factor(plot_data$pathway, 
                                 levels = rev(plot_data$pathway))
    
    p <- ggplot(plot_data, aes(x = NES, y = pathway, fill = direction)) +
      geom_col(width = 0.7) +
      geom_vline(xintercept = 0, color = "black") +
      scale_fill_manual(values = colors) +
      theme_bw() +
      labs(
        x = "Normalized Enrichment Score (NES)",
        y = NULL,
        title = "Death Pathway Enrichment (GSEA)",
        fill = "Direction"
      )
    
  } else {
    stop("enrich_result must be a death_ora or death_gsea object")
  }
  
  return(p)
}


#' Plot Enrichment Results - Dot Plot
#'
#' Create a dot plot showing enrichment results.
#'
#' @param enrich_result A death_ora or death_gsea object.
#' @param top_n Number of top pathways to show.
#' @param size_by Variable for point size: "overlap" or "fold".
#' @param color_by Variable for point color: "pvalue" or "enrichment".
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @export
#'
plot_enrichment_dot <- function(enrich_result,
                                 top_n = 10,
                                 size_by = "overlap",
                                 color_by = "pvalue") {
  
  if (inherits(enrich_result, "death_ora")) {
    plot_data <- head(enrich_result, top_n)
    plot_data$pathway <- factor(plot_data$pathway, 
                                 levels = rev(plot_data$pathway))
    
    p <- ggplot(plot_data, aes(x = fold_enrichment, y = pathway)) +
      geom_point(aes(size = n_overlap, color = -log10(p_adjust))) +
      scale_color_gradient(low = "blue", high = "red", 
                           name = "-log10(p.adj)") +
      scale_size_continuous(name = "Gene Count", range = c(3, 10)) +
      theme_bw() +
      labs(
        x = "Fold Enrichment",
        y = NULL,
        title = "Death Pathway Enrichment"
      )
    
  } else if (inherits(enrich_result, "death_gsea")) {
    plot_data <- head(enrich_result[order(enrich_result$pval), ], top_n)
    plot_data$pathway <- factor(plot_data$pathway, 
                                 levels = rev(plot_data$pathway))
    
    p <- ggplot(plot_data, aes(x = NES, y = pathway)) +
      geom_point(aes(size = size, color = -log10(padj))) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_color_gradient(low = "blue", high = "red", 
                           name = "-log10(p.adj)") +
      scale_size_continuous(name = "Gene Set\nSize", range = c(3, 10)) +
      theme_bw() +
      labs(
        x = "Normalized Enrichment Score",
        y = NULL,
        title = "Death Pathway Enrichment (GSEA)"
      )
  }
  
  return(p)
}


#' Plot GSEA Enrichment Plot for Single Pathway
#'
#' Create classic GSEA enrichment plot for a specific pathway.
#'
#' @param ranked_genes Named numeric vector of ranked genes.
#' @param pathway Character. Name of the pathway.
#' @param title Character. Plot title.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @export
#'
plot_gsea_curve <- function(ranked_genes,
                             pathway,
                             title = NULL) {
  
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required")
  }
  
  # Get pathway genes
  pathway_genes <- get_death_geneset(pathway, type = "all")
  
  # Sort ranked genes
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  # Calculate running sum
  gene_names <- names(ranked_genes)
  in_set <- gene_names %in% pathway_genes
  
  # Running enrichment score
  n <- length(ranked_genes)
  n_hit <- sum(in_set)
  n_miss <- n - n_hit
  
  if (n_hit == 0) {
    stop("No pathway genes found in ranked list")
  }
  
  # Calculate scores
  hit_score <- abs(ranked_genes[in_set])
  hit_score <- hit_score / sum(hit_score)
  miss_score <- 1 / n_miss
  
  running_score <- rep(0, n)
  for (i in 1:n) {
    if (in_set[i]) {
      running_score[i] <- hit_score[sum(in_set[1:i])]
    } else {
      running_score[i] <- -miss_score
    }
  }
  running_score <- cumsum(running_score)
  
  # Prepare plot data
  plot_data <- data.frame(
    rank = 1:n,
    score = running_score,
    in_set = in_set
  )
  
  if (is.null(title)) {
    title <- paste("GSEA:", pathway)
  }
  
  # Create plot
  p <- ggplot(plot_data, aes(x = rank, y = score)) +
    geom_line(color = "green4", linewidth = 1) +
    geom_hline(yintercept = 0, color = "black") +
    geom_segment(data = plot_data[plot_data$in_set, ],
                 aes(x = rank, xend = rank, y = -0.05, yend = 0),
                 color = "black", linewidth = 0.3) +
    theme_bw() +
    labs(
      x = "Rank in Ordered Dataset",
      y = "Enrichment Score",
      title = title
    )
  
  return(p)
}


#' Compare Death Pathway Enrichment Between Groups
#'
#' Compare enrichment of death pathways between two groups.
#'
#' @param expr Expression matrix.
#' @param group Factor or character vector defining groups.
#' @param pathways Pathways to analyze.
#' @param method DE method: "t.test" or "wilcox".
#'
#' @return A data frame with comparison results.
#'
#' @export
#'
death_compare_groups <- function(expr,
                                  group,
                                  pathways = "all",
                                  method = c("t.test", "wilcox")) {
  
  method <- match.arg(method)
  
  # Get gene sets
  genesets <- get_death_geneset(pathways, type = "all")
  
  if (!is.list(genesets)) {
    genesets <- list(genesets)
    names(genesets) <- pathways
  }
  
  # Get unique groups
  groups <- unique(group)
  if (length(groups) != 2) {
    stop("Exactly 2 groups required for comparison")
  }
  
  group1 <- groups[1]
  group2 <- groups[2]
  
  results <- lapply(names(genesets), function(pathway) {
    genes <- genesets[[pathway]]
    genes <- intersect(genes, rownames(expr))
    
    if (length(genes) < 3) {
      return(NULL)
    }
    
    # Calculate mean expression for each sample
    pathway_expr <- colMeans(expr[genes, , drop = FALSE])
    
    # Split by group
    expr1 <- pathway_expr[group == group1]
    expr2 <- pathway_expr[group == group2]
    
    # Statistical test
    if (method == "t.test") {
      test <- t.test(expr2, expr1)
      p_value <- test$p.value
      stat <- test$statistic
    } else {
      test <- wilcox.test(expr2, expr1)
      p_value <- test$p.value
      stat <- test$statistic
    }
    
    # Calculate fold change
    mean1 <- mean(expr1, na.rm = TRUE)
    mean2 <- mean(expr2, na.rm = TRUE)
    log2fc <- log2(mean2 / mean1)
    
    data.frame(
      pathway = pathway,
      n_genes = length(genes),
      mean_group1 = mean1,
      mean_group2 = mean2,
      log2FC = log2fc,
      statistic = stat,
      p_value = p_value,
      stringsAsFactors = FALSE
    )
  })
  
  results <- do.call(rbind, results[!sapply(results, is.null)])
  results$p_adjust <- p.adjust(results$p_value, method = "BH")
  results <- results[order(results$p_value), ]
  
  # Add significance
  results$significance <- ifelse(results$p_adjust < 0.001, "***",
                                  ifelse(results$p_adjust < 0.01, "**",
                                         ifelse(results$p_adjust < 0.05, "*", "")))
  
  # Add direction
  results$direction <- ifelse(results$log2FC > 0, 
                               paste0("Up in ", group2),
                               paste0("Up in ", group1))
  
  rownames(results) <- NULL
  names(results)[3:4] <- c(paste0("mean_", group1), paste0("mean_", group2))
  
  class(results) <- c("death_compare", "data.frame")
  
  return(results)
}


#' Plot Volcano for Group Comparison
#'
#' Create a volcano plot for death pathway comparison between groups.
#'
#' @param compare_result A death_compare object.
#' @param fc_cutoff Fold change cutoff for significance.
#' @param p_cutoff P-value cutoff for significance.
#' @param label_top Number of top pathways to label.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @export
#'
plot_compare_volcano <- function(compare_result,
                                  fc_cutoff = 0.5,
                                  p_cutoff = 0.05,
                                  label_top = 5) {
  
  plot_data <- compare_result
  plot_data$neg_log_p <- -log10(plot_data$p_adjust)
  
  # Define significance
  plot_data$sig <- "Not Significant"
  plot_data$sig[plot_data$log2FC > fc_cutoff & plot_data$p_adjust < p_cutoff] <- "Up"
  plot_data$sig[plot_data$log2FC < -fc_cutoff & plot_data$p_adjust < p_cutoff] <- "Down"
  
  colors <- c("Up" = "#E41A1C", "Down" = "#377EB8", "Not Significant" = "grey50")
  
  p <- ggplot(plot_data, aes(x = log2FC, y = neg_log_p, color = sig)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
    scale_color_manual(values = colors) +
    theme_bw() +
    labs(
      x = "log2(Fold Change)",
      y = "-log10(Adjusted p-value)",
      title = "Death Pathway Comparison",
      color = "Significance"
    )
  
  # Add labels for top pathways
  if (requireNamespace("ggrepel", quietly = TRUE) && label_top > 0) {
    top_pathways <- head(plot_data[order(plot_data$p_value), ], label_top)
    p <- p + ggrepel::geom_text_repel(
      data = top_pathways,
      aes(label = pathway),
      size = 3,
      max.overlaps = 20
    )
  }
  
  return(p)
}
