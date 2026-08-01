#' =============================================================================
#' Crosstalk-Aware Pathway Scoring Module for CellDeathAnalysis
#' =============================================================================
#'
#' Novel scoring algorithm that accounts for gene overlap between cell death
#' pathways. Uses specificity weighting (IDF-inspired) and residual debiasing
#' to produce pathway scores that reflect independent pathway activity.
#'
#' @author Keran Sun
#' =============================================================================


#' Calculate Gene Specificity Weights
#'
#' Compute specificity weights for genes based on how many pathways they belong
#' to. Genes that appear in fewer pathways receive higher weights (IDF-inspired).
#'
#' @param genesets Named list of character vectors (pathway -> gene symbols)
#'
#' @return A named numeric vector of specificity weights for all unique genes
#'
#' @keywords internal
#'
.calc_gene_specificity <- function(genesets) {

  # Count how many pathways each gene belongs to
  all_genes <- unique(unlist(genesets))
  n_pathways <- length(genesets)

  gene_pathway_count <- setNames(
    sapply(all_genes, function(g) {
      sum(sapply(genesets, function(gs) g %in% gs))
    }),
    all_genes
  )

  # IDF-inspired weight: log2(N / n_g) + 1
  # Genes in 1 pathway  -> log2(14/1) + 1 = 4.81 (high weight)
  # Genes in 14 pathways -> log2(14/14) + 1 = 1.00 (low weight)
  weights <- log2(n_pathways / gene_pathway_count) + 1

  return(weights)
}


#' Specificity-Weighted Pathway Scoring
#'
#' Calculate pathway scores using specificity-weighted mean expression.
#' Genes shared across many pathways contribute less to any single pathway's
#' score, reducing redundancy from inter-pathway gene overlap.
#'
#' @param expr Numeric expression matrix (genes as rows, samples as columns)
#' @param gs_list Named list of gene sets (pathway -> character vector of genes)
#' @param gene_weights Named numeric vector of gene specificity weights
#'
#' @return A data.frame with samples as rows and pathways as columns
#'
#' @keywords internal
#'
.specificity_weighted_score <- function(expr, gs_list, gene_weights) {

  score_mat <- sapply(gs_list, function(genes) {
    # Get genes present in expression data
    valid_genes <- intersect(genes, rownames(expr))

    if (length(valid_genes) == 0) {
      return(rep(NA, ncol(expr)))
    }

    # Get weights for these genes
    w <- gene_weights[valid_genes]
    w[is.na(w)] <- 1  # fallback for genes not in weight vector

    # Weighted mean: sum(w_i * x_i) / sum(w_i)
    expr_sub <- expr[valid_genes, , drop = FALSE]
    weighted_sum <- colSums(expr_sub * w, na.rm = TRUE)
    weight_total <- sum(w, na.rm = TRUE)

    return(weighted_sum / weight_total)
  })

  as.data.frame(score_mat)
}


#' Residual Debiasing of Pathway Scores
#'
#' Remove cross-pathway confounding by regressing each pathway's score on the
#' mean score of all other pathways, and using the residual as the
#' debiased score.
#'
#' @param sw_scores Data.frame of specificity-weighted scores (samples x pathways)
#'
#' @return A data.frame of residual-debiased scores (same dimensions)
#'
#' @keywords internal
#'
.residual_debias <- function(sw_scores) {

  pathway_names <- names(sw_scores)
  n_pathways <- length(pathway_names)

  debiased <- as.data.frame(lapply(pathway_names, function(p) {
    # Score for current pathway
    y <- sw_scores[[p]]

    # Mean score of all other pathways (confounding variable)
    others <- setdiff(pathway_names, p)
    x <- rowMeans(sw_scores[, others, drop = FALSE], na.rm = TRUE)

    # Linear regression: y ~ x
    fit <- lm(y ~ x)

    # Residual = pathway-specific signal after removing cross-pathway effect
    residuals(fit)
  }))

  names(debiased) <- pathway_names
  rownames(debiased) <- rownames(sw_scores)

  return(debiased)
}


#' Calculate Crosstalk-Aware Pathway Scores
#'
#' A novel scoring method that accounts for gene overlap between cell death
#' pathways. Uses two key innovations:
#'
#' \enumerate{
#'   \item \strong{Specificity weighting}: Genes that appear in fewer pathways
#'   receive higher weights (IDF-inspired). This reduces the influence of
#'   shared genes like CASP8 (in 4 pathways) relative to pathway-specific
#'   genes like GPX4 (ferroptosis only).
#'
#'   \item \strong{Residual debiasing}: After specificity weighting, each
#'   pathway's score is regressed on the mean of all other pathways' scores.
#'   The residual represents the pathway's independent activity, removing
#'   confounding from broadly active cell death processes.
#' }
#'
#' @param expr A numeric matrix of gene expression data with genes as rows
#'   and samples as columns. Row names should be gene symbols.
#' @param pathways Character vector of pathway names to score, or "all" for
#'   all available pathways. Default is "all".
#' @param min_genes Integer. Minimum number of pathway genes required in
#'   expression data. Pathways with fewer genes will be skipped. Default is 5.
#' @param debias Logical. Whether to apply residual debiasing step.
#'   Default is TRUE. Set to FALSE to use only specificity weighting.
#' @param verbose Logical. Whether to print progress messages. Default is TRUE.
#'
#' @return A data frame (with class "death_scores") with samples as rows and
#'   pathway scores as columns. Additional attributes:
#'   \itemize{
#'     \item method: "crosstalk"
#'     \item pathways: names of scored pathways
#'     \item gene_overlap: overlap info for each pathway
#'     \item gene_weights: the specificity weights used
#'     \item specificity_summary: summary of weight distribution per pathway
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(example_expr)
#'
#' # Calculate crosstalk-aware scores
#' scores_ct <- calculate_death_score_crosstalk(example_expr)
#'
#' # Compare with standard z-score
#' scores_z <- calculate_death_score(example_expr, method = "zscore")
#'
#' # Check correlation between methods
#' cor(scores_ct$ferroptosis, scores_z$ferroptosis)
#' }
#'
calculate_death_score_crosstalk <- function(expr,species = 'human',
                                             pathways = "all",
                                             min_genes = 5,
                                             debias = TRUE,
                                             verbose = TRUE) {

  # -------------------------------------------------------------------------
  # Input validation
  # -------------------------------------------------------------------------
  if (!is.matrix(expr) && !is.data.frame(expr)) {
    stop("expr must be a matrix or data frame")
  }

  if (is.data.frame(expr)) {
    expr <- as.matrix(expr)
  }

  if (is.null(rownames(expr))) {
    stop("expr must have gene symbols as row names")
  }

  # -------------------------------------------------------------------------
  # Get gene sets
  # -------------------------------------------------------------------------
  if(sp='human')
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

  # -------------------------------------------------------------------------
  # Check overlap with expression data
  # -------------------------------------------------------------------------
  expr_genes <- rownames(expr)

  overlap_info <- lapply(genesets, function(gs) {
    common <- intersect(gs, expr_genes)
    list(
      total = length(gs),
      found = length(common),
      genes = common
    )
  })

  if (verbose) {
    cat("Gene set overlap with expression data:\n")
    for (name in names(overlap_info)) {
      info <- overlap_info[[name]]
      cat(sprintf("  %s: %d/%d genes (%.1f%%)\n",
                  name, info$found, info$total,
                  100 * info$found / info$total))
    }
    cat("\n")
  }

  # Filter pathways with insufficient genes
  valid_pathways <- names(overlap_info)[
    sapply(overlap_info, function(x) x$found >= min_genes)
  ]

  if (length(valid_pathways) == 0) {
    stop("No pathways have enough genes (min_genes = ", min_genes,
         ") in expression data")
  }

  if (length(valid_pathways) < length(genesets) && verbose) {
    skipped <- setdiff(names(genesets), valid_pathways)
    warning("Skipped pathways with < ", min_genes, " genes: ",
            paste(skipped, collapse = ", "))
  }

  # Prepare gene set list with only genes found in expression data
  gs_list <- lapply(valid_pathways, function(p) overlap_info[[p]]$genes)
  names(gs_list) <- valid_pathways

  # -------------------------------------------------------------------------
  # Step 1: Calculate gene specificity weights
  # -------------------------------------------------------------------------
  if (verbose) cat("Step 1: Calculating gene specificity weights...\n")

  gene_weights <- .calc_gene_specificity(genesets)

  # Show weight summary
  if (verbose) {
    pathway_weight_summary <- lapply(valid_pathways, function(p) {
      genes <- gs_list[[p]]
      w <- gene_weights[genes]
      w[is.na(w)] <- 1
      list(
        pathway = p,
        mean_weight = mean(w),
        min_weight = min(w),
        max_weight = max(w),
        n_unique = sum(w > 3),      # genes in <= 2 pathways
        n_shared = sum(w <= 2)      # genes in >= 4 pathways
      )
    })

    cat("  Specificity weight summary:\n")
    for (ws in pathway_weight_summary) {
      cat(sprintf("    %s: mean=%.2f, unique_genes=%d, shared_genes=%d\n",
                  ws$pathway, ws$mean_weight, ws$n_unique, ws$n_shared))
    }
    cat("\n")
  }

  # -------------------------------------------------------------------------
  # Step 2: Specificity-weighted scoring
  # -------------------------------------------------------------------------
  if (verbose) cat("Step 2: Computing specificity-weighted scores...\n")

  sw_scores <- .specificity_weighted_score(expr, gs_list, gene_weights)

  # -------------------------------------------------------------------------
  # Step 3: Residual debiasing (optional)
  # -------------------------------------------------------------------------
  if (debias) {
    if (verbose) cat("Step 3: Applying residual debiasing...\n")

    # Check that we have enough pathways for debiasing
    if (length(valid_pathways) < 3) {
      warning("Need at least 3 pathways for residual debiasing. ",
              "Skipping debiasing step.")
      debiased_scores <- sw_scores
    } else {
      debiased_scores <- .residual_debias(sw_scores)
    }
  } else {
    debiased_scores <- sw_scores
  }

  # -------------------------------------------------------------------------
  # Build result
  # -------------------------------------------------------------------------
  scores <- debiased_scores

  # Compute specificity summary for output
  specificity_summary <- do.call(rbind, lapply(valid_pathways, function(p) {
    genes <- gs_list[[p]]
    w <- gene_weights[genes]
    w[is.na(w)] <- 1
    data.frame(
      pathway = p,
      n_genes = length(genes),
      mean_specificity = mean(w),
      n_pathway_specific = sum(w > 3),
      n_highly_shared = sum(w <= 2),
      stringsAsFactors = FALSE
    )
  }))

  # Set attributes
  attr(scores, "method") <- "crosstalk"
  attr(scores, "pathways") <- valid_pathways
  attr(scores, "gene_overlap") <- overlap_info[valid_pathways]
  attr(scores, "gene_weights") <- gene_weights
  attr(scores, "specificity_summary") <- specificity_summary
  attr(scores, "debias") <- debias

  class(scores) <- c("death_scores", "data.frame")

  if (verbose) {
    cat("\nCrosstalk-aware scoring complete!\n")
    cat(sprintf("  Pathways scored: %d\n", length(valid_pathways)))
    cat(sprintf("  Samples: %d\n", nrow(scores)))
    cat(sprintf("  Debiasing: %s\n", ifelse(debias, "Applied", "Not applied")))
  }

  return(scores)
}
