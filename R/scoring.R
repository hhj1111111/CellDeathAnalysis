#' Calculate Cell Death Pathway Scores
#'
#' Calculate enrichment scores for cell death pathways in expression data
#' using various methods including ssGSEA, GSVA, mean expression, and AUCell.
#'
#' @param expr A numeric matrix of gene expression data with genes as rows
#'   and samples as columns. Row names should be gene symbols.
#' @param pathways Character vector of pathway names to score, or "all" for
#'   all available pathways. Default is "all".
#' @param method Character. Scoring method to use:
#'   \itemize{
#'     \item "ssgsea": Single-sample GSEA (requires GSVA package)
#'     \item "gsva": Gene Set Variation Analysis (requires GSVA package)
#'     \item "mean": Mean expression of pathway genes
#'     \item "median": Median expression of pathway genes
#'     \item "aucell": AUCell method (requires AUCell package)
#'     \item "zscore": Z-score normalized mean expression
#'     \item "crosstalk": Crosstalk-aware scoring with specificity weighting
#'       and residual debiasing (see \code{\link{calculate_death_score_crosstalk}})
#'   }
#' @param min_genes Integer. Minimum number of pathway genes required in
#'   expression data. Pathways with fewer genes will be skipped. Default is 5.
#' @param scale Logical. Whether to scale scores across samples. Default is FALSE.
#' @param verbose Logical. Whether to print progress messages. Default is TRUE.
#'
#' @return A data frame with samples as rows and pathway scores as columns.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create example expression matrix
#' set.seed(123)
#' expr <- matrix(rnorm(1000 * 50), nrow = 1000, ncol = 50)
#' rownames(expr) <- paste0("Gene", 1:1000)
#' colnames(expr) <- paste0("Sample", 1:50)
#'
#' # Calculate scores using mean expression
#' scores <- calculate_death_score(expr, method = "mean")
#'
#' # Calculate scores for specific pathways
#' scores <- calculate_death_score(expr, 
#'                                  pathways = c("ferroptosis", "pyroptosis"),
#'                                  method = "zscore")
#' }
#'
calculate_death_score <- function(expr,
                                   pathways = "all",
                                   method = c("ssgsea", "gsva", "mean", "median", "aucell", "zscore", "crosstalk"),
                                   min_genes = 5,
                                   scale = FALSE,
                                   verbose = TRUE) {
  
  method <- match.arg(method)
  
  # 输入检查
  if (!is.matrix(expr) && !is.data.frame(expr)) {
    stop("expr must be a matrix or data frame")
  }
  
  if (is.data.frame(expr)) {
    expr <- as.matrix(expr)
  }
  
  if (is.null(rownames(expr))) {
    stop("expr must have gene symbols as row names")
  }
  
  # 获取基因集
  genesets <- get_death_geneset(pathways, type = "all")
  
  if (!is.list(genesets)) {
    genesets <- list(genesets)
    names(genesets) <- pathways
  }
  
  # 检查基因集与表达数据的重叠
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
  
  # 过滤基因数不足的通路
  valid_pathways <- names(overlap_info)[sapply(overlap_info, function(x) x$found >= min_genes)]
  
  if (length(valid_pathways) == 0) {
    stop("No pathways have enough genes (min_genes = ", min_genes, ") in expression data")
  }
  
  if (length(valid_pathways) < length(genesets) && verbose) {
    skipped <- setdiff(names(genesets), valid_pathways)
    warning("Skipped pathways with < ", min_genes, " genes: ", 
            paste(skipped, collapse = ", "))
  }
  
  # 准备基因集列表
  gs_list <- lapply(valid_pathways, function(p) overlap_info[[p]]$genes)
  names(gs_list) <- valid_pathways
  
  # 根据方法计算得分
  scores <- switch(method,
    "mean" = .score_mean(expr, gs_list),
    "median" = .score_median(expr, gs_list),
    "zscore" = .score_zscore(expr, gs_list),
    "ssgsea" = .score_ssgsea(expr, gs_list, verbose),
    "gsva" = .score_gsva(expr, gs_list, verbose),
    "aucell" = .score_aucell(expr, gs_list, verbose),
    "crosstalk" = {
      if (verbose) cat("Using crosstalk-aware scoring...\n")
      gene_weights <- .calc_gene_specificity(genesets)
      sw_scores <- .specificity_weighted_score(expr, gs_list, gene_weights)
      if (length(valid_pathways) >= 3) {
        .residual_debias(sw_scores)
      } else {
        if (verbose) warning("Need >= 3 pathways for debiasing. Using specificity-weighted scores only.")
        sw_scores
      }
    }
  )
  
  # 可选：标准化
  if (scale) {
    scores <- as.data.frame(scale(scores))
  }
  
  # 添加属性
  attr(scores, "method") <- method
  attr(scores, "pathways") <- valid_pathways
  attr(scores, "gene_overlap") <- overlap_info[valid_pathways]

  # 为 crosstalk 方法添加特有属性
  if (method == "crosstalk") {
    attr(scores, "gene_weights") <- gene_weights
    attr(scores, "specificity_summary") <- do.call(rbind, lapply(valid_pathways, function(p) {
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
  }
  
  class(scores) <- c("death_scores", "data.frame")
  
  return(scores)
}


# 内部函数：均值评分
.score_mean <- function(expr, gs_list) {
  scores <- sapply(gs_list, function(genes) {
    colMeans(expr[genes, , drop = FALSE], na.rm = TRUE)
  })
  as.data.frame(scores)
}

# 内部函数：中位数评分
.score_median <- function(expr, gs_list) {
  scores <- sapply(gs_list, function(genes) {
    apply(expr[genes, , drop = FALSE], 2, median, na.rm = TRUE)
  })
  as.data.frame(scores)
}

# 内部函数：Z-score评分
.score_zscore <- function(expr, gs_list) {
  # 标准化表达矩阵
  expr_scaled <- t(scale(t(expr)))
  
  scores <- sapply(gs_list, function(genes) {
    colMeans(expr_scaled[genes, , drop = FALSE], na.rm = TRUE)
  })
  as.data.frame(scores)
}

# 内部函数：ssGSEA评分
.score_ssgsea <- function(expr, gs_list, verbose) {
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("Package 'GSVA' is required for ssGSEA. ",
         "Install it with: BiocManager::install('GSVA')")
  }
  
  if (verbose) cat("Running ssGSEA...\n")
  
  # GSVA 1.46.0+ 使用新API
  if (packageVersion("GSVA") >= "1.46.0") {
    gsva_param <- GSVA::ssgseaParam(expr, gs_list)
    scores <- GSVA::gsva(gsva_param, verbose = verbose)
  } else {
    scores <- GSVA::gsva(expr, gs_list, method = "ssgsea", verbose = verbose)
  }
  
  as.data.frame(t(scores))
}

# 内部函数：GSVA评分
.score_gsva <- function(expr, gs_list, verbose) {
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("Package 'GSVA' is required for GSVA. ",
         "Install it with: BiocManager::install('GSVA')")
  }
  
  if (verbose) cat("Running GSVA...\n")
  
  if (packageVersion("GSVA") >= "1.46.0") {
    gsva_param <- GSVA::gsvaParam(expr, gs_list)
    scores <- GSVA::gsva(gsva_param, verbose = verbose)
  } else {
    scores <- GSVA::gsva(expr, gs_list, method = "gsva", verbose = verbose)
  }
  
  as.data.frame(t(scores))
}

# 内部函数：AUCell评分
.score_aucell <- function(expr, gs_list, verbose) {
  if (!requireNamespace("AUCell", quietly = TRUE)) {
    stop("Package 'AUCell' is required for AUCell. ",
         "Install it with: BiocManager::install('AUCell')")
  }
  
  if (verbose) cat("Running AUCell...\n")
  
  # 构建排名
  rankings <- AUCell::AUCell_buildRankings(expr, plotStats = FALSE, verbose = verbose)
  
  # 计算AUC
  auc <- AUCell::AUCell_calcAUC(gs_list, rankings, verbose = verbose)
  
  scores <- as.data.frame(t(AUCell::getAUC(auc)))
  
  return(scores)
}


#' Print Method for death_scores
#' @param x A death_scores object
#' @param ... Additional arguments (ignored)
#' @export
print.death_scores <- function(x, ...) {
  cat("\nCell Death Pathway Scores\n")
  cat("═════════════════════════\n")
  cat("Method:", attr(x, "method"), "\n")
  cat("Samples:", nrow(x), "\n")
  cat("Pathways:", ncol(x), "\n")
  cat("Pathways scored:", paste(names(x), collapse = ", "), "\n\n")
  
  cat("Score summary:\n")
  print(summary(x))
  
  invisible(x)
}


#' Calculate Score for Single Pathway
#'
#' A simplified function to calculate score for a single pathway.
#'
#' @param expr Expression matrix
#' @param pathway Character. Name of the pathway.
#' @param method Scoring method. Default is "zscore".
#'
#' @return A numeric vector of scores for each sample.
#'
#' @export
#'
score_pathway <- function(expr, pathway, method = "zscore") {
  scores <- calculate_death_score(expr, pathways = pathway, 
                                   method = method, verbose = FALSE)
  return(scores[[1]])
}


#' Classify Samples by Death Score
#'
#' Classify samples into high/low groups based on pathway scores.
#'
#' @param scores A death_scores object or numeric vector.
#' @param method Character. Classification method:
#'   \itemize{
#'     \item "median": Split by median value (default)
#'     \item "mean": Split by mean value
#'     \item "quantile": Split by specified quantile
#'     \item "optimal": Optimal cutpoint for survival (requires survival data)
#'   }
#' @param quantile Numeric. Quantile value for method="quantile". Default is 0.5.
#'
#' @return A factor with levels "Low" and "High".
#'
#' @export
#'
classify_by_score <- function(scores, 
                               method = c("median", "mean", "quantile", "optimal"),
                               quantile = 0.5) {
  
  method <- match.arg(method)
  
  if (is.data.frame(scores) && ncol(scores) == 1) {
    scores <- scores[[1]]
  }
  
  if (!is.numeric(scores)) {
    stop("scores must be numeric")
  }
  
  cutoff <- switch(method,
    "median" = median(scores, na.rm = TRUE),
    "mean" = mean(scores, na.rm = TRUE),
    "quantile" = quantile(scores, probs = quantile, na.rm = TRUE)
  )
  
  groups <- ifelse(scores >= cutoff, "High", "Low")
  groups <- factor(groups, levels = c("Low", "High"))
  
  attr(groups, "cutoff") <- cutoff
  attr(groups, "method") <- method
  
  return(groups)
}


#' Calculate Death Index
#'
#' Calculate a composite death index combining multiple pathways.
#'
#' @param scores A death_scores data frame.
#' @param weights Named numeric vector of weights for each pathway.
#'   If NULL, equal weights are used.
#' @param method Character. Method for combining scores:
#'   \itemize{
#'     \item "weighted_mean": Weighted mean of scores (default)
#'     \item "pca": First principal component
#'     \item "sum": Simple sum of scores
#'   }
#'
#' @return A numeric vector of composite death index values.
#'
#' @export
#'
calculate_death_index <- function(scores, 
                                   weights = NULL, 
                                   method = c("weighted_mean", "pca", "sum")) {
  
  method <- match.arg(method)
  
  if (!is.data.frame(scores)) {
    stop("scores must be a data frame")
  }
  
  score_mat <- as.matrix(scores)
  
  if (is.null(weights)) {
    weights <- rep(1, ncol(score_mat))
    names(weights) <- colnames(score_mat)
  }
  
  index <- switch(method,
    "weighted_mean" = {
      w <- weights[colnames(score_mat)]
      rowSums(score_mat * w, na.rm = TRUE) / sum(w)
    },
    "pca" = {
      pca <- prcomp(score_mat, scale. = TRUE)
      pca$x[, 1]
    },
    "sum" = {
      rowSums(score_mat, na.rm = TRUE)
    }
  )
  
  attr(index, "method") <- method
  attr(index, "weights") <- weights
  
  return(index)
}
