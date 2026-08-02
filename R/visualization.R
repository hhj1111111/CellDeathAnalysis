#' Plot Death Score Heatmap
#'
#' Create a heatmap visualization of cell death pathway scores.
#'
#' @param scores A death_scores object or matrix/data frame.
#' @param annotation_col Data frame of sample annotations for column annotation.
#' @param cluster_rows Logical. Whether to cluster pathways. Default is TRUE.
#' @param cluster_cols Logical. Whether to cluster samples. Default is TRUE.
#' @param scale Character. Scale by "row", "column", or "none". Default is "row".
#' @param colors Color palette. Default uses blue-white-red.
#' @param show_colnames Logical. Whether to show sample names. Default is FALSE.
#' @param ... Additional arguments passed to ComplexHeatmap::Heatmap or pheatmap.
#'
#' @return A heatmap object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic heatmap
#' plot_death_heatmap(scores)
#'
#' # With sample annotations
#' anno <- data.frame(Group = c(rep("A", 25), rep("B", 25)))
#' plot_death_heatmap(scores, annotation_col = anno)
#' }
#'
plot_death_heatmap <- function(scores,
                                annotation_col = NULL,
                                cluster_rows = TRUE,
                                cluster_cols = TRUE,
                                scale = c("row", "column", "none"),
                                colors = NULL,
                                show_colnames = FALSE,
                                ...) {
  
  scale <- match.arg(scale)
  
  # 转换为矩阵
  if (is.data.frame(scores)) {
    score_mat <- as.matrix(scores)
  } else {
    score_mat <- scores
  }
  
  # 转置使pathway为行
  score_mat <- t(score_mat)
  
  # 设置颜色
  if (is.null(colors)) {
    colors <- grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  }
  
  # 尝试使用ComplexHeatmap，否则使用pheatmap
  if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    
    # 准备注释
    if (!is.null(annotation_col)) {
      ha <- ComplexHeatmap::HeatmapAnnotation(
        df = annotation_col,
        show_legend = TRUE
      )
    } else {
      ha <- NULL
    }
    
    # 创建热图
    ht <- ComplexHeatmap::Heatmap(
      score_mat,
      name = "Score",
      col = colors,
      cluster_rows = cluster_rows,
      cluster_columns = cluster_cols,
      show_column_names = show_colnames,
      top_annotation = ha,
      row_names_gp = grid::gpar(fontsize = 10),
      heatmap_legend_param = list(title = "Score"),
      ...
    )
    
    return(ht)
    
  } else if (requireNamespace("pheatmap", quietly = TRUE)) {
    
    pheatmap::pheatmap(
      score_mat,
      scale = scale,
      cluster_rows = cluster_rows,
      cluster_cols = cluster_cols,
      show_colnames = show_colnames,
      annotation_col = annotation_col,
      color = colors,
      ...
    )
    
  } else {
    stop("Please install 'ComplexHeatmap' or 'pheatmap' for heatmap visualization")
  }
}


#' Plot Death Score Boxplot
#'
#' Create a boxplot comparing death scores between groups.
#'
#' @param scores A death_scores object or data frame.
#' @param group A factor or character vector defining groups.
#' @param pathways Character vector of pathways to plot. Default is all.
#' @param colors Named vector of colors for groups.
#' @param add_points Logical. Whether to add individual points. Default is TRUE.
#' @param add_pvalue Logical. Whether to add significance test p-values. Default is TRUE.
#' @param test Character. Statistical test: "wilcox" or "t.test". Default is "wilcox".
#' @param ncol Integer. Number of columns for facet wrap. Default is 3.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @export
#'
plot_death_boxplot <- function(scores,
                                group,
                                pathways = NULL,
                                colors = NULL,
                                add_points = TRUE,
                                add_pvalue = TRUE,
                                test = c("wilcox", "t.test"),
                                ncol = 3) {
  
  test <- match.arg(test)
  
  # 准备数据
  if (is.null(pathways)) {
    pathways <- names(scores)
  }
  
  plot_data <- tidyr::pivot_longer(
    cbind(scores[, pathways, drop = FALSE], Group = group),
    cols = -Group,
    names_to = "Pathway",
    values_to = "Score"
  )
  
  # 设置颜色
  if (is.null(colors)) {
    n_groups <- length(unique(group))
    colors <- scales::hue_pal()(n_groups)
    names(colors) <- unique(group)
  }
  
  # 基础图
  p <- ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    facet_wrap(~ Pathway, scales = "free_y", ncol = ncol) +
    scale_fill_manual(values = colors) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    labs(x = NULL, y = "Pathway Score")
  
  # 添加点
  if (add_points) {
    p <- p + geom_jitter(width = 0.2, alpha = 0.5, size = 1)
  }
  
  # 添加p值
  if (add_pvalue && requireNamespace("ggpubr", quietly = TRUE)) {
    p <- p + ggpubr::stat_compare_means(
      method = test,
      label = "p.format",
      label.x.npc = "center",
      label.y.npc = "top"
    )
  }
  
  return(p)
}


#' Plot Death Score Radar Chart
#'
#' Create a radar (spider) chart for pathway scores.
#'
#' @param scores A death_scores object, or a named numeric vector for single sample.
#' @param group Optional grouping variable for multiple samples.
#' @param pathways Character vector of pathways to include.
#' @param colors Colors for groups or samples.
#' @param fill_alpha Numeric. Alpha for fill. Default is 0.3.
#' @param line_width Numeric. Line width. Default is 1.5.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @export
#'
plot_death_radar <- function(scores,
                              group = NULL,
                              pathways = NULL,
                              colors = NULL,
                              fill_alpha = 0.3,
                              line_width = 1.5) {
  
  # 准备数据
  if (is.vector(scores) && !is.list(scores)) {
    # 单个样本
    if (is.null(pathways)) pathways <- names(scores)
    
    plot_data <- data.frame(
      Pathway = factor(pathways, levels = pathways),
      Score = scores[pathways],
      Group = "Sample"
    )
  } else {
    # 多个样本
    if (is.null(pathways)) pathways <- names(scores)
    
    if (!is.null(group)) {
      # 按组计算均值
      score_mat <- as.matrix(scores[, pathways])
      means <- aggregate(score_mat, by = list(Group = group), FUN = mean)
      
      plot_data <- tidyr::pivot_longer(
        means,
        cols = -Group,
        names_to = "Pathway",
        values_to = "Score"
      )
    } else {
      # 计算总体均值
      plot_data <- data.frame(
        Pathway = pathways,
        Score = colMeans(scores[, pathways]),
        Group = "Mean"
      )
    }
    
    plot_data$Pathway <- factor(plot_data$Pathway, levels = pathways)
  }
  
  # 设置颜色
  if (is.null(colors)) {
    n_groups <- length(unique(plot_data$Group))
    colors <- scales::hue_pal()(n_groups)
  }
  
  # 创建雷达图 (使用coord_polar)
  p <- ggplot(plot_data, aes(x = Pathway, y = Score, group = Group, 
                              color = Group, fill = Group)) +
    geom_polygon(alpha = fill_alpha) +
    geom_line(linewidth = line_width) +
    geom_point(size = 3) +
    coord_polar() +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      panel.grid.major = element_line(color = "grey80"),
      legend.position = "bottom"
    )
  
  return(p)
}


#' Plot Pathway Correlation
#'
#' Create a correlation heatmap between cell death pathways.
#'
#' @param scores A death_scores object or data frame.
#' @param method Correlation method: "pearson", "spearman", or "kendall".
#' @param colors Color palette for correlation values.
#' @param show_values Logical. Whether to show correlation values. Default is TRUE.
#' @param ... Additional arguments passed to heatmap function.
#'
#' @return A heatmap object.
#'
#' @export
#'
plot_pathway_correlation <- function(scores,
                                      method = c("pearson", "spearman", "kendall"),
                                      colors = NULL,
                                      show_values = TRUE,
                                      ...) {
  
  method <- match.arg(method)
  
  # 计算相关性
  cor_mat <- cor(scores, method = method, use = "pairwise.complete.obs")
  
  # 设置颜色
  if (is.null(colors)) {
    colors <- grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  }
  
  if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    
    # 添加相关系数文本
    cell_fun <- NULL
    if (show_values) {
      cell_fun <- function(j, i, x, y, width, height, fill) {
        grid::grid.text(sprintf("%.2f", cor_mat[i, j]), x, y, 
                        gp = grid::gpar(fontsize = 8))
      }
    }
    
    ht <- ComplexHeatmap::Heatmap(
      cor_mat,
      name = paste0(method, "\ncorrelation"),
      col = circlize::colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
      cell_fun = cell_fun,
      row_names_gp = grid::gpar(fontsize = 10),
      column_names_gp = grid::gpar(fontsize = 10),
      ...
    )
    
    return(ht)
    
  } else if (requireNamespace("pheatmap", quietly = TRUE)) {
    
    pheatmap::pheatmap(
      cor_mat,
      color = colors,
      display_numbers = show_values,
      number_format = "%.2f",
      ...
    )
    
  } else {
    # 基础R热图
    heatmap(cor_mat, col = colors, ...)
  }
}


#' Plot Gene Expression in Pathway
#'
#' Visualize expression of genes in a specific pathway.
#'
#' @param expr Expression matrix.
#' @param pathway Character. Name of the pathway.
#' @param group Optional grouping variable.
#' @param top_n Integer. Number of top variable genes to show. Default is 20.
#' @param scale Logical. Whether to scale expression. Default is TRUE.
#' @param ... Additional arguments passed to heatmap function.
#'
#' @return A heatmap object.
#'
#' @export
#'
plot_pathway_genes <- function(expr,species = 'human',
                                pathway,
                                group = NULL,
                                top_n = 20,
                                scale = TRUE,
                                ...) {
  
  # 获取通路基因
  genes <- get_death_geneset(pathway, type = "all")
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


  # 找到表达数据中存在的基因
  common_genes <- intersect(genes, rownames(expr))
  
  if (length(common_genes) == 0) {
    stop("No pathway genes found in expression data")
  }
  
  # 选择top变异基因
  if (length(common_genes) > top_n) {
    gene_var <- apply(expr[common_genes, ], 1, var, na.rm = TRUE)
    common_genes <- names(sort(gene_var, decreasing = TRUE))[1:top_n]
  }
  
  # 提取表达矩阵
  plot_mat <- expr[common_genes, , drop = FALSE]
  
  # 标准化
  if (scale) {
    plot_mat <- t(scale(t(plot_mat)))
  }
  
  # 准备注释
  annotation_col <- NULL
  if (!is.null(group)) {
    annotation_col <- data.frame(Group = group)
    rownames(annotation_col) <- colnames(plot_mat)
  }
  
  # 绘制热图
  plot_death_heatmap(
    t(plot_mat),
    annotation_col = annotation_col,
    scale = "none",
    ...
  )
}


#' Plot Death Score Distribution
#'
#' Visualize the distribution of pathway scores.
#'
#' @param scores A death_scores object or data frame.
#' @param pathways Character vector of pathways to plot.
#' @param type Character. Plot type: "density", "histogram", "violin", or "ridge".
#' @param group Optional grouping variable.
#' @param colors Colors for groups.
#' @param ncol Integer. Number of columns for facets.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @export
#'
plot_score_distribution <- function(scores,
                                     pathways = NULL,
                                     type = c("density", "histogram", "violin", "ridge"),
                                     group = NULL,
                                     colors = NULL,
                                     ncol = 3) {
  
  type <- match.arg(type)
  
  if (is.null(pathways)) {
    pathways <- names(scores)
  }
  
  # 准备数据
  if (is.null(group)) {
    plot_data <- tidyr::pivot_longer(
      scores[, pathways, drop = FALSE],
      cols = everything(),
      names_to = "Pathway",
      values_to = "Score"
    )
    plot_data$Group <- "All"
  } else {
    plot_data <- tidyr::pivot_longer(
      cbind(scores[, pathways, drop = FALSE], Group = group),
      cols = -Group,
      names_to = "Pathway",
      values_to = "Score"
    )
  }
  
  # 设置颜色
  if (is.null(colors)) {
    n_groups <- length(unique(plot_data$Group))
    colors <- scales::hue_pal()(n_groups)
  }
  
  # 根据类型绘图
  p <- ggplot(plot_data, aes(x = Score, fill = Group, color = Group))
  
  if (type == "density") {
    p <- p + 
      geom_density(alpha = 0.5) +
      facet_wrap(~ Pathway, scales = "free", ncol = ncol)
      
  } else if (type == "histogram") {
    p <- p + 
      geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
      facet_wrap(~ Pathway, scales = "free", ncol = ncol)
      
  } else if (type == "violin") {
    p <- ggplot(plot_data, aes(x = Pathway, y = Score, fill = Group)) +
      geom_violin(alpha = 0.7) +
      geom_boxplot(width = 0.1, position = position_dodge(0.9))
      
  } else if (type == "ridge") {
    if (!requireNamespace("ggridges", quietly = TRUE)) {
      stop("Package 'ggridges' required for ridge plots")
    }
    p <- ggplot(plot_data, aes(x = Score, y = Pathway, fill = Group)) +
      ggridges::geom_density_ridges(alpha = 0.5)
  }
  
  p <- p +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90"),
      legend.position = "bottom"
    )
  
  return(p)
}


#' Plot Pathway Network
#'
#' Create a network visualization of pathway relationships.
#'
#' @param scores A death_scores object.
#' @param cor_threshold Numeric. Minimum correlation to show edge. Default is 0.3.
#' @param layout Character. Network layout algorithm. Default is "fr".
#'
#' @return A ggplot object.
#'
#' @export
#'
plot_pathway_network <- function(scores,
                                  cor_threshold = 0.3,
                                  layout = "fr") {
  
  if (!requireNamespace("igraph", quietly = TRUE) ||
      !requireNamespace("ggraph", quietly = TRUE)) {
    stop("Packages 'igraph' and 'ggraph' are required for network plots")
  }
  
  # 计算相关性
  cor_mat <- cor(scores, method = "spearman")
  
  # 创建边列表
  edges <- which(abs(cor_mat) >= cor_threshold & upper.tri(cor_mat), arr.ind = TRUE)
  
  if (nrow(edges) == 0) {
    stop("No edges with correlation >= ", cor_threshold)
  }
  
  edge_df <- data.frame(
    from = rownames(cor_mat)[edges[, 1]],
    to = colnames(cor_mat)[edges[, 2]],
    weight = cor_mat[edges]
  )
  
  # 创建图
  g <- igraph::graph_from_data_frame(edge_df, directed = FALSE)
  
  # 绘制网络
  p <- ggraph::ggraph(g, layout = layout) +
    ggraph::geom_edge_link(aes(edge_alpha = abs(weight), 
                                edge_color = weight > 0),
                           edge_width = 1) +
    ggraph::geom_node_point(size = 8, color = "#3498db") +
    ggraph::geom_node_text(aes(label = name), repel = TRUE) +
    ggraph::scale_edge_color_manual(values = c("TRUE" = "#e74c3c", "FALSE" = "#2980b9"),
                                     name = "Correlation") +
    ggraph::theme_graph() +
    labs(title = "Cell Death Pathway Network")
  
  return(p)
}
