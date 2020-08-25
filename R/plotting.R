# 
#' @export
plotHaplotypeMap <- function(sigs) {
  df <- data.matrix(sigs$signatures)
  col_fun <- circlize::colorRamp2(c(0, 0.001, 1, 5), c("white", "yellow", "red", "black"))
  
  ha <- ComplexHeatmap::Heatmap(df,
                               name = "contribution",
                               border = TRUE,
                               column_title = "Haplotype contributions",
                               show_row_names = F,
                               cluster_rows = F,
                               show_column_dend = F,
                               show_row_dend = F,
                               cluster_columns = T,
                               col = col_fun)
  
  return(ha)
}

# 
#' @export
plotNumberHaplotyes <- function(gof) {
  m <- reshape2::melt(gof, id.vars = c("NumberHaplotyes", "Replicate"),
           measure.vars = c("ExplainedVariance"), variable.name = "stat")
  p <- ggplot2::ggplot(m, ggplot2::aes_string(x = "NumberHaplotyes", y = "value", group = "NumberHaplotyes"))
  p <- p + ggplot2::stat_summary(fun.y = mean, colour = "red", size = 2.5, geom = "point")
  p <- p + ggplot2::geom_point(color = "black", shape = 3)
  p <- p + ylim(0, 1)
  p <- p + ggplot2::facet_wrap(~stat, nrow = 2, scales = "free")
  p <- p + ggplot2::theme_bw() + ggplot2::xlab("Number of Haplotyes") + ggplot2::ylab("explained variance")
  return(p)
}


# 
#' @export
plotSamples <- function(s, normalize = FALSE, remove.sample.names = F, title = "") {
  h <- data.matrix(s$samples)
 # if(normalize)
  #  h = h / rowSums(h)
  w_df <- reshape2::melt(h, varnames = c("sample", "haplotypes"))
  w_df$haplotypes = factor(w_df$haplotypes)
  
  set.seed(42)
  palette <- randomcoloR::distinctColorPalette(decomposed$nHapotypes)
  
  p <- ggplot2::ggplot(w_df)
  p <- p + ggplot2::geom_bar(aes_string(x = "sample", y = "value", fill = "haplotypes"),
                   color = "black", size = 0.3, stat = "identity", position = "stack")
 # p <- p + ggplot2::scale_fill_manual(values = palette) 
  p <- p + ggplot2::coord_flip() + ggplot2::theme_bw()
  p <- p + ggplot2::xlab("") + ggplot2::ylab("Haplotype Contribution")
  p <- p + ggplot2::theme_bw() + theme(panel.border = element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             axis.line = element_line(color = "black"))
  p <- p + ggplot2::ggtitle(title)
  p <- p + ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
  if(remove.sample.names)
    p <- p + ggplot2::theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) 

  return(p)
}

# 
#' @export
plotSample <- function(s, normalize = T, sample = "reseq 1696", title = "") {
  h <- data.matrix(s$samples)
  if(normalize)
    h = h / rowSums(h)
  w_df <- reshape2::melt(h, varnames = c("sample", "signature"))
  w_df$signature <- factor(w_df$signature)
  w_df <- w_df[which(w_df$sample == sample),]
  
  set.seed(42)
  palette <- randomcoloR::distinctColorPalette(decomposed$nHapotypes)
  
  p <- ggplot(w_df, aes (x = sample, y = value, fill = signature))
  p <- p + geom_bar(size = 0, color = "black", stat = "identity",
                    position = "stack")
  p <- p + theme_bw() + theme(panel.border = element_blank(),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              axis.line = element_line(color = "black"))
  p <- p + scale_fill_manual(values = palette) + theme_minimal() 
  p <- p + xlab("") + ylab("Haplotype Contribution")
  p <- p + ggplot2::ggtitle(title)
  p <- p + scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
  return(p)
}

