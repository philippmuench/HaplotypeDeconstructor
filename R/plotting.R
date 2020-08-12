# 
#' @export
plotHaplotypeMap <- function(sigs) {
  df <- data.matrix(sigs$signatures)
  ha <- ComplexHeatmap::Heatmap(df,
                               border = TRUE,
                               column_title = "Haplotype contributions",
                               show_row_names = F,
                               show_column_dend = F,
                               show_row_dend = F,
                               cluster_columns = F)
  
  return(ha)
}

# 
#' @export
plotNumberHaplotyes <- function(gof) {
  m = reshape2::melt(gof, id.vars = c("NumberHaplotyes", "Replicate"),
           measure.vars = c("ExplainedVariance"), variable.name = "stat")
  p = ggplot(m, aes_string(x = "NumberHaplotyes", y = "value", group = "NumberHaplotyes"))
  p = p + stat_summary(fun.y = mean, colour = "red", size = 2.5, geom = "point")
  p = p + geom_point(color = "black", shape = 3)
  p = p + facet_wrap(~stat, nrow = 2, scales = "free")
  p = p + theme_bw() + xlab("Number of Haplotyes") + ylab("explained variance")
  return(p)
}


# 
#' @export
plotSamples <- function(s, normalize = FALSE, percent = FALSE) {
  h <- data.matrix(s$samples)
  if(normalize) {
    h = h / rowSums(h)
    if(percent) {
      h = h * 100
    }
  }
  w_df = reshape2::melt(h, varnames = c("sample", "signature"))
  w_df$signature = factor(w_df$signature)
  
  p <- ggplot(w_df)
  p <- p + geom_bar(aes_string(x = "sample", y = "value", fill = "signature"),
                   color = "black", size = 0.3, stat = "identity", position = "stack")
  p <- p + scale_fill_brewer(palette = "Set3") + coord_flip() + theme_bw()
  p <- p + xlab("") + ylab("Haplotype Contribution")

  return(p)
}




# 
#' @export
plotSamplesByGroup <- function(s, m, normalize = FALSE, percent = FALSE) {
  h <- data.matrix(s$samples)
  if(normalize) {
    h = h / rowSums(h)
    if(percent) {
      h = h * 100
    }
  }
  w_df = reshape2::melt(h, varnames = c("sample", "signature"))
  w_df$signature = factor(w_df$signature)
  
  w_df$group <- m[match(w_df$sample, m$sample),]$group
  w_df$day <- m[match(w_df$sample, m$sample),]$day
  
  p <- ggplot(w_df, aes (x = reorder(sample, day), y = value, fill = signature))
  p <- p + geom_bar( color = "black", size = 0.3, stat = "identity", position = "stack")
  p <- p + scale_fill_brewer(palette = "Set3") + theme_minimal() + coord_flip()
  p <- p + facet_grid(group ~ ., space = "free", scales = "free")
  p <- p + xlab("") + ylab("Haplotype Contribution")
  
  return(p)
}
