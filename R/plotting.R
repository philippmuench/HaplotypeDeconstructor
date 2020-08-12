plotHaplotypeMap <- function(sigs) {
  df <- data.matrix(sigs$signatures)
  ha <- ComplexHeatmap::Heatmap(df,
                               border = TRUE,
                               column_title = "Haplotype contributions",
                               show_row_names = F,
                               show_column_dend = F,
                               show_row_dend = F)
  
  return(ha)
}

plotNumberHaplotyes <- function(gof) {
  m = melt(gof, id.vars = c("NumberHaplotyes", "Replicate"),
           measure.vars = c("ExplainedVariance"), variable.name = "stat")
  p = ggplot(m, aes_string(x = "NumberHaplotyes", y = "value", group = "NumberHaplotyes"))
  p = p + stat_summary(fun.y = mean, colour = "red", size = 2.5, geom = "point")
  p = p + geom_point(color = "black", shape = 3)
  p = p + facet_wrap(~stat, nrow = 2, scales = "free")
  p = p + theme_bw() + xlab("Number of Haplotyes") + ylab("Statistic")
  return(p)
}