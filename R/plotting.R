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
  p = p + ggplot2::stat_summary(fun.y = mean, colour = "red", size = 2.5, geom = "point")
  p = p + ggplot2::geom_point(color = "black", shape = 3)
  p = p + ggplot2::facet_wrap(~stat, nrow = 2, scales = "free")
  p = p + ggplot2::theme_bw() + ggplot2::xlab("Number of Haplotyes") + ggplot2::ylab("explained variance")
  return(p)
}


# 
#' @export
plotSamples <- function(s, normalize = FALSE, remove.sample.names =F) {
  h <- data.matrix(s$samples)
  if(normalize)
    h = h / rowSums(h)
  w_df <- reshape2::melt(h, varnames = c("sample", "haplotypes"))
  w_df$haplotypes = factor(w_df$haplotypes)
  
  set.seed(42)
  palette <- randomcoloR::distinctColorPalette(decomposed$nHapotypes)
  
  p <- ggplot2::ggplot(w_df)
  p <- p + ggplot2::geom_bar(aes_string(x = "sample", y = "value", fill = "haplotypes"),
                   color = "black", size = 0.3, stat = "identity", position = "stack")
  p <- p + ggplot2::scale_fill_manual(values = palette) + ggplot2::coord_flip() + ggplot2::theme_bw()
  p <- p + ggplot2::xlab("") + ggplot2::ylab("Haplotype Contribution")
  p <- p + ggplot2::theme_bw() + theme(panel.border = element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             axis.line = element_line(color = "black"))
  p <- p + ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
  if(remove.sample.names)
    p <- p + ggplot2::theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) 

  return(p)
}

# 
#' @export
plotSamplesByGroup <- function(s, m, normalize = FALSE) {
  h <- data.matrix(s$samples)
  if(normalize)
    h = h / rowSums(h)
  w_df <- reshape2::melt(h, varnames = c("sample", "signature"))
  w_df$signature <- factor(w_df$signature)
  
  set.seed(42)
  palette <- randomcoloR::distinctColorPalette(decomposed$nHapotypes)

  w_df$group <- m[match(w_df$sample, m$sample),]$group
  w_df$day <- m[match(w_df$sample, m$sample),]$day
  w_df$study <- m[match(w_df$sample, m$sample),]$study
  w_df$mouse <- m[match(w_df$sample, m$sample),]$mouse
  
  
  p <- ggplot(w_df, aes (x = reorder(sample, day), y = value, fill = signature))
  p <- p + geom_bar(size = 0, color = "black", stat = "identity",
                    position = "stack")
  p <- p + scale_fill_manual(values = palette) + theme_minimal() + coord_flip()
  p <- p + facet_wrap(~ study + group,  scales = "free_x", ncol = 2, shrink = T,
                      drop =T, strip.position = "left")
  p <- p + xlab("") + ylab("Haplotype Contribution")
  p <- p + theme_bw() + theme(panel.border = element_blank(),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              axis.line = element_line(color = "black"))
  p <- p + scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
  
  return(p)
}

# 
#' @export
plotHaplotypeAnnotation <- function(decomposed, omm_snp_annotation,
                                    hide_hyp = T, sig_threshold = 0.005){
  
  sig <- data.matrix(decomposed$signatures)
  res <- list()
  for (i in seq_along(1:ncol(sig))){
    message(i)
    one_sig <- as.data.frame(sig[,i])
    one_sig$annotation <- omm_snp_annotation[match(rownames(one_sig),
                                                   omm_snp_annotation$id),]$description
    if (hide_hyp){
      one_sig <- one_sig[which(one_sig$annotation != "outside ORFs" &
                                 one_sig$annotation != "hypothetical protein"),]
    }
    
    colnames(one_sig) <-c("value", "annotation")
    one_sig <- one_sig[order(-one_sig$value),] 
    one_sig$id <- rownames(one_sig)
    one_sig <- one_sig[which(one_sig$value > sig_threshold),]
    one_sig$sig_num <- i
    res[[i]] <- one_sig
  }
  
  all_sig <- do.call(rbind, res)
 
  p <- ggplot(all_sig, aes(x= reorder(annotation, value), y = value, fill = factor(sig_num)))
  p <- p + geom_bar(stat = "identity")
  p <- p + coord_flip() + theme_minimal() + xlab("") 
  p <- p + scale_fill_brewer(palette = "Set3") 
  p <- p + facet_wrap(. ~ sig_num)
  
  return(p)
}
