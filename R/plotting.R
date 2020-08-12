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

# 
#' @export
plotHaplotypeAnnotation <- function(decomposed, omm_snp_annotation, sig_id = 2,
                                    hide_hyp = T, sig_threshold = 5){
  
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
