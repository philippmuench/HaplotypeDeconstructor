

# 
#' @export
plotSNPProfile <- function(profile, row_split, title) {
  library(ComplexHeatmap)
  col_fun = circlize::colorRamp2(c(0, 0.5, 1), c("white", "#F0C06A", "#7C3233"))
  ht <- ComplexHeatmap::Heatmap(profile,  name = "AF", col = col_fun, border = TRUE,
                                cluster_rows = F, cluster_columns = F,
                                show_row_names = F,
                                column_title = title,
                                row_gap = unit(1, "mm"),column_gap = unit(1, "mm"),
                                row_split = row_split,
                                row_title_gp = gpar(col = c("black") , fontsize = 6),
                                column_title_gp = gpar(col = c("black") , fontsize = 6),
                                column_names_gp = gpar(fontsize = 6),
                                row_names_gp = gpar(fontsize = 6))
  return(ht)
}


# 
#' @export
plotHaplotypeMap <- function(sigs, cluster_columns = F, cluster_rows = F) {
  df <- data.matrix(sigs$signatures)
  col_fun = circlize::colorRamp2(c(0, max(df)/2, max(df)), c("white", "grey", "black"))
  
  pal <- wesanderson::wes_palette("Zissou1", 21, type = "continuous")
  ha <- ComplexHeatmap::Heatmap(df,
                              column_title = "Haplotype contributions",
                               name = "Contribution",
                               border = TRUE,
                               show_row_names = F,
                               row_gap = unit(1, "mm"),column_gap = unit(1, "mm"),
                               cluster_rows = cluster_rows,
                               row_title_gp = gpar(col = c("black") , fontsize = 6),
                               column_title_gp = gpar(col = c("black") , fontsize = 6),
                               column_names_gp = gpar(fontsize = 6),
                               row_names_gp = gpar(fontsize = 6),
                               show_column_dend = F,
                               show_row_dend = F,
                               cluster_columns = cluster_columns,
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
  p <- p + ggplot2::ylim(0, 1)
  p <- p + ggplot2::facet_wrap(~stat, nrow = 2, scales = "free")
  p <- p + ggplot2::theme_bw() + ggplot2::xlab("Number of Haplotyes") + ggplot2::ylab("explained variance")
  return(p)
}


# 
#' @export
plotSamples <- function(s, normalize = FALSE, remove.sample.names = F, title = "") {
  h <- data.matrix(s$samples)
  if(normalize)
    h = h / rowSums(h)
  w_df <- reshape2::melt(h, varnames = c("sample", "haplotypes"))
  w_df$haplotypes = factor(w_df$haplotypes)
  
  set.seed(42)
  palette <- randomcoloR::distinctColorPalette(s$nHapotypes)
  
  p <- ggplot2::ggplot(w_df)
  p <- p + ggplot2::geom_bar(aes_string(x = "sample", y = "value", fill = "haplotypes"),
                   color = "black", size = 0.3, stat = "identity", position = "stack")
  p <- p + ggplot2::scale_fill_manual(values = palette) 
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
  set.seed(42)
  palette <- randomcoloR::distinctColorPalette(decomposed$nHapotypes)

  p <- ggplot(w_df, aes (x = reorder(sample, day), y = value, fill = signature))
  p <- p + geom_bar( color = "black", size = 0.3, stat = "identity", position = "stack")
  p <- p +  scale_fill_manual(values = palette) + theme_minimal() + coord_flip()
  p <- p + facet_grid(group ~ ., space = "free", scales = "free")
  p <- p + xlab("") + ylab("Haplotype Contribution")
  
  return(p)
}

plotSamplesByTime <- function(s, m, normalize = T, title = "") {
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
  
  w_df <- w_df[which(!is.na(w_df$mouse)),] 
  
  p <- ggplot(w_df, aes (x = day, y = value, fill = signature))
  p <- p + geom_bar(size = 0, color = "black", stat = "identity",
                    position = "stack")
  p <- p + facet_wrap(group ~  mouse,  scales = "free", shrink = T, ncol = 3,
                      drop =T)
  p <- p + theme_bw() + theme(panel.border = element_blank(),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              axis.line = element_line(color = "black"))
  p <- p + geom_vline(xintercept = c(4, 18, 53, 67))
  p <- p + ggtitle(title)
  p <- p + scale_fill_manual(values = palette) + theme_minimal() 
  p <- p + xlab("") + ylab("Haplotype Contribution")
  p <- p + scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
  return(p)
}


# 
#' @export
plotHaplotypeAnnotation <- function(decomposed, omm_snp_annotation,
                                    hide_hyp = F, sig_threshold = 0.005){
  
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
  set.seed(42)
  palette <- randomcoloR::distinctColorPalette(decomposed$nHapotypes)
  
  p <- ggplot(all_sig, aes(x= reorder(annotation, value), y = value, fill = factor(sig_num)))
  p <- p + geom_bar(stat = "identity")
  p <- p + coord_flip() + theme_minimal() + xlab("") 
  p <- p + scale_fill_manual(values = palette) 
  p <- p + facet_wrap(. ~ sig_num)
  p <- p + theme(text = element_text(size = 6))
  
  return(p)
}

plotSamples2 <- function(s) {
  s$combined_variance <- NULL
  s2 <- reshape2::melt(s, varnames = c("sample"))
  set.seed(42)
  palette <- randomcoloR::distinctColorPalette(length(unique(s2$variable)))
  library(ggplot2)
  p <- ggplot2::ggplot(s2, aes(x = sample, y = value, fill = variable)) 
  p <- p + ggplot2::geom_bar(color = "black", size = 0.3, stat = "identity", position = "stack")
  p <- p + ggplot2::scale_fill_manual(values = palette) 
  p <- p + ggplot2::coord_flip() + ggplot2::theme_bw()
  p <- p + ggplot2::xlab("") + ggplot2::ylab("explained variance")
  p <- p + ggplot2::theme_bw() + theme(panel.border = element_blank(),
                                       panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(),
                                       axis.line = element_line(color = "black"))
  return(p)
}