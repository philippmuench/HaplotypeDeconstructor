nmfDecomposition <- function(x, r, includeFit = FALSE) {
  y = NMF::nmf(x, r, method= "snmf/l", beta = 1e-3)
  w = NMF::basis(y) ## signatures x k
  h = t(NMF::coef(y)) ## samples x k
  ## order signatures
  ord = order(rowMax(t(w)), decreasing = TRUE)
  w = w[ ,ord]
  h = h[ ,ord]
  ## name signatures: S1, ..., Sn
  sig_names = paste0("H", 1:r)
  colnames(w) = colnames(h) = sig_names
  v = fitted(y)
  res = list(w = w, h = h, v = v, m = x, r = r)
  if(includeFit)
    res[["raw"]] = y
  return(res)
}

# # applies decompositon on x, for number of signatures r 
#' @export
findHaplotypes <- function(x, r) {
  x <- data.matrix(x)
  dc <- nmfDecomposition(x, r, includeFit = T)
  res <- list("signatures" = dc$w, "samples" = dc$h,
              "fitted" = dc$v, "observed" = dc$m, 
              "nHapotypes" = r)
  return(res)
}

# if ratio=T, the individual haplotype contributions are calculated as fraction to the total haplotype contribution 
HaplotypeEvar <- function(decomposed){
  res_per_sample_all <- list() #stores the % explained variance of all haplotypes on a sample
  observed <- as.data.frame(decomposed$observed)
  signatures <- as.data.frame(decomposed$signatures)
  h <- data.matrix(decomposed$samples)
  h <- h / rowSums(h)
  it <- 1
  # reconstruct using all haplotypes
  rec <- as.matrix(as.data.frame(decomposed$samples)) %*% t(as.matrix(signatures))
  # calculate total haplotype contribution
  NMF::evar(as.matrix(rec), as.matrix(t(observed)))
  it <- 1
  for (s in 1:ncol(decomposed$observed)){
    message(paste0("Processing Sample"), s)
    res_per_sample_all[[it]] <- data.frame(sample= colnames(decomposed$observed)[s], missing_variance = 1- NMF::evar(as.matrix(rec[s,]),
                                                                                                                     as.matrix(observed[, s])),
                                           combined_variance = NMF::evar(as.matrix(rec[s,]),
                                                                         as.matrix(observed[, s])))
    it <- it + 1
  }
  res <- do.call("rbind", res_per_sample_all)
  norm_h <- h * 100 * res$combined_variance
  results <- cbind(res, norm_h)
  results$missing_variance <- 100 * results$missing_variance
  return(results)
}