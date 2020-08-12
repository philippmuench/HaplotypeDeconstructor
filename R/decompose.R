# # applies decompositon on x, for number of signatures r 
#' @export
nmfDecomposition <- function(x, r, includeFit = FALSE) {
  y = NMF::nmf(x, r, method ='lee')
  w = basis(y) ## signatures x k
  h = t(coef(y)) ## samples x k
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
  dc <- nmfDecomposition(x, r, includeFit = FALSE)
  
  res <- list("signatures" = dc$w, "samples" = dc$h,
              "fitted" = dc$v, "observed" = dc$m, 
              "nHapotypes" = r)
  return(res)
}



