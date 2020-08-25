# 
#' @export
assessOneSignature <- function(m, n) {
  sigs <- findHaplotypes(m, n)
  fitted <- data.matrix(sigs$fitted)
  observed <- data.matrix(sigs$observed)
  exvar <- evar(fitted, observed)
  gof <- data.frame("NumberHaplotyes" = n,
    "Observed" = observed,
    "ExplainedVariance" = exvar)
  return(gof)
}

# 
#' @export
assessNumberHaplotyes <- function(m, nHaps, nReplicates = 1) {
  ## compute fit statistics for
  ## - outer :: number of signatures
  ## - inner :: replicates
  dev = lapply(nHaps, function(r, m) {
    d = lapply(1:nReplicates, function(i) {
      dev = assessOneSignature(m, r)
      dev$Replicate = i
      dev
    })
    return(do.call(rbind, d))
  }, m)
  
  ## merge results to data frame
  gof = do.call(rbind, dev)
  return(gof)
}
