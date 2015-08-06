#' @export
corr_sp = function (setup, cor = "truncpow", cor.args = NULL) 
  {
    if (setup$numlt == 0) return(spam::diag(setup$n))
    newentries <- do.call(cor, 
                          c(list(d = setup$dmat[1, ], cutoff = setup$cutoff[1]), cor.args))
    if (setup$p > 1) 
      for (i in 2:setup$p) newentries <- newentries * do.call(cor, 
                                                              c(list(d = setup$dmat[i, ], cutoff = setup$cutoff[i]), 
                                                                cor.args))
      sp <- new("spam", entries = newentries, colindices = setup$colinds, 
                rowpointers = setup$rowpointers, dimension = rep(setup$n, 
                                                                 2))
      sp <- sp + t(sp) + spam::diag(setup$n)
      return(sp)
  }
