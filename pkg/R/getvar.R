#' @export
getvar = function (name, gs, burnin = 0, drop = 1) 
{
  ok <- !sapply(gs, is.null)
  if (!(length(dim(gs[[1]][[name]])) >= 2)) {
    return(sapply(gs[ok][seq(burnin + 1, sum(ok), by = drop)], 
                  function(x) {x[[name]]}))
  }
  else {
    temp <- lapply(gs[ok][seq(burnin + 1, sum(ok), by = drop)], 
                   function(x) {x[[name]]})
    return(structure(unlist(temp), dim = c(dim(temp[[1]]), length(temp))))
  }
}
