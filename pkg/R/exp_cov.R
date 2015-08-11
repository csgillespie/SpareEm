#' @export
exp_cov = function (x, phi, pow = 2)
  exp(-phi * x^pow)
