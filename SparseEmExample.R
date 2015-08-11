library("lhs")
library("SparseEm")

## Test function from Morris, Mitchell, and Ylvisaker (1993)

borehole <- function(x){
  rw <- x[1] * (0.15 - 0.05) + 0.05
  r <-  x[2] * (50000 - 100) + 100
  Tu <- x[3] * (115600 - 63070) + 63070
  Hu <- x[4] * (1110 - 990) + 990
  Tl <- x[5] * (116 - 63.1) + 63.1
  Hl <- x[6] * (820 - 700) + 700
  L <-  x[7] * (1680 - 1120) + 1120
  Kw <- x[8] * (12045 - 9855) + 9855
  m1 <- 2 * pi * Tu * (Hu - Hl)
  m2 <- log(r / rw)
  m3 <- 1 + 2*L*Tu/(m2*rw^2*Kw) + Tu/Tl
  return(m1/m2/m3)
}

## Set up design & simulate the data
n <- 4000 # A moderate n for illustration
npred <- 500
dim <- 8
x <- randomLHS(n+npred, dim)
y <- apply(x, 1, borehole)
ypred.0 <- y[-(1:n)]; y <- y[1:n]
xpred <- x[-(1:n),]; x <- x[1:n,]


############################################################
## Method 1: With regression terms and sparse correlation ##
############################################################

## Specify the maximum sum of polynomial exponents and maximum number
## of terms in an interaction
## Ideally would choose this using a subset or cross-validation

degree <- 2; maxint <- 2 

## Specify the desired amount of sparsity (% off-diagonal elements that = 0)

sparsity <- 0.99
## Estimate of needed cutoff based on uniform sampling
mc <- find_tau(den = 1 - sparsity, dim = ncol(x)) * ncol(x) 

## Sample from tau|Data

B <- 2000
time1 <- system.time(tau <- mcmc_sparse(y, x, mc = mc, degree = degree, maxint = maxint, B = B))
## Note: an estimate of remaining time will be printed after 10 iterations
## Also, You may safely ignore any warnings from spam that look like
## Increased 'nnzcolindices' with 'NgPeyton' method

matplot(tau, type = "l") # Trace plots from MCMC
plot(apply(tau, 1, sum), type = "l")

## Make predictions

burnin <- 500
index <- seq(burnin+1, B, by = 10) # Remove burnin; sample based on every 10th tau
time2 <- system.time(ypred.sparse <- pred_sparse(tau[index,], x, y, xpred,
                                                 degree = degree, maxint = maxint))

################################
## Method 2: "Standard model" ##
################################

## Sample from phi|Data

time3 <- system.time(phi <- mcmc_nonsparse(y, x, B = B))
## Note: a time estimate will be printed after 10 iterations
## If n is large, this can be very slow
matplot(phi, type = "l") # Trace plots from MCMC

## Make predictions

index <- seq(burnin+1, B, by = 10) # Remove burnin; sample based on every 10th tau
time4 <- system.time(ypred.nonsparse <- pred_sparse(phi[index,], x, y, xpred))

#############################
## Compare the two methods ##
#############################

## Method 1

plot(ypred.0, ypred.sparse$mean)
points(ypred.0, ypred.sparse$mean + 2 * sqrt(ypred.sparse$var), col = 2, pch = "-")
points(ypred.0, ypred.sparse$mean - 2 * sqrt(ypred.sparse$var), col = 2, pch = "-")
abline(0, 1)

1 - sum((ypred.sparse$mean - ypred.0)^2) / sum((mean(y) - ypred.0)^2) # NSE
mean(ypred.0 > ypred.sparse$mean - 2*sqrt(ypred.sparse$var) & # Coverage across the function
       ypred.0 < ypred.sparse$mean + 2*sqrt(ypred.sparse$var))

## Method 2

plot(ypred.0, ypred.nonsparse$mean)
points(ypred.0, ypred.nonsparse$mean + 2 * sqrt(ypred.nonsparse$var), col = 2, pch = "-")
points(ypred.0, ypred.nonsparse$mean - 2 * sqrt(ypred.nonsparse$var), col = 2, pch = "-")
abline(0, 1)

1 - sum((ypred.nonsparse$mean - ypred.0)^2) / sum((mean(y) - ypred.0)^2) # NSE
mean(ypred.0 > ypred.nonsparse$mean - 2*sqrt(ypred.nonsparse$var) & # Coverage across the function
       ypred.0 < ypred.nonsparse$mean + 2*sqrt(ypred.nonsparse$var))