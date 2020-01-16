# Clear all
rm(list = ls())

# Load the shared library
if (.Platform$OS.type == "windows"){
  # Windows
  dyn.load("likelihoodCalc.dll")

} else if (.Platform$OS.type == "unix"){
  # Unix
  dyn.load("likelihoodCalc.so")

} else {
  stop("Dynamic library can only be loaded in Windows or Unix OS")
}

# Load magittr library for data piping and dplyr for tidy data
library("magrittr")
library("dplyr")
library("stats4")

# Load the plosdata
load("optLinearBoth.RData")
#load("randomVarIn.RData")

# Source the llCalc.R function - Contains the loaded dynamic library and the called C function
# Makes llHoodFxn available
source("llCalc.R",echo = TRUE)

source("regressionFunction_LinearBoth.R",echo = TRUE)

# Load the data as tibble
plosdata <- readr::read_csv("plosData_300519.csv")

###################### Optimisation

# Starting guess informed from previous optimisation
#x00 = xxopt$par

# Run the model - This is 10 times faster than the origin Matlab - Largely because the matrix computations and sol are done in C
#out <- regressionFunction_LinearBoth(x00,plosdata)

# Start the timing
#start_time <- Sys.time()

# Run the optimiser
#xxopt <- optim(par = x00,fn = regressionFunction_LinearBoth,method = "L-BFGS-B",lower = c(0,0,rep(-Inf,15)),upper = c(Inf,2,rep(Inf,15)),control = list(trace = 10, maxit = 1000, factr = 1e-10),hessian = TRUE,plosdata = plosdata)

# Stop timing and calculate difference
#end_time <- Sys.time()
#end_time - start_time
###################### Optimisation

# Do the hypothesis test - This is the pseudo Wald's W test

# Derivative of the nested hypothesis test
H1 = matrix(0,nrow = 2,ncol = length(xxopt$par)-2)
H2 = diag(length(xxopt$par)-2)
H0 <- rbind(H1,H2)

# Compute the inverse of the observed fisher information matrix
Iinv <- solve(solve(xxopt$hessian))

# Households in data
#m = max(plosdata$HH)
m = length(xxopt$par) - 2

# Compute the hypothesis test for the PseudoMLE estimator
h <- 0
for (i in 1:length(xxopt$par)-2)
{
  # h[i] = xxopt$par[i+3] - xxopt$par[3]
  h[i] = xxopt$par[i+2]
}
h = as.matrix(h)

# Calculate the chisquare statistic
chiQ <- m * t(h) %*% solve(t(H0) %*% Iinv %*% H0) %*% h

# Compute the p-value
pval <- pchisq(chiQ,length(xxopt$par)-2,lower.tail = FALSE)

# For each co-variate, calculate the p-value
chiq = 0;
for (i in 1:(length(xxopt$par)-2))
{
  chiqScore = m * xxopt$par[i+2] * (solve(H0[,i] %*% Iinv %*% H0[,i])) * xxopt$par[i+2]

  chiq[i] = pchisq(chiqScore,1,lower.tail = FALSE)
}

# Calculate AIC
AIC = -2*(-xxopt$value) + 2*(length(xxopt$par)-2)

# Compute the confidence interval. If you are maximising the likelihood, then the covariance matrix of the estimates iS asymptotically the inverse of the negative of the Hessian
#fisher_info <- solve(new_hess)
prop_sigma <- sqrt(diag(new_hess))
upper <- xxopt$par + 1.96/prop_sigma
lower <- xxopt$par - 1.96/prop_sigma

cidata <- rbind(upper,xxopt$par,lower)
View(cidata)

# Compute the confidence interval. If you are maximising the likelihood, then the covariance matrix of the estimates iS asymptotically the inverse of the negative of the Hessian
#timD <- diag(abs(diag(xxopt$hessian)))
#ei <- eigen(xxopt$hessian)
#new_hess <- ei$vectors %*% timD %*% t(ei$vectors)
#fisher_info <- solve(new_hess)
#prop_sigma <- sqrt(diag(fisher_info))
#upper <- xxopt$par + 1.96*prop_sigma
#lower <- xxopt$par - 1.96*prop_sigma

# Unload the shared library
if (.Platform$OS.type == "windows"){
  # Windows
  dyn.unload("likelihoodCalc.dll")

} else if (.Platform$OS.type == "unix"){
  # Unix
  dyn.unload("likelihoodCalc.so")

} else {
  stop("Dynamic library can only be unloaded in Windows or Unix OS")
}

# Run the multi-level model - trial code - Random slope and random intercept
#out2 <- glmer(secondary_case_with_ILI ~ Vaccinated + (1 + Vaccinated|HH), data = plosdata, family = binomial(link = "cloglog"))

# Run random slope
#out1 <- glmer(secondary_case_with_ILI ~ Vaccinated + (1|HH), data = plosdata, family = binomial(link = "cloglog"))
