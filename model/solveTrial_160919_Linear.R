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
load("optLinear.RData")
#load("randomVarIn.RData")

# Source the llCalc.R function - Contains the loaded dynamic library and the called C function
# Makes llHoodFxn available
source("llCalc.R",echo = TRUE)

source("regressionFunction_Linear.R",echo = TRUE)

# Load the data as tibble
plosdata <- readr::read_csv("plosData_300519.csv")

# Run the model - This is 10 times faster than the origin Matlab - Largely because the matrix computations and sol are done in C
#xx0 <- xxopt$par
#out <- regressionFunction_Linear(xx0,plosdata)

############### Do the optimisation
# Start the timing
#start_time <- Sys.time()

# Run the optimiser
#xxopt <- optim(par = x0,fn = regressionFunction_Linear,method = "L-BFGS-B",lower = c(0.001,0,0,rep(-Inf,15)),upper = c(Inf,Inf,2,rep(Inf,15)),control = list(trace = 10, maxit = 100, factr = 1e-10),hessian = TRUE,plosdata = plosdata)

#end_time <- Sys.time()
#end_time - start_time
############### End of optimisation

############### Modify the Hessian if it has non-positive diagonal elements
timD <- diag(abs(diag(xxopt$hessian)))
ei <- eigen(xxopt$hessian)
new_hess <- ei$vectors %*% timD %*% t(ei$vectors)

######### Hypothesis test - This is the pseudo Wald's W test

# Derivative of the nested hypothesis test
H1 = matrix(0,nrow = 3,ncol = length(xxopt$par)-3)
H3 = diag(length(xxopt$par)-3)
H0 <- rbind(H1,H3)
# Compute the inverse of the observed fisher information matrix
Iinv <- solve(new_hess)

# Households in data
#m = max(plosdata$HH)
m = length(xxopt$par) - 3

# Compute the hypothesis test for the PseudoMLE estimator
# h <- 0
# for (i in 1:length(xxopt$par)-3)
# {
#   h[i] = xxopt$par[i+3]
# }
# h = as.matrix(h)
h = t(H0) %*% xxopt$par

# Calculate the chisquare statistic
chiQ <- m * t(h) %*% solve(t(H0) %*% Iinv %*% H0) %*% h

# Compute the p-value
pval <- pchisq(chiQ,length(xxopt$par)-3,lower.tail = FALSE)

# For each co-variate, calculate the p-value

chiq = 0;
for (i in 1:(length(xxopt$par)-3))
{
  chiqScore = m * xxopt$par[i+3] * (solve(H0[,i] %*% Iinv %*% H0[,i])) * xxopt$par[i+3]
  
  chiq[i] = pchisq(chiqScore,1,lower.tail = FALSE)
}

# Calculate AIC
AIC = -2*(-xxopt$value) + 2*(length(xxopt$par)-3)

# Compute the confidence interval. If you are maximising the likelihood, then the covariance matrix of the estimates iS asymptotically the inverse of the negative of the Hessian
#fisher_info <- solve(new_hess)
prop_sigma <- sqrt(diag(new_hess))
upper <- xxopt$par + 1.96/prop_sigma
lower <- xxopt$par - 1.96/prop_sigma

cidata <- rbind(upper,xxopt$par,lower)
View(cidata)

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
