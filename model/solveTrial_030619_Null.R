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
load("optLinearNull.RData")

# Source the llCalc.R function - Contains the loaded dynamic library and the called C function
# Makes llHoodFxn available
source("llCalc.R",echo = TRUE)

source("regressionFunction_LinearNull.R",echo = TRUE)

# Load the data as tibble
plosdata <- readr::read_csv("plosData_300519.csv")

# Run the model - This is 10 times faster than the origin Matlab - Largely because the matrix computations and sol are done in C
out <- regressionFunction_LinearNull(xxopt$par,plosdata)

# Start the timing
#start_time <- Sys.time()

# Run optimizer
#xxopt <- optim(par = x0,fn = regressionFunction_LinearB,method = "L-BFGS-B",lower = c(0.01,0,0,rep(-Inf,15)),upper = c(1,Inf,2,rep(Inf,15)),control = list(trace = 10, maxit = 1000, factr = 1e-15),hessian = TRUE,plosdata = plosdata)

# Stop timing and calculate difference
#end_time <- Sys.time()
#end_time - start_time


# Compute the inverse of the observed fisher information matrix
Iinv <- solve(xxopt$hessian)

# Calculate AIC
AIC = -2*(-xxopt$value) + 2*(length(xxopt$par))

# Compute the confidence interval. If you are maximising the likelihood, then the covariance matrix of the estimates iS asymptotically the inverse of the negative of the Hessian
#fisher_info <- solve(new_hess)
prop_sigma <- sqrt(diag(new_hess))
upper <- xxopt$par + 1.96/prop_sigma
lower <- xxopt$par - 1.96/prop_sigma

# Compute CIs
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
