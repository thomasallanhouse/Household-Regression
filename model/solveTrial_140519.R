# Clear the workspace
rm(list = ls())

# Source the llCalc.R function - Contains the loaded dynamic library and the called C function
# Makes llHoodFxn available
source("llCalc.R",echo = T)

# Set up the data
n = matrix(c(1,1,1,1,1),nrow = 1)
B = matrix(c(0.7,0.4,0.6,0.4,0.3),nrow = 1)
scale = 8.7
beta = matrix(0.0327,nrow = 5,ncol = 5)
HHdata = as.matrix(rbinom(5,1,0.5))

# Generate Addy's/Ball's matrix and solve it - a wrapper for the C function
out <- llHoodFxn(n, B, scale, beta, HHdata)

# Extract the HH configurations
indexSub <- matrix(out$indexS,ncol = 5)

P = out$P
sum(P)

# check out
plot(c(1:length(P)),P,type = "l",xlab = "Type configuration",ylab = "Probability")
points(out$loc,out$lHood,pch = 19,cex = 1.5,col = "red")

# Clean up
if (.Platform$OS.type == "windows"){
  # Windows
  dyn.unload("likelihoodCalc.dll")

} else if (.Platform$OS.type == "unix"){
  # Unix
  dyn.unload("likelihoodCalc.so")

} else {
  stop("Dynamic library can only be unloaded if loaded in Windows or Unix OS")
}
