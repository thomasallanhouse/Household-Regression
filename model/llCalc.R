# Compile C code- This needs to be done only once for any operating system - Checked for windows and unix-like OSes
# Compile using the following as we are now using lapack library for matrix solution
# R CMD SHLIB -lRlapack likelihoodCalc.c - in windows
# R CMD SHLIB -llapack likelihoodCalc.c - in Linux (Done in Ubuntu)
#
# To see the options flags for compiling do: R CMD config LAPACK_LIBS

# Loads the shared library containing the llHoodFxn
#if (.Platform$OS.type == "windows"){
  # Windows
#  dyn.load("likelihoodCalc.dll")

#} else if (.Platform$OS.type == "unix"){
  # Unix
#  dyn.load("likelihoodCalc.so")

#} else {
#  stop("Dynamic library can only be loaded in Windows or Unix OS")
#}

# A wrapper to simplify calling the C function from R
llHoodFxn <- function(n, B, scale, beta, HHdata){

  # Set up the relevant data
  # Length of N
  lenN = length(n)

  # Each type has one single index case
  m = matrix(1,ncol = lenN)

  # All the occurences of possible combination
  allComb = prod(n+1)

  J = matrix(0,nrow = 1, ncol = length(n))

  omega = J

  cumP = c(1,cumprod(n[c(1:length(n)-1)] + 1))

  # Initialise matrix that stores the co-efficients
  coefMatrix = matrix(0,nrow = allComb,ncol = allComb)

  # All possible combinations are stored in indexSub. Corresponds to the sub-scripts of P i.e. P_indexSub
  indexSub = matrix(0,nrow = allComb, ncol = lenN)
  
  # Initialise P the probabilities
  P = matrix(0,nrow = allComb)
  P[1] = 1
  
  # Likelihood is returned here
  lHood = 0;

  # Call the C function
  out1 <- .C("llHoodFxn", as.integer(n), as.double(beta), as.double(B), loc = as.double(scale), as.integer(lenN), as.integer(m), as.integer(allComb), as.integer(J), as.integer(omega), as.integer(cumP), coefMat = as.double(coefMatrix), indexS = as.integer(indexSub), P = as.double(P), as.integer(HHdata), lHood = as.double(lHood))
}
