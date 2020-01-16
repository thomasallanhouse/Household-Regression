# Script checks if the solution is achievalbe with the matrix computed

# Clear the workspace
rm(list = ls())

# Set up the data
n = matrix(c(1,1,1,1),nrow = 1)
B = matrix(c(0.7,0.6,0.5,0.3),nrow = 1)
scale = 8.7
beta = matrix(0.0327,nrow = 4,ncol = 4)

# Load relevant functions and packages
source("finalSizeProbGen.R")
#library("prac")

# Run the function
output = finalSizeProbGen(n,beta,B,scale)
coefMatrix = output$coefMatrix
indexSub = output$indexSub

ones <- matrix(1,nrow = length(coefMatrix[,1]),ncol = 1)
P = solve(coefMatrix, ones)
sum(P)

# check out
plot(P,type = "l")
