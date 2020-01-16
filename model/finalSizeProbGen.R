# Code generates the relevant matrix for the system of equations given by
# Addy 1991
# Input
#   beta: Transmission parameter - As coded currently, it can only be a
#   vector or a scalar. TO DO: Code for matrix
#   B: Probability that a susceptible of type i will escape infection from
#   outside of the population during the entire course of the epidemic.
#   B is a column vector
#   Scale: Parameters for the Gamma distribution. Can be coded to
#   be any arbitrary distribution
#
# Written by Tim Kinyanjui
# University of Manchester
# on 9th May 2019 - week
#
# Might become a package later or a ShinyApp

# Define function  and its inputs
finalSizeProbGen <- function(n,beta,B,scale){

  # Length of N
  lenN = length(n)

  # Each type has one single index case
  m = matrix(1,ncol = lenN)

  # Laplace transform of Gamma distribution.
  phi = function(x) 1/((1+scale*x)**(1/scale))

  # All the occurences of possible combination
  allComb = prod(n+1)

  # Initialise matrix that stores the co-efficients
  coefMatrix = matrix(0,nrow = allComb,ncol = allComb)

  # All possible combinations are stored in indexSub. Corresponds to the sub-scripts of P i.e. P_indexSub
  indexSub = matrix(0,nrow = allComb, ncol = lenN)

  # CumP
  cumP = c(1,cumprod(n[c(1:length(n)-1)] + 1))
  J = matrix(0,nrow = 1, ncol = length(n))

  # Construct loop
  for (indexAll in 1:allComb) {

    # Determine the storage location on the x location
    for (i in lenN:2){

      # Loop
      if (J[i] > n[i]){
        J[i] = 0
        J[i-1] = J[i-1] + 1
      }
      else {
        # The vector should be fine as it is
      }
    }

    ii = sum(cumP * J) + 1

    # Store the configurations
    indexSub[ii,c(1:lenN)] = J

    # Determine the storage locations
    omega = matrix(0,nrow = 1,ncol = length(n))

    for (w in 1:prod(J+1)) {

      for (i in lenN:2){

        if (omega[i] > J[i]){

          omega[i] = 0
          omega[i-1] = omega[i-1] + 1
        }
        else {
          # Should remain the same
        }
      }
      jj = sum(cumP * omega) + 1

      # Calculate the co-efficients here. Phi is the Laplace transform of arbitrary infectious perios distribution
      prodDD = 1

      for (cont in 1:length(omega)) {

        prodDD[cont] = phi(sum((n - J)*beta[,cont]))**(omega[cont] + m[cont])
      }

      prodD = prod(prodDD)

      # Matrix from Frank Ball 1986
      # coefMatrix[ii,jj] = prod(choose(n - omega,J - omega))/(prodD * prod(choose(n,J)) * prod(B**(n-J)))

      # Matrix from Addy 1991
      coefMatrix[ii,jj] = prod(choose(J,omega))/(prodD * prod(B**(n-J)))

      # Increase omega counter
      omega[lenN] = omega[lenN] + 1
    }

    # Increment J
    J[lenN] = J[lenN] + 1
  }

  # End of function and return value
  output = list(coefMatrix = coefMatrix, indexSub = indexSub)
  return(output)
}
