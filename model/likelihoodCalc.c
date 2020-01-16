/* Function generated the relevant matrix for the system of equations given by Addy
 Input
 beta: Transmission parameter - As coded currently, it can only be a
 vector or a scalar. TO DO: Code for matrix
 B: Probability that a susceptible of type i will escape infection from
 outside of the population during the entire course of the epidemic.
 B is a column vector
 Scale and Shape: Parameters for the Gamma distribution. Can be coded to
 be any arbitrary distribution
 
 Written by Tim Kinyanjui
 University of Manchester
 Started on 14th May 2019
 
 */

#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <R_ext/Lapack.h> // Include to solve linear system of equations

// Declare own functions - consider shoving them inside a header file instead

//  This returns a 1 if two vectors (arrays) are same i.e. all elements are the same
int areSame (int *tempHH, int *HHdata, int *lenN);

// This returns the product of nchoosek's
int nchoosekP (int *N, int *K, int *lenN);

// This returns n choose k
int nchoosek (int n, int k);

void llHoodFxn(int *n, double *beta, double *B, double *scale, int *lenN, int *m, int *allComb, int *J, int *omega, int *cumP, double *coefMatrix, int *indexSub, double *P, int *HHdata, double *lHood){
  
  // Declaration of variables here
  int indexAll,i1,tt,ii = 0,jj = 0,w,prodJ = 1,i2,tt1,tt2,bi = 0,tt3,prodChoose = 0,ipiv[*allComb],info,nrhs = 1,tempHH[*lenN],sameOut = 99;
  double prodnJ[*lenN], bbeta[*lenN], sumnJ = 0, prodD = 1, prodBnj = 1, Q[*allComb][*allComb], BB[*allComb];
  
  // Initialise Q - output matrix
  for (tt = 0; tt <= *allComb -1; tt++){
    for (tt1 = 0; tt1 <= *allComb - 1; tt1++){
      Q[tt][tt1] = 0;
    }
  }
  
  for (indexAll = 1; indexAll <= *allComb; indexAll++){
    
    // Determine the storage location on the x location
    for (i1 = *lenN; i1 >= 2; i1--){
      
      if (J[i1-1] > n[i1-1]){
        
        J[i1-1] = 0;
        
        J[i1-2] = J[i1-2] + 1;
      }
      
    }
    
    // For loop to calculate the storage location
    for (tt = 0; tt <= *lenN-1; tt++){
      ii = ii + (cumP[tt] * J[tt]);
    }
    
    // Store the configurations - assumes the matrix from R is expanded column-wise i.e. stacked columns
    for (tt = 0; tt <= *lenN - 1; tt++){
      indexSub[ii + (tt * (*allComb))] = J[tt];
    }
    
    // Determine the stopping rule for the y-location loop below
    for (tt = 0; tt <= *lenN-1; tt++){
      prodJ = prodJ * (J[tt] + 1);
      omega[tt] = 0;
    }
    
    // Determine storage location in orthorgonal direction i.e. y direction
    for (w = 1; w <= prodJ; w++){
      
      for (i2 = *lenN; i2 >= 2; i2--){
        
        if (omega[i2-1] > J[i2-1]){
          
          omega[i2-1] = 0;
          
          omega[i2-2] = omega[i2-2] + 1;
        }
      }
      
      // For loop to calculate the storage location
      for (tt = 0; tt <= *lenN-1; tt++){
        jj = jj + (cumP[tt] * omega[tt]);
      }
      
      // Calculate the co-efficients here
      for (tt = 0; tt <= *lenN-1; tt++){
        
        // Extract the beta column of interest
        bi = 0;
        for (tt3 = tt*(*lenN); tt3 <= ((tt+1)*(*lenN))-1; tt3++){
          bbeta[bi] = beta[tt3];
          bi++;
        }
        
        for (tt1 = 0; tt1 <= *lenN-1; tt1++){
          sumnJ = sumnJ + ((n[tt1] - J[tt1]) * bbeta[tt1]);
        }
        
        // We'ne hard coded the infectious time distribution. An alternative is to allow a user defined dsistribution
        // Calculate the laplace transform of sumnJ
        prodnJ[tt] = (pow((1 / pow((1 + ((*scale) * sumnJ)),(1/(*scale)))),(omega[tt] + m[tt])));
        sumnJ = 0;
      }
      
      for (tt = 0; tt<= *lenN - 1; tt++){
        prodD = prodD * prodnJ[tt];
      }
      
      // Product for nchoosek
      prodChoose = nchoosekP(J,omega,lenN);
      
      for (tt = 0; tt <= *lenN - 1; tt++){
        prodBnj = prodBnj * (pow(B[tt],(n[tt] - J[tt])));
      }
      
      // Store the result in a two-D array - This is Addy's matrix
      Q[ii][jj] = prodChoose / (prodD * prodBnj);
      
      // Increment
      omega[*lenN-1] = omega[*lenN-1] + 1;
      
      // Reset the variables to re-use here
      jj = 0; prodD = 1; prodBnj = 1;
      for (tt = 0; tt<= *lenN - 1; tt++){
        prodnJ[tt] = 1;
      }
    }
    
    // Increment
    J[*lenN-1] = J[*lenN-1] + 1;
    
    // Reset the variables to re-use
    ii = 0;
    prodJ = 1;
    for (tt = 0; tt<= *lenN - 1; tt++){
      omega[tt] = 0;
    }
  }
  
  // Put matrix elements in an array to return to R or use with the dgesv fortran routine
  tt2 = 0;
  for (tt = 0; tt <= *allComb - 1; tt++){
    
    // Also generate the right hand side of the linear system
    BB[tt] = 1;
    
    for (tt1 = 0; tt1 <= *allComb - 1; tt1++){
      
      // Extract column-wise
      coefMatrix[tt2] = Q[tt1][tt];
      tt2++;
    }
  }
  
  // Compute the solution to a real system of linear equations using FORTRAN routine
  // A * X = B, using the DGESV fortran routine
  // A = coefMatrix; B = BB (input) and the solution is returned in BB as well
  
  // Use the R macro provided to call the DGESV fortran routine
  F77_CALL(dgesv)(allComb, &nrhs, coefMatrix, allComb, ipiv, BB, allComb, &info);
  
  // This call works too but could be platform dependent
  // dgesv_(allComb, &nrhs, coefMatrix, allComb, ipiv, BB, allComb, &info);
  
  // Check if solution exists
  if (info != 0){
    
    // If solution failed, return the original P vector
    // This is important especially during optimisation when the proposed parameter values make the matrix badly scaled
    error("Solution failed - original vector returned");
    
  } else if (info == 0) {
    
    // Store solution in P
    for (tt = 0; tt <= *allComb - 1; tt++){
      P[tt] = BB[tt];
    }
    
  }
  
  // Compute the likelihood
  
  tt = 0; // Initialise outer counter
  while(tt <= *allComb - 1){
    
    tt1 = 0; // Initialise inner counter
    while(tt1 <= *lenN - 1){
      
      tempHH[tt1] = indexSub[tt + (tt1 * (*allComb))];
      
      // Increment the inner counter
      tt1++;
    }
    
    // Check if the two vectors are the same. If same output = 1 else output = 0
    sameOut = areSame(tempHH, HHdata, lenN);
    
    if (sameOut == 1){
      
      // Vectors are the same - Exit and return tt which is the location of the likelihood
      // printf("P = %f, Loc: %d\n",P[tt],tt);
      *lHood = P[tt];
      
      // Return the location of likelihood in scale
      *scale = tt + 1; // R indexes from location 1 so add 1
      
      // Set tt so that they exit the loop without necessarily going through all elements
      tt = *allComb + 1;
    } else {
      
      // Increment the outer counter
      tt++;
    }
  }
  if (sameOut != 1){
    error("Configuration not found. Should really never be here. Contact me");
  }
}

/*
 ********************************************************************
 * These are helper functions which are called within the llHoodFxn *
 ********************************************************************
 */

// A function that checks if two vectors are the same
int areSame (int *tempHH, int *HHdata, int *lenN){
  
  // Declare and initialise variables here
  int out = 2,i=0;
  
  for (i = 0; i <= *lenN - 1; i++){
    
    // printf("%d||%d-->",tempHH[i],HHdata[i]);
    if (tempHH[i] == HHdata[i]){
      
      out = 1;
    } else{
      
      out = 0;
      break;
    }
  }
  
  // Retrurn value
  return out;
}

// A function that calculates the product of n choose k (combinatorial)
int nchoosekP (int *N, int *K, int *lenN){
  
  // Declare variables
  int ii, prodChoose = 1;
  
  // Pass them
  for (ii = 0; ii <= *lenN - 1; ii++){
    prodChoose = prodChoose * nchoosek(N[ii],K[ii]);
  }
  
  // Return
  return prodChoose;
}

// A function that calculates n choose k (combinatorial)
int nchoosek (int n, int k){
  
  // Declare variable
  int i, prodN = 1, prodK = 1, prodNK = 1, nCombk = 0;
  
  // Factorial n
  if (n == 0){
    prodN = 1;
  } else{
    for (i = 1; i <= n; i++){
      prodN = prodN * i;
    }
  }
  
  
  // Factorial k
  if (k == 0){
    prodK = 1;
  } else{
    for (i = 1; i <= k; i++){
      prodK = prodK * i;
    }
  }
  
  // Factorial n-k
  if ((n-k) == 0){
    prodNK = 1;
  } else{
    for (i = 1; i <= (n-k); i++){
      prodNK = prodNK * i;
    }
  }
  
  // Compute the nCk
  nCombk = prodN / (prodNK * prodK);
  
  // Return
  return nCombk;
}
