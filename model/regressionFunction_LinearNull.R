# Calculates the likelihood of the model - Within household transmission
#
# Define function  and its inputs
regressionFunction_LinearNull <- function(xopt,plosdata){
  
  # save(list = c("plosdata","xopt"),file = "optLinear.RData")
  # Load saved parameter
  # load("optLinear.RData")
  
  # Extract the required parameters
  print(xopt)
  probB = xopt[1]
  scale = xopt[2]
  alpha = xopt[3]
  b = xopt[4]
  
  # Extract a vector of unique households (IDs) - All members of the same household have the same ID
  uniqueHH <- unique(plosdata$HH)
  lHood = rep(0,length(uniqueHH))
  
  # Loop through to process all households
  for (ii in 1:length(uniqueHH)){
    
    # Extract a subset and drop HH (use minus sign)
    plosdata %>% dplyr::filter(HH == uniqueHH[ii]) %>% select(-HH) -> subHH
    
    # Extract the finalSize
    finalSize <- subHH$secondary_case_with_ILI
    
    # Drop the dependent variable
    subHH %>% select(-secondary_case_with_ILI) -> subHH
    
    # Linear regression to estimate tranmsmission
    lambda = b * rep(1,nrow(subHH))
    
    B = probB * rep(1,nrow(subHH))
    
    # Household size
    n = rep(1,length(finalSize))
    
    # Set up the data
    if (length(finalSize) == 1){
      
      # For trivial household
      beta = lambda[1] * matrix(1,nrow = length(finalSize),ncol = length(finalSize))
      
    } else {
      
      # For Household with more than 1 member
      beta = t((lambda / ((length(finalSize) - 1) ** alpha)) * matrix(1, nrow = length(finalSize), ncol = length(finalSize)))
      
    }
    
    # Generate Addy's/Ball's matrix and solve it and calculate the likelihood - a wrapper for the C function
    out <- llHoodFxn(n, B, scale, beta, finalSize)
    
    # Save the likelihood
    lHood[ii] <- out$lHood
    
    
    ############################# Testing part #############################
    # If you run into c stack overflow error increase the size from terminal using ulimit -s 65535 to make 64MB (Cstack_info())
    #P = out$P
    #sum(P)
    #plot(c(1:length(P)),P,type = "l",xlab = "Type configuration",ylab = "Probability")
    #points(out$loc,out$lHood,pch = 19,cex = 1.5,col = "red")
    ################################ Remove ################################
  }
  
  # For debugging optimisation
  print(-1*sum(log(lHood)))
  
  # End of function
  return(-1*sum(log(lHood)))
}
