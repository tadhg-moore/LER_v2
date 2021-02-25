qual_fun <- function(O, P){
  
  # remove datetime column
  O <- as.matrix(O[, -1])
  P <- as.matrix(P[, -1])
  
  # rmse
  rmse <- sqrt(mean((O - P)^2, na.rm = TRUE))
  
  
  # nash sutcliff
  nse <- 1 - sum((O - P)^2, na.rm = TRUE)/sum((O - mean(O, na.rm=TRUE))^2, na.rm = TRUE)
  
  # pearson corelation coef
  r <- sum((O - mean(O, na.rm = TRUE))*(P - mean(P, na.rm = TRUE)),
           na.rm = TRUE)/sqrt(sum((O - mean(O, na.rm = TRUE))^2, na.rm = TRUE)*
                                sum((P - mean(P, na.rm = TRUE))^2, na.rm = TRUE))
  
  # relative error
  re <- mean((P - O)/O, na.rm = TRUE)
  
  # mean absolute error
  mae <- mean(abs(O - P), na.rm = TRUE)
  
  # normalised mean absolute error
  LL <- -2 * sum(dnorm(O, mean = P, sd = 0.75, log = TRUE), na.rm = TRUE)
  
  qual <- data.frame(rmse = rmse, nse = nse, r = r, re = re, mae = mae, LL = LL)
  
  return(qual)
}