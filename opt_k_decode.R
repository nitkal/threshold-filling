###################################
###Function for Threshold-K filling
##################################


opt_k_decode <- function(y, thresh_k) {

  no_of_interim_zeros = 0
  length_of_seq = length(y)
  ec_daughter_chromatin = vector(length = length_of_seq)
  ec_daughter_chromatin[1] = 1
  
  for (i in 1:length(y)) {
    if (y[i] == 0) {
      no_of_interim_zeros  = no_of_interim_zeros + 1
    }
    
    else {

      if (no_of_interim_zeros > 0) {
  
        if (thresh_k >= no_of_interim_zeros) {
            for (g in 1:no_of_interim_zeros) {
            ec_daughter_chromatin[i - g] = 1
          }
        }
        else {
          #Leave all intermediate zeros as 0's
          for (g in 1:no_of_interim_zeros) {
            ec_daughter_chromatin[i - g] = 0
          }
        }
        ec_daughter_chromatin[i] = 1
        no_of_interim_zeros = 0
      }
      else {
        # consecutuve 1's in the y sequence
        ec_daughter_chromatin[i] = 1
        
      }
      
      
    }
  }
  return(ec_daughter_chromatin)
}
