########################
##Helper Functions#######
#########################

generate_mother <- function(alpha,beta,length_of_seq) {
  mother_chromatin = vector(length = length_of_seq)
  
mother_chromatin[1] = 1
for (i in 2:length_of_seq) {
  u = runif(1, 0, 1)
  if (mother_chromatin[i - 1] == 1) {
    mother_chromatin[i] = if (u <= alpha)
      1
    else
      0
  }
  else {
    mother_chromatin[i] = if (u <= (beta))
      0
    else
      1
  }
}
return(mother_chromatin)
}

calculate_0_length <- function(seq) {
  j = 1
  gap = vector();
  no_of_inter_zeros = 0
  for (i in 1: length(seq)) {
    if(seq[i] == 0) {
      no_of_inter_zeros = no_of_inter_zeros+1
      
    }
    else if(seq[i] == 1) {
      if(no_of_inter_zeros > 0) {
        gap[j] = no_of_inter_zeros
        no_of_inter_zeros = 0
        j = j + 1
      }
    }
  }
  return(mean(gap))
}

calculate_1_length <- function(seq) {
  j = 1
  gap = vector();
  no_of_inter_ones = 0
  for (i in 1: length(seq)) {
    if(seq[i] == 1) {
      no_of_inter_ones = no_of_inter_ones+1
      
    }
    else if(seq[i] == 0) {
      if(no_of_inter_ones > 0) {
        gap[j] = no_of_inter_ones
        no_of_inter_ones = 0
        j = j + 1
      }
    }
  }
  return(mean(gap))
}

check_if_region_b <- function(a,b) {
  if((a < 2*b) && (((a^2)/4) > ((1-a)*(1-b)/2))) {
    return(TRUE)
  }
    else {
      return(FALSE)
    }
  }
