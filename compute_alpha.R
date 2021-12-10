##################################################
#####Helper functions to compute alpha, beta for a sequence
#################################################

compute_alpha <- function(mother_chromatin) {
  
  length_of_seq = length(mother_chromatin)
  n11 = 0
  n10=0
  n01 = 0
  n00 = 0
  for(i in 2:length_of_seq) {
    if ((mother_chromatin[i-1] == 1 ) && (mother_chromatin[i]==1) ) {
      n11 = n11 + 1
    }
    else if( (mother_chromatin[i-1] == 1 ) && (mother_chromatin[i]==0)) {
      n10 = n10 + 1
    }
    else if( (mother_chromatin[i-1] == 0 ) && (mother_chromatin[i]==1)) {
      n01 = n01 + 1
    }
    else if( (mother_chromatin[i-1] == 0 ) && (mother_chromatin[i]==0)) {
      n00 = n00 + 1
    }
    
  }
  ones = 0
  for( i in 1: length_of_seq) {
    if(mother_chromatin[i] == 1) {
      ones = ones + 1
    }
  }
  
  alpha = n11/(n11 + n10 )
  beta = n00/(n01 + n00)
  return( list("a"=alpha, "b"= beta))
}

compute_interim_zeros <- function(chromatin) {
  length_of_seq = length(chromatin)
  inter_zeros_vector = vector()
  no_of_inter_zeros =0
  b = 1
  for (i in 2:length_of_seq) {
    
    if (chromatin[i] == 0) {
      no_of_inter_zeros = no_of_inter_zeros + 1
    }
    else {
      if(no_of_inter_zeros != 0) {
        inter_zeros_vector[b] = no_of_inter_zeros
        b = b + 1
        no_of_inter_zeros = 0
      }
      else {
        next()
      }
    }
  }
  return(inter_zeros_vector)
}

compute_interim_ones <- function(chromatin) {
  length_of_seq = length(chromatin)
  inter_ones_vector = vector()
  no_of_inter_ones =0
  b = 1
  for (i in 2:length_of_seq) {
    
    if (chromatin[i] == 1) {
      no_of_inter_ones = no_of_inter_ones + 1
    }
    else {
      if(no_of_inter_ones != 0) {
        inter_ones_vector[b] = no_of_inter_ones
        b = b + 1
        no_of_inter_ones = 0
      }
      else {
        next()
      }
    }
  }
  return(inter_ones_vector)
}

seq = c(1,0,1,1,1,0,0,1,1,1,0,0,0,1,1)
ab = compute_alpha(seq)