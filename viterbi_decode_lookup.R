source('viterbi_decode_0_islands.R')
lookup_pattern_from_map <- function(no_of_inter_zeros,alpha,beta,pattern.df) {
  # if(!exists("pattern.df")) {
  #   pattern.df <- data.frame(numeric(),numeric(),numeric(),character(),stringsAsFactors = FALSE)
  #   colnames(pattern.df) <- c("no_0","alpha","beta","pattern_int")
  # }
  
  opt_pattern = pattern.df[which( pattern.df$no_0==no_of_inter_zeros & pattern.df$alpha==alpha & pattern.df$beta==beta) , "pattern_int"]
  if(length(opt_pattern >0)) {
    return(opt_pattern)
  }
  else {
    return(NA)
  }
}


viterbi_decode_lookup <- function(daughter_chromatin,alpha,beta,viterbi_pattern.df) {
  if(is.null(viterbi_pattern.df)) {
    viterbi_pattern.df <- data.frame(numeric(),numeric(),numeric(),character(),stringsAsFactors = FALSE)
    colnames(viterbi_pattern.df) <- c("no_0","alpha","beta","pattern_int")
  }
  y = daughter_chromatin
  length_of_seq = length(daughter_chromatin)
  no_of_interim_zeros = 0
  ec_daughter_chromatin = vector(length = length_of_seq)
  ec_daughter_chromatin[1] = 1
  
  for ( i in 1:length(y)) {
    #print(y[i])
    if(y[i] == 0) {
      no_of_interim_zeros  = no_of_interim_zeros + 1
    }
    
    else {
      #print("its a 1")
      
      if (no_of_interim_zeros > 0) {
       # print(paste("Looking up best pattern for ", no_of_interim_zeros, " zeros, ",alpha," alpha ,",beta," beta"))
        
        best_pattern = lookup_pattern_from_map(no_of_interim_zeros,alpha,beta,viterbi_pattern.df)
        if(!is.na(best_pattern)) {
          best_pattern = unlist(strsplit(best_pattern,""))
          
         # print(best_pattern)
        }
        else  {
        #  print(paste("computing best pattern for ", no_of_interim_zeros, " zeros, ",alpha," alpha ,",beta," beta"))
          best_pattern = viterbi_decode_0_islands(no_of_interim_zeros,alpha,beta)
          viterbi_pattern.df[nrow(viterbi_pattern.df)+1,]<- c(no_of_interim_zeros,alpha,beta,(gsub(", ","",toString(best_pattern))))
          
        #  print(viterbi_pattern.df)
          
        }
        
        ec_daughter_chromatin[i] = 1
        
        for (p in 1:no_of_interim_zeros) {
          ec_daughter_chromatin[i -no_of_interim_zeros +p - 1] = best_pattern[p]
        }
        no_of_interim_zeros = 0
      }
      else {
        # consecutuve 1's in the y sequence
        ec_daughter_chromatin[i] = 1
      }
    }
  }
  
  ec_daughter_chromatin = tryCatch({as.integer(unlist(strsplit(ec_daughter_chromatin, "")))}, error=function(e){print(ec_daughter_chromatin)}) 
  #print(paste("Computed output for alpha= ",alpha," beta=",beta))
  return(list("y"=ec_daughter_chromatin,"viterbi.pattern.df"=viterbi_pattern.df))
  
}

