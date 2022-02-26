compute_acc_path_metric <- function(possible_1_paths,no_of_interim_zeros,alpha, beta) 
{
  #yn = 0,xn = 0, xn-1 = 1
  branch_metric_1_0 = (1-alpha)
  #yn = 0,xn = 1, xn-1 = 1
  branch_metric_1_1 = 0.5 * alpha
  #yn = 0,xn = 0, xn-1 = 0
  branch_metric_0_0 = beta
  #yn = 0,xn = 1, xn-1 = 0
  branch_metric_0_1 = 0.5 * (1-beta)
  
  if(length(possible_1_paths) == 1) {
    no_of_1_paths = 1
  }
  else {
    no_of_1_paths = dim(possible_1_paths)[1]
  }
  acc_path_metric = 0
  #print(paste("no of paths"), no_of_1_paths)
  #print(paste("possible_1_paths"), possible_1_paths)
  
  
  for(v in 1: no_of_1_paths) {
    path_metric = 0
    for(w in 1: no_of_interim_zeros) {
      #1st bit calculation
      if (w == 1 ) {
        if ( possible_1_paths[w] == 1) {
             path_metric = branch_metric_1_1
        }
        else if (possible_1_paths[w] == 0) {
            path_metric = branch_metric_1_0
        }
      }
      # Last bit 
      else if (w == no_of_interim_zeros) {
        if ( possible_1_paths[w] == 1) {
          path_metric = path_metric * branch_metric_1_1
        }
        else if (possible_1_paths[w] == 0) {
          path_metric = path_metric * branch_metric_0_1
        }
      }
      # Remaining bits between first and last
      
      else  {
        if((possible_1_paths[w-1] == 1) && (possible_1_paths[w] == 1)) {
          path_metric = path_metric * branch_metric_1_1
        }
        else if((possible_1_paths[w-1] == 1) && (possible_1_paths[w] == 0)) {
          path_metric = path_metric * branch_metric_1_0 
          
        }
        else if((possible_1_paths[w-1] == 0) && (possible_1_paths[w] == 1)) {
          path_metric = path_metric * branch_metric_0_1 
  
        }
        else if((possible_1_paths[w-1] == 0) && (possible_1_paths[w] == 0)) {
          path_metric = path_metric * branch_metric_0_0 
  
        }
      }
    }
    acc_path_metric = path_metric + acc_path_metric
    path_metric = 0
    
  
  }
return(acc_path_metric)  
  
  
}

bcjr_decode <- function(daughter_chromatin,alpha,beta) {
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
       # print("INtermediate zeros")
       # print(no_of_interim_zeros)
        poss_combinations = expand.grid(replicate(no_of_interim_zeros, 0:1, simplify = FALSE))
      #  print(poss_combinations)
        for (k in 1: no_of_interim_zeros) {
          # kth bit
            possible_1_paths = poss_combinations[poss_combinations[,k] == 1,]
            acc_path_metric_1 = compute_acc_path_metric(possible_1_paths,no_of_interim_zeros,alpha,beta)
              
            possible_0_paths = poss_combinations[poss_combinations[,k] == 0,]
            acc_path_metric_0 = compute_acc_path_metric(possible_0_paths,no_of_interim_zeros,alpha,beta)
            
            # boundary conditions in 0 island
            ec_daughter_chromatin[i - no_of_interim_zeros -1] = 1
            ec_daughter_chromatin[i] = 1
            
            if(acc_path_metric_0 > acc_path_metric_1) {
              ec_daughter_chromatin[i- no_of_interim_zeros+ k -1] = 0
            }
            else {
              ec_daughter_chromatin[i- no_of_interim_zeros+ k -1] = 1
              
            }
          # print(paste(k,"th bit decision: "),ec_daughter_chromatin[i- no_of_interim_zeros+ k -1])
            
        }
        no_of_interim_zeros = 0
        }
       else {
        # consecutuve 1's in the y sequence
        ec_daughter_chromatin[i] = 1
       }
    }
  }
  no_of_interim_zeros = 0
  
  return(ec_daughter_chromatin)
  }


