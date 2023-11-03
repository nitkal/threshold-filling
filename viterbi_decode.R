viterbi_decode <- function(daughter_chromatin,alpha,beta) {
  p=0.5
  branch_metric = vector(length=4)
  
  
  #yn = 0,xn = 0, xn-1 = 1
  branch_metric[1] = (1-alpha)
  #yn = 0,xn = 1, xn-1 = 1
  branch_metric[2] = p * alpha
  #yn = 0,xn = 0, xn-1 = 0
  branch_metric[3] = beta
  #yn = 0,xn = 1, xn-1 = 0
  branch_metric[4] = p * (1-beta)
  
  #for (k in 1:no_of_daughters) {
  y = daughter_chromatin
  length_of_seq = length(daughter_chromatin)
  no_of_interim_zeros = 0
  ec_daughter_chromatin = vector(length = length_of_seq)
  ec_daughter_chromatin[1] = 1
  
  #print(y)
  
  for ( i in 1:length(y)) {
    #print(y[i])
    if(y[i] == 0) {
      no_of_interim_zeros  = no_of_interim_zeros + 1
      
    }
    
    else {
       #print("its a 1")
      
      if (no_of_interim_zeros > 0) {
       # print("intermediate zeroes present, algo starts")
        node_path_metric_xn_0 = array(dim = c(no_of_interim_zeros + 1, 2))
        
        node_path_metric_xn_1 = array(dim = c(no_of_interim_zeros + 1, 2))

        #For each 0, find the node metrics
        for (q in 1:(no_of_interim_zeros + 1)) {
          #print("q")
          #print(q)
          if(q == 1) {
            
            # Previously decoded starting border as 1
            node_path_metric_xn_0[q,1] = branch_metric[1]
            node_path_metric_xn_1[q,1] = branch_metric[2]
            
            node_path_metric_xn_0[q,2] = 1
            node_path_metric_xn_1[q,2] = 1
            
            
          }
          
          else {
            # Current node is xn = 0. Choose the best of the incoming path
            path_metric_1_0 = node_path_metric_xn_1[q-1,1]*branch_metric[1]
            path_metric_0_0 = node_path_metric_xn_0[q-1,1]*branch_metric[3]
            
            if(path_metric_1_0 >= path_metric_0_0) {
              
              node_path_metric_xn_0[q,1] = path_metric_1_0
              node_path_metric_xn_0[q,2] = 1
              
            }
            else {
              node_path_metric_xn_0[q,1] = path_metric_0_0
              node_path_metric_xn_0[q,2] = 0
            }
            
            # Current node is xn = 1. Choose best of incoming path
            path_metric_0_1 = node_path_metric_xn_0[q-1,1]*branch_metric[4]
            path_metric_1_1 = node_path_metric_xn_1[q-1,1]*branch_metric[2]
            
            if(path_metric_1_1 >= path_metric_0_1) {
              
              node_path_metric_xn_1[q,1] = path_metric_1_1
              node_path_metric_xn_1[q,2] = 1
              
            }
            else {
              node_path_metric_xn_1[q,1] = path_metric_0_1
              node_path_metric_xn_1[q,2] = 0
            }
            
          }
          
        }
      #  print("node_path_metric_xn_0")
       # print(node_path_metric_xn_0)
       # print("node_path_metric_xn_1")
       # print(node_path_metric_xn_1)
        
       # print(no_of_interim_zeros)
        ec_daughter_chromatin[i] = 1 # Border 1
        for ( q in (no_of_interim_zeros+1) :2) {
      #    print("q")
          #print(q)
         
          r = no_of_interim_zeros + 2 - q 
          #print("r")
          #print(r)
          

          if(q == (no_of_interim_zeros+1)) {
            #Border condition, yn = 1
            ec_daughter_chromatin[i-r] = node_path_metric_xn_1[q,2]
          }
          else {
            #If the next bit has a path from 0, choose 0 path, else 1 path
            if (ec_daughter_chromatin[i-r + 1] == 0) {
              ec_daughter_chromatin[i-r] = node_path_metric_xn_0[q,2]
            }
            else if(ec_daughter_chromatin[i-r + 1] == 1){
              ec_daughter_chromatin[i-r] = node_path_metric_xn_1[q,2]
            }
            
          }
          #print("ec_daughter_chromatin")
          #print(ec_daughter_chromatin)
          
        }
        no_of_interim_zeros = 0
      }
     
      
      else {
        # consecutuve 1's in the y sequence
        ec_daughter_chromatin[i] = 1
        
      }
      #print("final ec_daughter_chromatin")
      #print(ec_daughter_chromatin)
      
    
    
  }
  }
  corrected_chromatin = ec_daughter_chromatin
}

viterbi_decode_bias <- function(daughter_chromatin,alpha,beta,p,mu) {
  
  branch_metric = vector(length=4)
  
  
  #yn = 0,xn = 0, xn-1 = 1
  branch_metric[1] = (1-alpha)*(1-((1-p)*mu))
  #yn = 0,xn = 1, xn-1 = 1
  branch_metric[2] = p * alpha
  #yn = 0,xn = 0, xn-1 = 0
  branch_metric[3] = beta*(1-((1-p)*mu))
  #yn = 0,xn = 1, xn-1 = 0
  branch_metric[4] = p * (1-beta)
  
  #for (k in 1:no_of_daughters) {
  y = daughter_chromatin
  length_of_seq = length(daughter_chromatin)
  no_of_interim_zeros = 0
  ec_daughter_chromatin = vector(length = length_of_seq)
  ec_daughter_chromatin[1] = 1
  
  #print(y)
  
  for ( i in 1:length(y)) {
    #print(y[i])
    if(y[i] == 0) {
      no_of_interim_zeros  = no_of_interim_zeros + 1
      
    }
    
    else {
      #print("its a 1")
      
      if (no_of_interim_zeros > 0) {
        # print("intermediate zeroes present, algo starts")
        node_path_metric_xn_0 = array(dim = c(no_of_interim_zeros + 1, 2))
        
        node_path_metric_xn_1 = array(dim = c(no_of_interim_zeros + 1, 2))
        
        #For each 0, find the node metrics
        for (q in 1:(no_of_interim_zeros + 1)) {
          #print("q")
          #print(q)
          if(q == 1) {
            
            # Previously decoded starting border as 1
            node_path_metric_xn_0[q,1] = branch_metric[1]
            node_path_metric_xn_1[q,1] = branch_metric[2]
            
            node_path_metric_xn_0[q,2] = 1
            node_path_metric_xn_1[q,2] = 1
            
            
          }
          
          else {
            # Current node is xn = 0. Choose the best of the incoming path
            path_metric_1_0 = node_path_metric_xn_1[q-1,1]*branch_metric[1]
            path_metric_0_0 = node_path_metric_xn_0[q-1,1]*branch_metric[3]
            
            if(path_metric_1_0 >= path_metric_0_0) {
              
              node_path_metric_xn_0[q,1] = path_metric_1_0
              node_path_metric_xn_0[q,2] = 1
              
            }
            else {
              node_path_metric_xn_0[q,1] = path_metric_0_0
              node_path_metric_xn_0[q,2] = 0
            }
            
            # Current node is xn = 1. Choose best of incoming path
            path_metric_0_1 = node_path_metric_xn_0[q-1,1]*branch_metric[4]
            path_metric_1_1 = node_path_metric_xn_1[q-1,1]*branch_metric[2]
            
            if(path_metric_1_1 >= path_metric_0_1) {
              
              node_path_metric_xn_1[q,1] = path_metric_1_1
              node_path_metric_xn_1[q,2] = 1
              
            }
            else {
              node_path_metric_xn_1[q,1] = path_metric_0_1
              node_path_metric_xn_1[q,2] = 0
            }
            
          }
          
        }
        #  print("node_path_metric_xn_0")
        # print(node_path_metric_xn_0)
        # print("node_path_metric_xn_1")
        # print(node_path_metric_xn_1)
        
        # print(no_of_interim_zeros)
        ec_daughter_chromatin[i] = 1 # Border 1
        for ( q in (no_of_interim_zeros+1) :2) {
          #    print("q")
          #print(q)
          
          r = no_of_interim_zeros + 2 - q 
          #print("r")
          #print(r)
          
          
          if(q == (no_of_interim_zeros+1)) {
            #Border condition, yn = 1
            ec_daughter_chromatin[i-r] = node_path_metric_xn_1[q,2]
          }
          else {
            #If the next bit has a path from 0, choose 0 path, else 1 path
            if (ec_daughter_chromatin[i-r + 1] == 0) {
              ec_daughter_chromatin[i-r] = node_path_metric_xn_0[q,2]
            }
            else if(ec_daughter_chromatin[i-r + 1] == 1){
              ec_daughter_chromatin[i-r] = node_path_metric_xn_1[q,2]
            }
            
          }
          #print("ec_daughter_chromatin")
          #print(ec_daughter_chromatin)
          
        }
        no_of_interim_zeros = 0
      }
      
      
      else {
        # consecutuve 1's in the y sequence
        ec_daughter_chromatin[i] = 1
        
      }
      #print("final ec_daughter_chromatin")
      #print(ec_daughter_chromatin)
      
      
      
    }
  }
  corrected_chromatin = ec_daughter_chromatin
}

