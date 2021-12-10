########################
#Program to compute the opt k decoding for HeLA data
# HeLda data sent by Google Drive
########################


library("vioplot")
source("../compute_alpha.R")
source("../opt_k_decode/opt_k_decode.R")
smoothen_values <- function(doi) {
  vals = as.double(doi[,4])
  length_of_seq = dim(doi)[1]
  avg_curr_nuc = vector(length = length_of_seq)
  for ( i in 1:length_of_seq) {
    prev_nuc = 0
    next_nuc = 0
    if ((i - 1) >= 1) {
      prev_nuc = vals[i-1]
    }
    curr_nuc = vals[i]
    if ((i + 1)<= length_of_seq) {
      next_nuc = vals[i+1]
    }
    
    avg_curr_nuc[i] = (prev_nuc + curr_nuc + next_nuc)/3
  }
  
  return(avg_curr_nuc)
  
}

#H3K27me3
#h3_k27me3_bed <- read.table("/Users/nitrama/research/data/Hela_groth_data/GSM2988386_H3K27me3_ChIP-seq_Parental_RPM_rep1.bedgraph",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="" )
#H3K4me3
h3_k27me3_bed <- read.table("/Users/nitrama/research/data/Hela_groth_data/GSM3227882_H3K4me3_ChIP-seq_Parental_RPM_rep1.bedgraph",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="" )
#H3k36me3
#h3_k27me3_bed <- read.table("/Users/nitrama/research/data/Hela_groth_data/GSM3227892_H3K36me3_ChIP-seq_Parental_RPM_rep1.bedgraph",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="" )
#H3K79me3
#h3_k27me3_bed <- read.table("/Users/nitrama/research/data/Hela_groth_data/GSM3227896_H3K79me3_ChIP-seq_Parental_RPM_rep1.bedgraph",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="" )




###Control H3


h3_bed <- read.table("/Users/nitrama/research/data/Hela_groth_data/GSM2988395_H3_ChIP-seq_Parental_RPM_rep1.bedgraph",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="" )

chr1_h3_k27me3_bed =h3_k27me3_bed[which(h3_k27me3_bed[,1] == "chr19"),]
#data_of_interest_m = chr1_h3_k27me3_bed[which((as.integer(chr1_h3_k27me3_bed[,2]) > 56484323 )),]
data_of_interest_m = chr1_h3_k27me3_bed[which((as.integer(chr1_h3_k27me3_bed[,2]) > 7047900 )),]

doi_m = data_of_interest_m[which(as.integer(data_of_interest_m[,3]) < 8713217),]

#doi_m = doi_m[12800:12899,]
doi_m1 = doi_m[,4]/max(doi_m[,4]) 



#doi_m = chr1_bed

avg_curr_nuc_m = smoothen_values(doi_m)
chr1_h3_bed =h3_bed[which(h3_bed[,1] == "chr19"),]
#data_of_interest_c = chr1_h3_bed[which((as.integer(chr1_h3_bed[,2]) > 56484323 )),]
data_of_interest_c = chr1_h3_bed[which((as.integer(chr1_h3_bed[,2]) > 7047900 )),]

doi_c = data_of_interest_c[which(as.integer(data_of_interest_c[,3]) < 8713217),]

#doi_c = doi_c[12800:12899,]

doi_c1 = doi_c[,4]/max(h3_bed[,4]) 
#doi_c = chr1_h3_bed
avg_curr_nuc_c = smoothen_values(doi_c)

h3k27me3_doi = data.frame("chr"=doi_m[,1],"start"= as.double(doi_m[,2]),end=as.double(doi_m[,3]),"val"=avg_curr_nuc_m)
h3_doi = data.frame("chr"=doi_c[,1],"start"= as.double(doi_c[,2]),end=as.double(doi_c[,3]),"val"=avg_curr_nuc_c)

mod_vals_ratio = h3k27me3_doi[,4]/max(h3_doi[,4])
control_vals_ratio = h3_doi[,4]/max(h3_doi[,4])
mean_h3_ratio = mean(control_vals_ratio)

length_of_seq = dim(h3_doi)[1]


thresh_k = seq(from = 1,to = 9, length.out = 9)
no_of_mothers = 200
no_of_daughters_per_mother = 100
error_bc = array(numeric(), c(no_of_mothers,no_of_daughters_per_mother))
error_ac = array(numeric(), c(no_of_mothers,no_of_daughters_per_mother,length(thresh_k)))
mothers = array(numeric(), c(no_of_mothers,length_of_seq))
alphas = vector(length=no_of_mothers)
betas = vector(length=no_of_mothers)


daughters = array(numeric(), c(no_of_mothers, no_of_daughters_per_mother,length_of_seq))

ec_daughters = array(numeric(), c(no_of_mothers, no_of_daughters_per_mother,length(thresh_k),length_of_seq))

for (mo in 1:no_of_mothers ) {
  
  #Mother chromatin
  
  mother_seq = vector(length = length_of_seq)
  random_seq = runif(length_of_seq, 0,1)
  
  
  k = 1
  for (j in 1: length_of_seq) {
    
    if ( (h3_doi[j,4] !=0)) {
      if (mod_vals_ratio[j] > mean_h3_ratio) {
        mother_seq[k] = 1
      }
      else {
        if (random_seq[j] <= (mod_vals_ratio[j]/mean_h3_ratio)) {
          mother_seq[k] = 1
        }
        else {
          mother_seq[k] = 0
        }
      }
      k = k+1
    }
  }
  
  
  mothers[mo,] = mother_seq
  
  mother_ab= compute_alpha(mother_seq)
  #print("Alpha - ")
  #print(ab$alpha)
  #print("Beta - ")
  #print(ab$beta)
  alphas[mo] = mother_ab$a
  betas[mo] = mother_ab$b
  
  length_of_seq = length(mother_seq)
  
  
  # Daughter Chromatin
  
  for (da in 1:no_of_daughters_per_mother) {
    prob_of_copy = rbinom(length_of_seq,1,0.5)
    #prob_of_daughter_histone = rbinom(length_of_seq,1,0.5)
    
    daughter_chromatin = vector(length = length_of_seq)
    #daughter_chromatin[1] = mother_seq[1]
    for (i in 2:length_of_seq) {
      
      daughter_chromatin[i] = ifelse(((mother_seq[i] == 1) && (prob_of_copy[i] == 1)), mother_seq[i] ,0)
      
    }
    daughters[mo,da,] = daughter_chromatin
    error_bc[mo,da] = sum(xor(mother_seq,daughter_chromatin))/length_of_seq
    
    
    
    #Opt K decoding
    alpha = mother_ab$a
    beta = mother_ab$b
    for (kk in 1: length(thresh_k)) {
      opt_k = thresh_k[kk]
      
      y = daughter_chromatin
      ec_daughter_chromatin = vector(length = length_of_seq)
      
      ec_daughter_chromatin = opt_k_decode(y,opt_k)
      
      
      
      
      ec_daughters[mo,da,kk,] = ec_daughter_chromatin
      error_ac[mo,da,kk] = sum(xor(mother_seq,ec_daughter_chromatin))/length_of_seq
    }
    
  }
}
##Violin plot
x0 = error_bc
x1 = error_ac[,,1]
x2 = error_ac[,,2]
x3 = error_ac[,,3]
x4 = error_ac[,,4]
x5 = error_ac[,,5]
x6 = error_ac[,,6]
x7 = error_ac[,,7]
x8 = error_ac[,,8]
x9 = error_ac[,,9]

# plot(x=0:1,                  # Make x range (0,1)
#      y=0:1,                     # Make y range (0,1)
#      type='n',                  # Don't actually plot anything
#      xlim=c(0,10.2),           # Weird xlim settings for two groups
#      ylim=c(min(c(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9)-0.005),max(c(x1,x2,x3,x4,x5,x6,x7,x8,x9))), # Y-axis limits based on data
#      axes=FALSE, ann=FALSE
# )
x0 = as.vector(x0)
x1 = as.vector(x1)
x2 = as.vector(x2)
x3 = as.vector(x3)
x4 = as.vector(x4)
x5 = as.vector(x5)
x6 = as.vector(x6)
x7 = as.vector(x7)
x8 = as.vector(x8)
x9 = as.vector(x9)



par(cex.lab=1.8, cex.axis=1.2,  cex.sub=1.5,mgp=c(2.5,1,0))

vioplot(lty=1,wex=1.1,col="darkgreen",add=FALSE, x0,x1,x2,x3,x4,x5,x6,x7,x8,x9, xlab=expression(italic("k"[t])),ylab="Mean error",names=c('0','1','2','3', '4','5','6','7','8','9'))
 axis(side=2, las=2)
 axis(side=1, at=c(0,1,2,3,4,5,6,7,8,9), c('0','1','2','3', '4','5','6','7','8','9'))

## Mean plots

mean_mother = apply(mothers,c(2),mean)
mean_daughter = apply(daughters,c(3),mean)
mean_ec_daughter_2 = apply(ec_daughters,c(3,4),mean)
error_ac_k = apply(error_ac,c(3),mean)

#par(mfrow= c(4,1),mai = c(0.2, 0.4, 0.4, 0.4),mar=c(3,4.5,3, 0.2))
par(mfrow= c(4,1),mar=c(3,4,2, 0.2))
kt=paste(expression(italic("k"[t])),"")
plot(type="h",xlab="",ylim=c(0,1),main=expression("Mother Sequence"),ylab="",mean_mother,col="red",xaxt="n",cex.lab=2.6)
#plot(type="h",xlab="",ylim=c(0,1),main=expression("Daughter Sequence"),ylab="",mean_daughter,col="yellow",xaxt="n",cex.lab=2.6)
plot(type="h",xlab="",ylim=c(0,1),main=expression("Corrected Daughter,"~ italic("k"[t])~"=4,Error=0.0204"),mean_ec_daughter_2[4,],xaxt="n",col="darkgreen",ylab="",cex.lab=1.4)
plot(type="h",xlab="",ylim=c(0,1),mean_ec_daughter_2[5,],xaxt="n",col="purple",ylab="",cex.lab=1.4,main=expression("Corrected Daughter,"~ italic("k"[t])~"=5, Error=0.0200"))
plot(type="h",xlab="",ylim=c(0,1),mean_ec_daughter_2[6,],col="lightblue",xaxt="n",ylab="",cex.lab=1.4,main=expression("Corrected Daughter,"~ italic("k"[t])~"=6 , Error=0.0203"))

mtext("Nucleosome Positions in chr19 (bp:7,047,900-8,713,217)",side=1,las=1,line=1,cex=0.8)

mtext("H3K4me3 Occupancy",side=2,las=3,line=-1.2, cex=1,outer=TRUE)

