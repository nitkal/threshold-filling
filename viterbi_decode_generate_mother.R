###########################
#Program for SMAP decoding using viterbi 
# for generated mothers (of different alpha, beta)
# and computing mean error
# Generates heatmap and also error plots in supplementary
#######################

source("../compute_alpha.R")
source("viterbi_decode_lookup.R")
library("caroline")
library("Hmisc")
library("tidyverse")
alphas = seq(from = 0.1, to = 0.9, length.out = 9 )
betas =  seq(from = 0.1, to = 0.9, length.out = 9 )
no_of_mothers = 300
no_of_daughters = 200

mean_error = matrix(ncol = 9, nrow = 9)
mean_ec_error = matrix(ncol = 9, nrow = 9)
sd_ec_error = matrix(ncol = 9, nrow = 9)
errors_bc = array(numeric(),c(no_of_mothers,no_of_daughters))
errors_ac = array(numeric(),c(no_of_mothers,no_of_daughters))
if(file.exists("viterbi_pattern_df.rda")) {
  load(file="viterbi_pattern_df.rda")
  
}

#library(binaryLogic)
for (h in 1:length(alphas) ) {
  for (g in 1:length(betas) ) {
    
    alpha = alphas[h]
    beta = betas[g]
    length_of_seq = 100
    #alpha = alphas[h]
    
    for (mo in 1: no_of_mothers) {
      mother_chromatin = vector(length=length_of_seq)
      error = vector(length=no_of_daughters)
      ec_error = vector(length = no_of_daughters)
      #no_of_zeros_correctly_filled  = vector(length = no_of_iterations)
      
      
      mother_chromatin[1] = 1
      for (i in 2:length_of_seq) {
        u = runif(1,0,1)
        if(mother_chromatin[i-1] == 1) {
          mother_chromatin[i] = if(u <= alpha) 1 else 0
        }
        else {
          mother_chromatin[i] = if(u <= (beta)) 0 else 1
        }
      }
      #compute_alpha(mother_chromatin)
      daughter_chromatins = matrix(ncol=length_of_seq)
      daughter_chromatins = daughter_chromatins[-1,]
      
      
      
      for (k in 1:no_of_daughters) {
        prob_of_copy = rbinom(length_of_seq,1,0.5)
        
        daughter_chromatin = vector(length = length_of_seq)
        daughter_chromatin[1] = mother_chromatin[1]
        for (i in 2:length_of_seq) {
          
          daughter_chromatin[i] = ifelse(((mother_chromatin[i] == 1) && (prob_of_copy[i] == 1)), mother_chromatin[i] ,0)
          
        }
        daughter_chromatins = rbind(daughter_chromatins,daughter_chromatin)
        errors_bc[mo,k] = sum(xor(mother_chromatin,daughter_chromatins[k,]))/length_of_seq
        
      }
      
      #Viterbi Decoding
      
      ec_daughter_chromatins = matrix(ncol=length_of_seq)
      ec_daughter_chromatins =ec_daughter_chromatins[-1,]
      
      
      
      for (k in 1:no_of_daughters) {
        y = daughter_chromatins[k,]
        if(exists("viterbi.pattern.df")) {
          yy =  viterbi_decode_lookup(y,alpha,beta,viterbi.pattern.df)
        }
        else {
          yy =  viterbi_decode_lookup(y,alpha,beta,NULL)
        }
        
        ec_daughter_chromatin = yy$y
        viterbi.pattern.df = yy$viterbi.pattern.df
        ec_daughter_chromatins = rbind(ec_daughter_chromatins,ec_daughter_chromatin)
        errors_ac[mo,k] = sum(xor(mother_chromatin,ec_daughter_chromatins[k,]))/length_of_seq
        # no_of_zeros_correctly_filled[k] = 0
        # for( t in 1:length_of_seq) {
        #   if ((mother_chromatin[t] == 1) && (daughter_chromatins[k,t] ==0) && (ec_daughter_chromatins[k,t] == 1))
        #   {
        #     no_of_zeros_correctly_filled[k] = no_of_zeros_correctly_filled[k] + 1
        #   }
        # }
      }
    }
    
    mean_error[h,g] = mean(errors_bc)
    mean_ec_error[h,g] = mean(errors_ac)
    sd_ec_error[h,g] = sd(errors_ac)
    # mean_no_of_zeros_correctly_filled[h,g] = mean(no_of_zeros_correctly_filled)
  }
}

#save(mean_ec_error,file="map_decode_mean_ec_error.rda")

# ymin = min(c(mean_ec_error[1,],mean_ec_error[3,],mean_ec_error[5,],mean_ec_error[7,],mean_ec_error[9,]))
# ymax = max(c(mean_ec_error[1,],mean_ec_error[3,],mean_ec_error[5,],mean_ec_error[7,],mean_ec_error[9,])) +0.1
# x = alphas
# y1 = mean_ec_error[,1]
# y2 = mean_ec_error[,3]
# y3 = mean_ec_error[,5]
# y4 = mean_ec_error[,7]
# y5 = mean_ec_error[,9]
# lo1 = smooth.spline(x,y1,spar=0.3)
# lo2 = smooth.spline(x,y2,spar=0.3)
# lo3 = smooth.spline(x,y3,spar=0.3)
# lo4 = smooth.spline(x,y4,spar=0.3)
# lo5 = smooth.spline(x,y5,spar=0.3)
# plot(alphas,y1,ylim=c(ymin,ymax),type="n",pch = 1, xlab=expression(alpha), ylab = "Mean error after correction", cex.lab=1.7, cex.axis=1.5, cex.main=1.2, cex.sub=1.5,col="red")
# lines(lo1,lwd="5",col="red",type="b")
# lines(lo2,lwd="5", col="blue", pch=1,type = "b")
# lines((lo3),lwd="5", col="darkorchid", pch=1,type="b")
# lines((lo4), lwd="5",col="brown", pch=1,type="b")
# lines((lo5), lwd="5",col="orange", pch=1,type="b")
# 
# 
# #legend(2,0.5, "Beta=0.5",box.lwd = 0,box.col = "white",cex=1.2)
# legend(0.7,0.4, legend=c(expression(paste(beta,"=0.1")),expression(paste(beta,"=0.3")),expression(paste(beta,"=0.5")),expression(paste(beta,"=0.7")),expression(paste(beta,"=0.9"))),col=c("red", "blue","darkorchid","brown","orange"),box.lwd = 0,box.col = "white",cex=0.8,lwd="5")
# box(lwd=5)
mean_ec_error1 = mean_ec_error

#betas = rev(betas)
#mean_ec_error1 = mean_ec_error1[nrow(mean_ec_error1):1, ]
rownames(mean_ec_error) <- alphas
colnames(mean_ec_error) <- betas
colours <- colorRampPalette(c("red","blue"))

#levelplot(mean_ec_error1, col.regions=colours,xlab=list(paste(expression(alpha) ,":", expression(P(m[i]),"|",expression(m[i-1])," = 1")),cex=1.8),ylab=list(expression(paste(beta," - P(m_i|m_(i-1)) = 0")),cex=1.8),scales=list(x=list(cex=1.2),y=list(cex=1.2)))

#xlabel = bquote(.sym(alpha)~"Probability of maintianing modified state")
#ylabel = bquote(.sym(beta)~"Probability of maintianing unmodified state")

xlabel = expression(paste(alpha," - Probability of maintaining modified state"))
ylabel = expression(paste(beta," - Probability of maintaining unmodified state"))

levelplot(mean_ec_error, col.regions=colours,xlab=list(xlabel,cex=1.5),ylab=list(ylabel,cex=1.5),scales=list(x=list(cex=1.2),y=list(cex=1.2)))



sem_ec_error = sd_ec_error / sqrt(no_of_mothers*no_of_daughters)
#sem_ec_error = sd_ec_error
par( oma = c(0.2,0.2,0.1,0.1),mfrow=c(1,2))
y1 = mean_ec_error[,2]
y2 = mean_ec_error[,4]
y3 = mean_ec_error[,6]
y4 = mean_ec_error[,8]
y5 = mean_ec_error[,9]

sem1 = sem_ec_error[,2]
sem2 = sem_ec_error[,4]
sem3 = sem_ec_error[,6]
sem4 = sem_ec_error[,8]
sem5 = sem_ec_error[,9]

ymin = min(c(y1,y2,y3,y4,y5))
ymax = max(c(y1,y2,y3,y4,y5)+0.1)
x = alphas

plot(alphas,y1,ylim=c(ymin,ymax),type="n",pch = 0, xlab=xlabel, ylab = "Mean error after correction", cex.lab=1.4, cex.axis=1.6,  cex.sub=1.6,col="red")
lines(alphas,y1,lwd="2",col="red",type="b",pch=19,lty=1)
lines(alphas,y2,lwd=2, col="blue", pch=19,type = "b")
lines(alphas,y3,lwd=2, col="darkorchid", pch=19,type="b")
lines(alphas,y4, lwd=2,col="brown", pch=19,type="b")
lines(alphas,y5, lwd=2,col="darkgreen", pch=19,type="b")
par(fg="red")
Hmisc::errbar(alphas,y1,y1+sem1,y1-sem1,add=TRUE,errbar.col="red",lwd=3)
par(fg="blue")
Hmisc::errbar(alphas,y2,y2+sem2,y2-sem2,add=TRUE,errbar.col="blue",lwd=3)
par(fg="darkorchid")
Hmisc::errbar(alphas,y3,y3+sem3,y3-sem3,add=TRUE,errbar.col="darkorchid",lwd=3)
par(fg="brown")
Hmisc::errbar(alphas,y4,y4+sem4,y4-sem4,add=TRUE,errbar.col="brown",lwd=3)
par(fg="darkgreen")
Hmisc::errbar(alphas,y5,y5+sem5,y5-sem5,add=TRUE,errbar.col="darkgreen",lwd=3)

par(fg="black")
#legend(2,0.5, "Beta=0.5",box.lwd = 0,box.col = "white",cex=1.2)
legend(0.1,0.435, legend=c(expression(paste(beta," = 0.2")),expression(paste(beta," = 0.4")),expression(paste(beta," = 0.6")),expression(paste(beta," = 0.8")),expression(paste(beta," = 0.9"))),col=c("red", "blue","darkorchid","brown","darkgreen"),box.lwd = 0,box.col = "white",cex=0.9,lwd="2")
par(xpd=TRUE)
legend(-0.1,0.5,legend="(C)",cex=1.6,box.col="white")
box(lwd=3)



y11 = mean_ec_error[2,]
y22 = mean_ec_error[4,]
y33 = mean_ec_error[6,]
y44 = mean_ec_error[8,]
y55 = mean_ec_error[9,]

sem11 = sem_ec_error[2,]
sem22 = sem_ec_error[4,]
sem33 = sem_ec_error[6,]
sem44 = sem_ec_error[8,]
sem55 = sem_ec_error[9,]

ymin = min(c(y11,y22,y33,y44,y55))
ymax = max(c(y11,y22,y33,y44,y55)+0.1)
x = alphas

plot(alphas,y11,ylim=c(ymin,ymax),type="n",pch = 0, xlab=ylabel, ylab = "Mean error after correction", cex.lab=1.4, cex.axis=1.6,  cex.sub=1.6,col="red")
lines(alphas,y11,lwd="2",col="red",type="b",pch=19,lty=1)
lines(alphas,y22,lwd=2, col="blue", pch=19,type = "b")
lines(alphas,y33,lwd=2, col="darkorchid", pch=19,type="b")
lines(alphas,y44, lwd=2,col="brown", pch=19,type="b")
lines(alphas,y55, lwd=2,col="darkgreen", pch=19,type="b")
par(fg="red")
Hmisc::errbar(alphas,y11,y11+sem11,y11-sem11,add=TRUE,errbar.col="red",lwd=3)
par(fg="blue")
Hmisc::errbar(alphas,y22,y22+sem22,y22-sem22,add=TRUE,errbar.col="blue",lwd=3)
par(fg="darkorchid")
Hmisc::errbar(alphas,y33,y33+sem33,y33-sem33,add=TRUE,errbar.col="darkorchid",lwd=3)
par(fg="brown")
Hmisc::errbar(alphas,y44,y44+sem44,y44-sem44,add=TRUE,errbar.col="brown",lwd=3)
par(fg="darkgreen")
Hmisc::errbar(alphas,y55,y55+sem55,y55-sem55,add=TRUE,errbar.col="darkgreen",lwd=3)

par(fg="black")
#legend(2,0.5, "Beta=0.5",box.lwd = 0,box.col = "white",cex=1.2)
legend(0.5,0.435, legend=c(expression(paste(alpha," = 0.2")),expression(paste(alpha," = 0.4")),expression(paste(alpha," = 0.6")),expression(paste(alpha," = 0.8")),expression(paste(alpha," = 0.9"))),col=c("red", "blue","darkorchid","brown","darkgreen"),box.lwd = 0,box.col = "white",cex = 0.9,lwd="2")
par(xpd=TRUE)
legend(-0.1,0.5,legend="(D)",cex=1.6,box.col="white")
box(lwd=3)

