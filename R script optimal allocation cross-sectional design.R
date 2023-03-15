### load alabama library
library(alabama)

########################################################################################################
### function to calculate variance of treatment effect and functions for equality and inequality constraints
########################################################################################################

### Function to calculate variance of treatment effect
f.vartreat=function(x)
{
  nr.seq=length(x) # number of sequences
  nr.per=nr.seq+1 # number of periods
  exponent <- abs(matrix(1:nr.per - 1, nrow = nr.per, ncol = nr.per, byrow = TRUE) - (1:nr.per - 1))
  VV <- ICC*CAC^exponent + diag(1-ICC,nrow=nr.per)/nr.subj # correlation matrix for cross-sectional
  VVinv=solve(VV)
  XX=cbind(rep(1,nr.per),diag(1,nrow=nr.per)[,2:nr.per],rep(1,nr.per))
  inf.mat=matrix(0,nr.per+1,nr.per+1)
  for(ii in 1:nr.seq)
  {
    XX[ii,(nr.per+1)] <- 0
    inf.mat=inf.mat+x[ii]*t(XX)%*%VVinv%*%XX
  }
   cov.mat=solve(inf.mat)
   cov.mat[(nr.per+1),(nr.per+1)]
}

### Equality constraint function
### Proportions should sum to 1
heq=function(x)
  sum(x)-1

### Inequality constraint function
### All proportions between 0 and 1
hin <- function(x)
{
  nr.seq <- length(x)
  h <- rep(NA, 2*nr.seq)
  h[1:nr.seq] <- 1-x            
  h[(nr.seq+1):(2*nr.seq)] <- x
  h
}

########################################################################################################
### plots for three sequences
########################################################################################################
nr.seq=3
nr.per=nr.seq+1
nr.subj=5
ICC=0.0125
  
CAC.vec=seq(0,1,l=51)
  
seq1=rep(NA,51)
seq2=rep(NA,51)
seq3=rep(NA,51)
RE=rep(NA,51)
  
for(ii in 1:51)
  {
    CAC=CAC.vec[ii] 
    start=rep(1/nr.seq,nr.seq)
    solution=constrOptim.nl(par=start,fn=f.vartreat,heq=heq,hin=hin)
    seq1[ii]=solution$par[1]
    seq2[ii]=solution$par[2]
    seq3[ii]=solution$par[3]
    RE[ii]=f.vartreat(solution$par)/f.vartreat(start)
    }
  

plot(CAC.vec,seq1,ylim=c(0,1),cex.axis=1.25,cex.lab=1.25,cex=2,xlab="Cluster autocorrelation",ylab="Proportion",type="l",lty=1) 
lines(CAC.vec,seq2,lty=2)
lines(CAC.vec,seq3,lty=1)
legend(0,1,legend=c("Sequences 1 and 3","Sequence 2"),lty=c(1,2),box.lty=0,cex=1.25)
plot(CAC.vec,RE,ylim=c(0,1),cex.axis=1.25,cex.lab=1.25,cex=2,xlab="Cluster autocorrelation",ylab="Relative efficiency",type="l") 


########################################################################################################
### plots for four sequences
########################################################################################################
nr.seq=4
nr.per=nr.seq+1
nr.subj=5
ICC=0.05
  
CAC.vec=seq(0,1,l=51)
  
seq1=rep(NA,51)
seq2=rep(NA,51)
seq3=rep(NA,51)
seq4=rep(NA,51)
RE=rep(NA,51)
  
for(ii in 1:51)
  {
    CAC=CAC.vec[ii] 
    start=rep(1/nr.seq,nr.seq)
    solution=constrOptim.nl(par=start,fn=f.vartreat,heq=heq,hin=hin)
    seq1[ii]=solution$par[1]
    seq2[ii]=solution$par[2]
    seq3[ii]=solution$par[3]
    seq4[ii]=solution$par[4]
    RE[ii]=f.vartreat(solution$par)/f.vartreat(start)
    }
  

plot(CAC.vec,seq1,ylim=c(0,1),cex.axis=1.25,cex.lab=1.25,cex=2,xlab="Cluster autocorrelation",ylab="Proportion",type="l",lty=1) 
lines(CAC.vec,seq2,lty=2)
lines(CAC.vec,seq3,lty=2)
lines(CAC.vec,seq4,lty=3)
legend(0,1,legend=c("Sequences 1 and 4","Sequences 2 and 3"),lty=c(1,2),box.lty=0,cex=1.25)
plot(CAC.vec,RE,ylim=c(0,1),cex.axis=1.25,cex.lab=1.25,cex=2,xlab="Cluster autocorrelation",ylab="Relative efficiency",type="l") 



########################################################################################################
### plots for five sequences
########################################################################################################
nr.seq=5
nr.per=nr.seq+1
nr.subj=5
ICC=0.05
  
CAC.vec=seq(0,1,l=51)
  
seq1=rep(NA,51)
seq2=rep(NA,51)
seq3=rep(NA,51)
seq4=rep(NA,51)
seq5=rep(NA,51)
RE=rep(NA,51)
  
for(ii in 1:51)
  {
    CAC=CAC.vec[ii] 
    start=rep(1/nr.seq,nr.seq)
    solution=constrOptim.nl(par=start,fn=f.vartreat,heq=heq,hin=hin)
    seq1[ii]=solution$par[1]
    seq2[ii]=solution$par[2]
    seq3[ii]=solution$par[3]
    seq4[ii]=solution$par[4]
    seq5[ii]=solution$par[5]
    RE[ii]=f.vartreat(solution$par)/f.vartreat(start)
    }
  

plot(CAC.vec,seq1,ylim=c(0,1),cex.axis=1.25,cex.lab=1.25,cex=2,xlab="Cluster autocorrelation",ylab="Proportion",type="l",lty=1) 
lines(CAC.vec,seq2,lty=2)
lines(CAC.vec,seq3,lty=3)
lines(CAC.vec,seq4,lty=2)
lines(CAC.vec,seq5,lty=1)
legend(0,1,legend=c("Sequences 1 and 5","Sequences 2 and 4","Sequence 3"),lty=c(1,2,3),box.lty=0,cex=1.25)
plot(CAC.vec,RE,ylim=c(0,1),cex.axis=1.25,cex.lab=1.25,cex=2,xlab="Cluster autocorrelation",ylab="Relative efficiency",type="l") 


########################################################################################################
### plots for six sequences
########################################################################################################
nr.seq=6
nr.per=nr.seq+1
nr.subj=25
ICC=0.05
  
CAC.vec=seq(0,1,l=51)
  
seq1=rep(NA,51)
seq2=rep(NA,51)
seq3=rep(NA,51)
seq4=rep(NA,51)
seq5=rep(NA,51)
seq6=rep(NA,51)
RE=rep(NA,51)
  
for(ii in 1:51)
  {
    CAC=CAC.vec[ii] 
    start=rep(1/nr.seq,nr.seq)
    solution=constrOptim.nl(par=start,fn=f.vartreat,heq=heq,hin=hin)
    seq1[ii]=solution$par[1]
    seq2[ii]=solution$par[2]
    seq3[ii]=solution$par[3]
    seq4[ii]=solution$par[4]
    seq5[ii]=solution$par[5]
    seq6[ii]=solution$par[6]
    RE[ii]=f.vartreat(solution$par)/f.vartreat(start)
    }
  

plot(CAC.vec,seq1,ylim=c(0,1),cex.axis=1.25,cex.lab=1.25,cex=2,xlab="Cluster autocorrelation",ylab="Proportion",type="l",lty=1) 
lines(CAC.vec,seq2,lty=2)
lines(CAC.vec,seq3,lty=3)
lines(CAC.vec,seq4,lty=3)
lines(CAC.vec,seq5,lty=2)
lines(CAC.vec,seq6,lty=1)
legend(0,1,legend=c("Sequences 1 and 6","Sequences 2 and 5","Sequences 3 and 4"),lty=c(1,2,3),box.lty=0,cex=1.25)
plot(CAC.vec,RE,ylim=c(0,1),cex.axis=1.25,cex.lab=1.25,cex=2,xlab="Cluster autocorrelation",ylab="Relative efficiency",type="l") 







