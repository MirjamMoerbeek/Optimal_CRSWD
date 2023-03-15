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
  VV <- ICC*CAC^exponent + ((1-ICC)*IAC^exponent)/nr.subj # correlation matrix for cohort design
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
### Upper and lower bounds for proportions 
hin <- function(x)
{
  nr.seq <- length(x)
  h <- rep(NA, 2*nr.seq)
  h[1:nr.seq] <- 1-x            
  h[(nr.seq+1):(2*nr.seq)] <- x
  h
}


########################################################################################################
### contour plots for three sequences
########################################################################################################
nr.seq=3
nr.per=nr.seq+1
nr.subj=50
ICC=0.05
  
CAC.vec=seq(0,1,l=51)
IAC.vec=seq(0,1,l=51)
  
seq1=matrix(NA,51,51)
seq2=matrix(NA,51,51)
seq3=matrix(NA,51,51)
RE=matrix(NA,51,51)
  
for(ii in 1:51)
  {
  for(jj in 1:51)
    {
    print(c(ii,jj))
    CAC=CAC.vec[ii] 
    IAC=IAC.vec[jj] 
     
    start=rep(1/nr.seq,nr.seq)
    solution=constrOptim.nl(par=start,fn=f.vartreat,heq=heq,hin=hin)
    seq1[ii,jj]=solution$par[1]
    seq2[ii,jj]=solution$par[2]
    seq3[ii,jj]=solution$par[3]
    RE[ii,jj]=f.vartreat(solution$par)/f.vartreat(start)
    }
  }

contour(CAC.vec,IAC.vec,seq1,levels=seq(0,1,l=51),cex.axis=1.25,cex.lab=1.25,cex=2,labcex=1.125,xlab="Cluster autocorrelation",ylab="Individual autocorrelation",main="Proportion in seqences 1 and 3")
contour(CAC.vec,IAC.vec,seq2,levels=seq(0,1,l=51),cex.axis=1.25,cex.lab=1.25,cex=2,labcex=1.125,xlab="Cluster autocorrelation",ylab="Individual autocorrelation",main="Proportion in seqence 2")
contour(CAC.vec,IAC.vec,RE,levels=seq(0,1,l=51),cex.axis=1.25,cex.lab=1.25,cex=2,labcex=1.125,xlab="Cluster autocorrelation",ylab="Individual autocorrelation",main="Relative efficiency of uniform design")

########################################################################################################
### contour plots for four sequences
########################################################################################################

nr.seq=4
nr.per=nr.seq+1
nr.subj=50
ICC=0.05

CAC.vec=seq(0,1,l=51)
IAC.vec=seq(0,1,l=51)

seq1=matrix(NA,51,51)
seq2=matrix(NA,51,51)
seq3=matrix(NA,51,51)
seq4=matrix(NA,51,51)
RE=matrix(NA,51,51)

for(ii in 1:51)
{
  for(jj in 1:51)
  {
    print(c(ii,jj))
    CAC=CAC.vec[ii] 
    IAC=IAC.vec[jj] 
    
    start=rep(1/nr.seq,nr.seq)
    solution=constrOptim.nl(par=start,fn=f.vartreat,heq=heq,hin=hin)
    seq1[ii,jj]=solution$par[1]
    seq2[ii,jj]=solution$par[2]
    seq3[ii,jj]=solution$par[3]
    seq4[ii,jj]=solution$par[4]
    RE[ii,jj]=f.vartreat(solution$par)/f.vartreat(start)
  }
}

contour(CAC.vec,IAC.vec,seq1,levels=seq(0,1,l=51),cex.axis=1.25,cex.lab=1.25,cex=2,labcex=1.125,xlab="Cluster autocorrelation",ylab="Individual autocorrelation",main="Proportion in seqences 1 and 4")
contour(CAC.vec,IAC.vec,seq2,levels=seq(0,1,l=51),cex.axis=1.25,cex.lab=1.25,cex=2,labcex=1.125,xlab="Cluster autocorrelation",ylab="Individual autocorrelation",main="Proportion in seqences 2 and 3")
contour(CAC.vec,IAC.vec,RE,levels=seq(0,1,l=51,cex.axis=1.25,cex.lab=1.25),cex=2,labcex=1.125,xlab="Cluster autocorrelation",ylab="Individual autocorrelation",main="Relative efficiency of uniform design")

########################################################################################################
### contour plots for five sequences
########################################################################################################

nr.seq=5
nr.per=nr.seq+1
nr.subj=25
ICC=0.025

CAC.vec=seq(0,1,l=51)
IAC.vec=seq(0,1,l=51)

seq1=matrix(NA,51,51)
seq2=matrix(NA,51,51)
seq3=matrix(NA,51,51)
seq4=matrix(NA,51,51)
seq5=matrix(NA,51,51)
RE=matrix(NA,51,51)

for(ii in 1:51)
  {
  for(jj in 1:51)
    {
    print(c(ii,jj))
    CAC=CAC.vec[ii] 
    IAC=IAC.vec[jj] 
    
    start=rep(1/nr.seq,nr.seq)
    solution=constrOptim.nl(par=start,fn=f.vartreat,heq=heq,hin=hin)
    seq1[ii,jj]=solution$par[1]
    seq2[ii,jj]=solution$par[2]
    seq3[ii,jj]=solution$par[3]
    seq4[ii,jj]=solution$par[4]
    seq5[ii,jj]=solution$par[5]
    RE[ii,jj]=f.vartreat(solution$par)/f.vartreat(start)
  }
}

contour(CAC.vec,IAC.vec,seq1,levels=seq(0,1,l=51),cex.axis=1.25,cex.lab=1.25,cex=2,labcex=1.125,xlab="Cluster autocorrelation",ylab="Individual autocorrelation",main="Proportion in seqences 1 and 5")
contour(CAC.vec,IAC.vec,seq2,levels=seq(0,1,l=51),cex.axis=1.25,cex.lab=1.25,cex=2,labcex=1.125,xlab="Cluster autocorrelation",ylab="Individual autocorrelation",main="Proportion in seqences 2 and 4")
contour(CAC.vec,IAC.vec,seq3,levels=seq(0,1,l=51),cex.axis=1.25,cex.lab=1.25,cex=2,labcex=1.125,xlab="Cluster autocorrelation",ylab="Individual autocorrelation",main="Proportion in seqence 3")
contour(CAC.vec,IAC.vec,RE,levels=seq(0,1,l=51),cex.axis=1.25,cex.lab=1.25,cex=2,labcex=1.125,xlab="Cluster autocorrelation",ylab="Individual autocorrelation",main="Relative efficiency of uniform design")

########################################################################################################
### contour plots for six sequences
########################################################################################################

nr.seq=6
nr.per=nr.seq+1
nr.subj=25
ICC=0.05

CAC.vec=seq(0,1,l=51)
IAC.vec=seq(0,1,l=51)

seq1=matrix(NA,51,51)
seq2=matrix(NA,51,51)
seq3=matrix(NA,51,51)
seq4=matrix(NA,51,51)
seq5=matrix(NA,51,51)
seq6=matrix(NA,51,51)
RE=matrix(NA,51,51)

for(ii in 1:51)
  {
  for(jj in 1:51)
   {
    print(c(ii,jj))
    CAC=CAC.vec[ii] 
    IAC=IAC.vec[jj] 
    
    start=rep(1/nr.seq,nr.seq)
    solution=constrOptim.nl(par=start,fn=f.vartreat,heq=heq,hin=hin)
    seq1[ii,jj]=solution$par[1]
    seq2[ii,jj]=solution$par[2]
    seq3[ii,jj]=solution$par[3]
    seq4[ii,jj]=solution$par[4]
    seq5[ii,jj]=solution$par[5]
    seq6[ii,jj]=solution$par[6]
    RE[ii,jj]=f.vartreat(solution$par)/f.vartreat(start)
  }
}


contour(CAC.vec,IAC.vec,seq1,levels=seq(0,1,l=51),cex.axis=1.25,cex.lab=1.25,cex=2,labcex=1.125,xlab="Cluster autocorrelation",ylab="Individual autocorrelation",main="Proportion in seqences 1 and 6")
contour(CAC.vec,IAC.vec,seq2,levels=seq(0,1,l=51),cex.axis=1.25,cex.lab=1.25,cex=2,labcex=1.125,xlab="Cluster autocorrelation",ylab="Individual autocorrelation",main="Proportion in seqences 2 and 5")
contour(CAC.vec,IAC.vec,seq3,levels=seq(0,1,l=51),cex.axis=1.25,cex.lab=1.25,cex=2,labcex=1.125,xlab="Cluster autocorrelation",ylab="Individual autocorrelation",main="Proportion in seqences 3 and 4")
contour(CAC.vec,IAC.vec,RE,levels=seq(0,1,l=51),cex.axis=1.25,cex.lab=1.25,cex=2,labcex=1.125,xlab="Cluster autocorrelation",ylab="Individual autocorrelation",main="Relative efficiency of uniform design")


