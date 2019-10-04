library(lavaan)
library(regem)
library(lavaan.survey)
library(psych)

library(seminr)
library(plspm)
#df of lavaan mod

load('C:/Users/Administrator/Desktop/sem/src/50_2436_2000_20_0_11.RData')
load('C:/Users/Administrator/Desktop/sem/src/50_2436_2000_20_1_11.RData')
load('C:/Users/Administrator/Desktop/sem/src/50_8796_2000_20_0_12.RData')
load('C:/Users/Administrator/Desktop/sem/src/50_8796_2000_20_1_12.RData')

mod_num=12



mod=c()
mod[1]<-"F1=~A1+A2+A3+A4+A5+C1+C2+C3+C4+C5+E1+E2+E3+E4+E5+N1+N2+N3+N4+N5+O1+O2+O3+O4+O5
  F2=~C1+C2+C3+C4+C5+A1+A2+A3+A4+A5+E1+E2+E3+E4+E5+N1+N2+N3+N4+N5+O1+O2+O3+O4+O5
  F3=~E1+E2+E3+E4+E5+A1+A2+A3+A4+A5+C1+C2+C3+C4+C5+N1+N2+N3+N4+N5+O1+O2+O3+O4+O5
  F4=~N1+N2+N3+N4+N5+A1+A2+A3+A4+A5+C1+C2+C3+C4+C5+E1+E2+E3+E4+E5+O1+O2+O3+O4+O5
  F5=~O1+O2+O3+O4+O5+A1+A2+A3+A4+A5+C1+C2+C3+C4+C5+E1+E2+E3+E4+E5+N1+N2+N3+N4+N5"
dat=list()
data(bfi)
dat[[1]]=bfi[complete.cases(bfi[,1:25]),1:25]


mod[2] <- "
  math =~ PV1MATH1 + PV1MATH2 + PV1MATH3 + PV1MATH4
  neg.efficacy =~ ST31Q01 + ST31Q02 + ST31Q03 + ST31Q04 +
  ST31Q05 + ST31Q06 + ST31Q07 + ST31Q08
  neg.selfconcept =~ ST32Q02 + ST32Q04 + ST32Q06 + ST32Q07 + ST32Q09
  neg.selfconcept ~ neg.efficacy + ESCS + male
  neg.efficacy ~ neg.selfconcept + school.type + ESCS + male
  math ~ neg.selfconcept + neg.efficacy + school.type + ESCS + male
  "
data(pisa.be.2003)
dat[[2]]=pisa.be.2003
if(mod_num==11){
  lmod<<-sem(mod[1],data=dat[[1]],meanstructure=TRUE)
}else{
  lmod<<-sem(mod[2], data = pisa.be.2003,
             meanstructure = TRUE)
}

#lmod@Fit@test[[1]]$df

#df lavaan.df + number of zeros in coefficients


#p=model$nvar or width of sample cov matrix
#nvar = lmod@pta$nvar[[1]][1]
cbic_real<-function(a,n,mode=-1){
  seq0=exp(  seq(log(0.5),log(0.01),length=30 ) )
  precision=5e-4
  x=a$coef
  if(length(x)==0){
    return (Inf)
  }
  x[a$pen][abs(x[a$pen])<precision]=0
  #nvar=a$nx-sum( abs(x[a$pen])==0)
  nvar=length(a$pen)-sum( abs(x[a$pen])<5e-4)
  npenvar=length(a$pen)-sum( abs(x[a$pen])==0)
  if(mode==-1){
    return( 2*n*a$fit+log(n)*nvar )
  }
  if(mode==0){
    return( n*a$fit  +log(n)*nvar )
  }
  if(mode>0){
    return( 2*n*(a$fit-0.5*seq0[mode]*sum(abs(x[a$pen])))  +log(n)*nvar )
  }
}
 
if(mod_num==11){
  NDATA=2436
  bicvec=c()
  seq0=exp(  seq(log(0.5),log(0.01),length=30 ) )
  for(i in c(1:6,8)){
    if(i<7){
      for(j in 1:30){
      
        bicvec=c(bicvec,cbic_real(result1[[i]][[j]],NDATA,-1))
      }
    }
    if(i==7){
      for(j in 1:30){
        
        bicvec=c(bicvec,cbic_real(result1[[i]][[j]],NDATA,1))
      }
    }
    if(i==8){
      for(j in 1:30){
        
        bicvec=c(bicvec,cbic_real(result1[[i]][[j]],NDATA,0))
      }
    }
  }
}

if(mod_num==12){
  NDATA=8796
  bicvec=c()
  seq0=exp(  seq(log(0.5),log(0.01),length=30 ) )
  for(i in 1:6){
    
      for(j in 1:30){
        
        bicvec=c(bicvec,cbic_real(result1[[i]][[j]],NDATA,-1))
      
    }
  }
  
}

mat=matrix(bicvec,nrow=30)

N= NDATA

#fit =  0.5*(log(det(ImpCov)) + trace(SampCov %*% solve(ImpCov)) - log(det(SampCov))  - p)

chisq = fit*N*2

fm=fitmeasures(lmod)

baseline.chisq = fm["baseline.chisq"]

baseline.df = fm["baseline.df"]

fm['unrestricted.logl']



lavaan.df= fm['df']

df = lavaan.df + sum(abs(x[pars_pen_mod])<5e-4)


chisq = fit*N*2


#ncp
d = function(chisq,df,N){
  max(0,(chisq -df)/(N-1))
}
ncp = d(chisq,df,N)

rmsea = function(ncp,df){
  sqrt(ncp/df)
}
rmsea = max(rmsea(ncp,df),0)


ncp.null = d(baseline.chisq,baseline.df,N)
CFI <- function(ncp.null,ncp){
  (ncp.null - ncp) / ncp.null
}

cfi=CFI(ncp.null,ncp)

TLI <- function(baseline.chisq,baseline.df,chisq,df){#NNFI
  (baseline.chisq/baseline.df - chisq/df) / (baseline.chisq/baseline.df - 1)
}
 TLI(baseline.chisq,baseline.df,chisq,df)
 
npar=length(coef) - sum(abs(x[pars_pen_mod])<5e-4)
bic=2*(N*fit) + log(N)*npar

aic=2*(N*fit) + 2*npar


imp = cov2cor(ImpCov);obs = cov2cor(SampCov2)
lobs <-  obs[!lower.tri(obs)]
limp <-  imp[!lower.tri(imp)]
ret["srmr"] <- sqrt(mean((limp - lobs)^2))




sig=r1$Imp_Cov
S=r1$SampCov
dim(S)
log(det(sig))-
log(det(S))+
sum(diag(S%*%solve(sig)))-
nrow(sig)

A=S%*%solve(sig)
sum(diag(A))-log(det(A))-nrow(A)
rcond(sig)
svd(sig)$d
err=solve(sig)%*%sig-diag(18)
sum(abs(err))
