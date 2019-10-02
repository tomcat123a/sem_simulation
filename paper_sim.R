args <- commandArgs(trailingOnly=TRUE)
print(args)
MYSEED=as.integer(args[1])
NDATA=as.integer(args[2])
maxit=as.integer(args[3])
mod_num=as.integer(args[4])
VFAC=as.double(args[5])
simumode=as.integer(args[6])
#MYSEED=0
#NDATA=300
#maxit=2
#mod_num=3
#VFAC=1
#simumode=2
SEQU=exp(  seq(log(0.5),log(0.01),length=30 ) )
const_para_gradsamp_ra=1e-3

set.seed(MYSEED)

randcoef=(2*rbinom(size=1,n=21,prob=0.5)-1) * runif(21,min=0.5,max=1)

initcoef=(2*rbinom(size=1,n=25,prob=0.5)-1) * runif(25,min=0.5,max=1)



covcoef=runif(30,min=0.1,max=0.5)
if(mod_num==3){
  randcoef[c(14,16,18)]=(2*rbinom(size=1,n=3,prob=0.5)-1)*runif(3,min=0.5,max=0.8)
}
if(mod_num==6){
  randcoef[c(14,16)]=(2*rbinom(size=1,n=2,prob=0.5)-1)*runif(2,min=0.5,max=0.8)
}


if(simumode==0){
  Typechar="lasso"
}
if(simumode==1){
  Typechar="scad"
}
if(simumode==2){
  Typechar="mcp"
}

#MYSEED=0
#NDATA=200
#maxit=1000
#mod_num=3
#VFAC=1
#simumode=0



if(!require(lavaan)){
  install_version("lavaan", version = "0.5-22", repos = "http://cran.us.r-project.org")
  library(lavaan)
}

if(!require(regsem)){
  #install_version("regsem", version = "0.9.2", repos = "http://cran.us.r-project.org")
  install.packages("regsem",repos="http://cran.us.r-project.org")
  library(regsem)
}



if(!require(lbfgs)){
  install.packages("lbfgs")
  library(lbfgs)
}
if(!require(methods)){
  install.packages("methods")
  library(methods)
}
if(!require(Rcpp)){
  install.packages("Rcpp")
  library(Rcpp)
}
if(!require(Rsolnp)){
  install.packages("Rsolnp")
  library(Rsolnp)
}
if(!require(psych)){
  install.packages("psych")
  library(psych)
}
if(!require(lsl)){
  install.packages("lsl")
  library(lsl)
}
if(!require(lslx)){
  install.packages("lslx")
  library(lslx)
}
if(!require(lavaan.survey)){
  install.packages("lavaan.survey")
  library(lavaan.survey)
}
#setwd("C:/Users/Administrator/Desktop/sem/src")
sourceCpp("my_constr.cpp")
##########functions start



initx<-function(){
  #set.seed(MYSEED)
  #sx=runif(n,-0.5,0.5)
  #zv=sample(v,ceiling(length(v)*0.8),replace=FALSE)
  #sx[zv]=0
  #sx[v]=0
  #sx[-v]=runif(n-length(v),0.4,0.6)
  #sx
  as.numeric(extractMatrices(lmod)$parameters)
}

SCAD<-function(x){
  if(abs(x)<=l){
    return (abs(x))
  }else{
    if(abs(x)>l & abs(x) < 3.7*l){
      return (  (2*3.7*l*abs(x)-x^2-l^2  )/2/(2.7)  )
    }else{
      return( 2.35*l^2)
    }
  }

}
P1<-function(x,alpha){
  sum(unlist(  lapply(x,SCAD) ) )
}

SCADderi<-function(x){
  if(abs(x)<=l){
    return (sign(x))
  }else{
    if(abs(x)>l & abs(x) < 3.7*l){
      return (  (3.7*l*sign(x)-x  )/(2.7)  )
    }else{
      return( 0)
    }
  }

}
gradscad<-function(x,alpha){
  unlist(  lapply(x,SCADderi) )
}



shuffle<-function(x){
  newx=vector("numeric",length(x))
  newx[pars_pen_mod]=x[1:length(pars_pen_mod)]
  newx[-pars_pen_mod]=x[-c(1:length(pars_pen_mod))]
  newx
}
rev_shuffle<-function(x){
  newx=vector("numeric",length(x))
  newx[1:length(pars_pen_mod)]=x[pars_pen_mod]
  newx[-c(1:length(pars_pen_mod))]=x[-pars_pen_mod]
  newx
}


fobjml<-function(x){
  ###object without penalty
  x=as.double(x)
  mult = regsem:::rcpp_RAMmult(par=x,A,S,S_fixed,A_fixed,A_est,S_est,F_b,I_b)
  pen_vec = x[pars_pen_mod]
  ret=tryCatch(
    regsem:::rcpp_fit_fun(ImpCov=mult$ImpCov,dataS_m,type2=0,l,3.7,pen_vec,pen_diff=0,e_alpha=1)
    ,error=function(error){
      Inf
    })
  if(is.na(ret)){
    ret=Inf
  }
  ret
}

fobjml_0pen<-function(x){
  ###object without penalty
  x=as.double(x)
  mult = regsem:::rcpp_RAMmult(par=x,A,S,S_fixed,A_fixed,A_est,S_est,F_b,I_b)
  pen_vec = x[pars_pen_mod]
  ret=tryCatch(
    regsem:::rcpp_fit_fun(ImpCov=mult$ImpCov,dataS_m,type2=0,0,3.7,pen_vec,pen_diff=0,e_alpha=1)
    ,error=function(error){
      Inf
    })
  if(is.na(ret)){
    ret=Inf
  }
  ret
}
gradml<-function(x){
  ###object without penalty
  x=as.double(x)
  mult = regsem:::rcpp_RAMmult(par=x,A,S,S_fixed,A_fixed,A_est,S_est,F_b,I_b)
  pen_vec = x[pars_pen_mod]
  #mult = regsem:::rcpp_RAMmult(par=x0,A,S,F,A_fixed,A_est,S_fixed,S_est)
  regsem:::rcpp_grad_ram(par=x,ImpCov=mult$ImpCov,SampCov=dataS_m,Areg = mult$A_est22,
                         Sreg=mult$S_est22,A,S,
                         F_b,lambda=l,type2=0,pen_vec=pen_vec,diff_par=0)
}

fl<-function(x){
  #for lbfgs ,shuffle
  x=shuffle(x)
  fobjml(x)
}


gl<-function(x){
  #for lbfgs ,shuffle
  x=shuffle(x)
  gradml(x)
}





Soft<-function(a,b){
  if( b<0 )stop("error in Soft the second parameter must be non-negative")
  if( b>abs(a) ){
    return (0)
  }else{
    return ( a - sign(a)*b )
  }
}

computebic<-function(fit,x,precision=5e-4,special=FALSE){
  #this applies to lbfgs,gradsamp,regsem
  x[pars_pen_mod][abs(x[pars_pen_mod])<precision]=0
  nvar=nx-sum( abs(x[pars_pen_mod])==0)
  npenvar=length(pars_pen_mod)-sum( abs(x[pars_pen_mod])==0  )
  #return( 2*NDATA*fit+(log(NDATA))*nvar+2*log(choose(length(pars_pen_mod),npenvar))  )
  if(special==FALSE){
    return( 2*NDATA*fit+(log(NDATA))*npenvar )
  }else{
    return(  NDATA*fit+(log(NDATA))*npenvar )
  }
}

####################start mod script
L1<-list()#true model
L2<-list()#generate data

L3<-list()#for estimation

#mod3 with f4~f5,f5~f6,f6~f4 in true model
#mod4 semi-confirm f5~f2, f4~f5 unpenalized
#mod7 no f4~f5,f5~f6,f6~f4 in true model
#mod8 semi-confirm f5~f2 unpenalized
#mod9 semi-confirm f5~f2, f4~f5 unpenalized
#mod10
L1[[3]]=paste0("f1=~x1+start(",randcoef[1],")*x2+start(",randcoef[2],")*x3+start(0)*x4+start(0)*x5+start(0)*x6+start(0)*x7+start(0)*x8+start(0)*x9 \n
               f2=~x4+start(",randcoef[3],")*x5+start(",randcoef[4],")*x6 + start(0)*x1 + start(0)*x2 + start(0)*x3+start(0)*x7+start(0)*x8+start(0)*x9 \n
               f3=~x7 + start(",randcoef[5],")*x8 + start(",randcoef[6],")*x9 + start(0)*x1 + start(0)*x2 + start(0)*x3+start(0)*x4+start(0)*x5+start(0)*x6 \n
               f4=~x10 + start(",randcoef[7],")*x11 + start(",randcoef[8],")*x12+start(0)*x13+start(0)*x14+start(0)*x15+start(0)*x16+start(0)*x17+start(0)*x18 \n
               f5=~x13 + start(",randcoef[9],")*x14 + start(",randcoef[10],")*x15+start(0)*x10+start(0)*x11+start(0)*x12+start(0)*x16+start(0)*x17+start(0)*x18 \n
               f6=~x16 + start(",randcoef[11],")*x17 + start(",randcoef[12],")*x18+start(0)*x10+start(0)*x11+start(0)*x12+start(0)*x13+start(0)*x14+start(0)*x15 \n
               f4~start(",randcoef[13],")*f1+start(0)*f2+start(0)*f3+ start(",randcoef[14],")*f5+start(0)*f6\n
               f5~start(0)*f1+start(",randcoef[15],")*f2+start(0)*f3+start(0)*f4+ start(",randcoef[16],")*f6\n
               f6~start(0)*f1+start(0)*f2+start(",randcoef[17],")*f3+ start(",randcoef[18],")*f4+start(0)*f5\n
               x1~~start(",VFAC,")*x1 \n
               x2~~start(",VFAC,")*x2 \n
               x3~~start(",VFAC,")*x3 \n
               x4~~start(",VFAC,")*x4 \n
               x5~~start(",VFAC,")*x5 \n
               x6~~start(",VFAC,")*x6 \n
               x7~~start(",VFAC,")*x7 \n
               x8~~start(",VFAC,")*x8 \n
               x9~~start(",VFAC,")*x9 \n
               x10~~start(",VFAC,")*x10 \n
               x11~~start(",VFAC,")*x11 \n
               x12~~start(",VFAC,")*x12 \n
               x13~~start(",VFAC,")*x13 \n
               x14~~start(",VFAC,")*x14 \n
               x15~~start(",VFAC,")*x15 \n
               x16~~start(",VFAC,")*x16 \n
               x17~~start(",VFAC,")*x17 \n
               x18~~start(",VFAC,")*x18 \n
               f1~~start(",VFAC,")*f1 \n
               f2~~start(",VFAC,")*f2 \n
               f3~~start(",VFAC,")*f3 \n
               f4~~start(",VFAC,")*f4 \n
               f5~~start(",VFAC,")*f5 \n
               f6~~start(",VFAC,")*f6 \n
               f1~~start(0.5*",VFAC,")*f2 \n
               f2~~start(-0.5*",VFAC,")*f3 \n
               f1~~start(0)*f3 \n")

L2[[3]]=paste0("f1=~x1+",randcoef[1],"*x2+",randcoef[2],"*x3
               f2=~x4+",randcoef[3],"*x5+",randcoef[4],"*x6
               f3=~x7 + ",randcoef[5],"*x8 + ",randcoef[6],"*x9
               f4=~x10 + ",randcoef[7],"*x11 + ",randcoef[8],"*x12
               f5=~x13 + ",randcoef[9],"*x14 + ",randcoef[10],"*x15
               f6=~x16 + ",randcoef[11],"*x17 + ",randcoef[12],"*x18
               f4~",randcoef[13],"*f1 + ",randcoef[14],"*f5
               f5~",randcoef[15],"*f2 + ",randcoef[16],"*f6
               f6~",randcoef[17],"*f3 + ",randcoef[18],"*f4
               ","x1~~",VFAC,"*x1 \n
               x2~~",VFAC,"*x2 \n
               x3~~",VFAC,"*x3 \n
               x4~~",VFAC,"*x4 \n
               x5~~",VFAC,"*x5 \n
               x6~~",VFAC,"*x6 \n
               x7~~",VFAC,"*x7 \n
               x8~~",VFAC,"*x8 \n
               x9~~",VFAC,"*x9 \n
               x10~~ ",VFAC," *x10 \n
               x11~~ ",VFAC," *x11 \n
               x12~~ ",VFAC," *x12 \n
               x13~~ ",VFAC," *x13 \n
               x14~~ ",VFAC," *x14 \n
               x15~~ ",VFAC," *x15 \n
               x16~~ ",VFAC," *x16 \n
               x17~~ ",VFAC," *x17 \n
               x18~~ ",VFAC," *x18 \n

               f1~~ ",VFAC," *f1 \n
               f2~~ ",VFAC," *f2 \n
               f3~~ ",VFAC," *f3 \n
               f4~~ ",VFAC," *f4 \n
               f5~~ ",VFAC," *f5 \n
               f6~~ ",VFAC," *f6 \n
               f1~~ 0.5*",VFAC," *f2 \n
               f2~~ -0.5 *",VFAC,"*f3 \n
               f1~~ 0 *f3 \n"
)
L3[[3]]=paste0("f1=~x1+start(0)*x2+start(0)*x3+start(0)*x4+start(0)*x5+start(0)*x6+start(0)*x7+start(0)*x8+start(0)*x9 \n
               f2=~x4+start(0)*x5+start(0)*x6 + start(0)*x1 + start(0)*x2 + start(0)*x3+start(0)*x7+start(0)*x8+start(0)*x9 \n
               f3=~x7 + start(0)*x8 + start(0)*x9 + start(0)*x1 + start(0)*x2 + start(0)*x3+start(0)*x4+start(0)*x5+start(0)*x6 \n
               f4=~x10 + start(0)*x11 + start(0)*x12+start(0)*x13+start(0)*x14+start(0)*x15+start(0)*x16+start(0)*x17+start(0)*x18 \n
               f5=~x13 + start(0)*x14 + start(0)*x15+start(0)*x10+start(0)*x11+start(0)*x12+start(0)*x16+start(0)*x17+start(0)*x18 \n
               f6=~x16 + start(0)*x17 + start(0)*x18+start(0)*x10+start(0)*x11+start(0)*x12+start(0)*x13+start(0)*x14+start(0)*x15 \n
               f4~start(0)*f1+start(0)*f2+start(0)*f3+ start(0)*f5+start(0)*f6\n
               f5~start(0)*f1+start(0)*f2+start(0)*f3+start(0)*f4+ start(0)*f6\n
               f6~start(0)*f1+start(0)*f2+start(0)*f3+ start(0)*f4+start(0)*f5\n
               x1~~start(",covcoef[1],")*x1 \n
               x2~~start(",covcoef[2],")*x2 \n
               x3~~start(",covcoef[3],")*x3 \n
               x4~~start(",covcoef[4],")*x4 \n
               x5~~start(",covcoef[5],")*x5 \n
               x6~~start(",covcoef[6],")*x6 \n
               x7~~start(",covcoef[7],")*x7 \n
               x8~~start(",covcoef[8],")*x8 \n
               x9~~start(",covcoef[9],")*x9 \n
               x10~~start(",covcoef[10],")*x10 \n
               x11~~start(",covcoef[11],")*x11 \n
               x12~~start(",covcoef[12],")*x12 \n
               x13~~start(",covcoef[13],")*x13 \n
               x14~~start(",covcoef[14],")*x14 \n
               x15~~start(",covcoef[15],")*x15 \n
               x16~~start(",covcoef[16],")*x16 \n
               x17~~start(",covcoef[17],")*x17 \n
               x18~~start(",covcoef[18],")*x18 \n
               f1~~start(",covcoef[19],")*f1 \n
               f2~~start(",covcoef[20],")*f2 \n
               f3~~start(",covcoef[21],")*f3 \n
               f4~~start(",covcoef[22],")*f4 \n
               f5~~start(",covcoef[23],")*f5 \n
               f6~~start(",covcoef[24],")*f6 \n
               f1~~start(0)*f2 \n
               f2~~start(0)*f3 \n
               f1~~start(0)*f3 \n")




L1[[4]]=paste0("f1=~x1+start(",randcoef[1],")*x2+start(",randcoef[2],")*x3+start(0)*x4+start(0)*x5+start(0)*x6+start(0)*x7+start(0)*x8+start(0)*x9 \n
               f2=~x4+start(",randcoef[3],")*x5+start(",randcoef[4],")*x6 + start(0)*x1 + start(0)*x2 + start(0)*x3+start(0)*x7+start(0)*x8+start(0)*x9 \n
               f3=~x7 + start(",randcoef[5],")*x8 + start(",randcoef[6],")*x9 + start(0)*x1 + start(0)*x2 + start(0)*x3+start(0)*x4+start(0)*x5+start(0)*x6 \n
               f4=~x10 + start(",randcoef[7],")*x11 + start(",randcoef[8],")*x12+start(0)*x13+start(0)*x14+start(0)*x15+start(0)*x16+start(0)*x17+start(0)*x18 \n
               f5=~x13 + start(",randcoef[9],")*x14 + start(",randcoef[10],")*x15+start(0)*x10+start(0)*x11+start(0)*x12+start(0)*x16+start(0)*x17+start(0)*x18 \n
               f6=~x16 + start(",randcoef[11],")*x17 + start(",randcoef[12],")*x18+start(0)*x10+start(0)*x11+start(0)*x12+start(0)*x13+start(0)*x14+start(0)*x15 \n
               f4~start(",randcoef[13],")*f1+start(0)*f2+start(0)*f3+ start(",randcoef[14],")*f5+start(0)*f6\n
               f5~start(0)*f1+start(",randcoef[15],")*f2+start(0)*f3+start(0)*f4+ start(",randcoef[16],")*f6\n
               f6~start(0)*f1+start(0)*f2+start(",randcoef[17],")*f3+ start(",randcoef[18],")*f4+start(0)*f5\n
               x1~~start(",VFAC,")*x1 \n
               x2~~start(",VFAC,")*x2 \n
               x3~~start(",VFAC,")*x3 \n
               x4~~start(",VFAC,")*x4 \n
               x5~~start(",VFAC,")*x5 \n
               x6~~start(",VFAC,")*x6 \n
               x7~~start(",VFAC,")*x7 \n
               x8~~start(",VFAC,")*x8 \n
               x9~~start(",VFAC,")*x9 \n
               x10~~start(",VFAC,")*x10 \n
               x11~~start(",VFAC,")*x11 \n
               x12~~start(",VFAC,")*x12 \n
               x13~~start(",VFAC,")*x13 \n
               x14~~start(",VFAC,")*x14 \n
               x15~~start(",VFAC,")*x15 \n
               x16~~start(",VFAC,")*x16 \n
               x17~~start(",VFAC,")*x17 \n
               x18~~start(",VFAC,")*x18 \n
               f1~~start(",VFAC,")*f1 \n
               f2~~start(",VFAC,")*f2 \n
               f3~~start(",VFAC,")*f3 \n
               f4~~start(",VFAC,")*f4 \n
               f5~~start(",VFAC,")*f5 \n
               f6~~start(",VFAC,")*f6 \n
               f1~~start(0.5*",VFAC,")*f2 \n
               f2~~start(-0.5*",VFAC,")*f3 \n
               f1~~start(0)*f3 \n")

L2[[4]]=paste0("f1=~x1+",randcoef[1],"*x2+",randcoef[2],"*x3
               f2=~x4+",randcoef[3],"*x5+",randcoef[4],"*x6
               f3=~x7 + ",randcoef[5],"*x8 + ",randcoef[6],"*x9
               f4=~x10 + ",randcoef[7],"*x11 + ",randcoef[8],"*x12
               f5=~x13 + ",randcoef[9],"*x14 + ",randcoef[10],"*x15
               f6=~x16 + ",randcoef[11],"*x17 + ",randcoef[12],"*x18
               f4~",randcoef[13],"*f1 + ",randcoef[14],"*f5
               f5~",randcoef[15],"*f2 + ",randcoef[16],"*f6
               f6~",randcoef[17],"*f3 + ",randcoef[18],"*f4
               ","x1~~",VFAC,"*x1 \n
               x2~~",VFAC,"*x2 \n
               x3~~",VFAC,"*x3 \n
               x4~~",VFAC,"*x4 \n
               x5~~",VFAC,"*x5 \n
               x6~~",VFAC,"*x6 \n
               x7~~",VFAC,"*x7 \n
               x8~~",VFAC,"*x8 \n
               x9~~",VFAC,"*x9 \n
               x10~~ ",VFAC," *x10 \n
               x11~~ ",VFAC," *x11 \n
               x12~~ ",VFAC," *x12 \n
               x13~~ ",VFAC," *x13 \n
               x14~~ ",VFAC," *x14 \n
               x15~~ ",VFAC," *x15 \n
               x16~~ ",VFAC," *x16 \n
               x17~~ ",VFAC," *x17 \n
               x18~~ ",VFAC," *x18 \n

               f1~~ ",VFAC," *f1 \n
               f2~~ ",VFAC," *f2 \n
               f3~~ ",VFAC," *f3 \n
               f4~~ ",VFAC," *f4 \n
               f5~~ ",VFAC," *f5 \n
               f6~~ ",VFAC," *f6 \n
               f1~~ 0.5 *",VFAC,"*f2 \n
               f2~~ -0.5 *",VFAC,"*f3 \n
               f1~~ 0 *f3 \n"
)
L3[[4]]=paste0("f1=~x1+start(",initcoef[1],")*x2+start(",initcoef[2],")*x3+start(0)*x4+start(0)*x5+start(0)*x6+start(0)*x7+start(0)*x8+start(0)*x9 \n
               f2=~x4+start(",initcoef[3],")*x5+start(",initcoef[4],")*x6 + start(0)*x1 + start(0)*x2 + start(0)*x3+start(0)*x7+start(0)*x8+start(0)*x9 \n
               f3=~x7 + start(",initcoef[5],")*x8 + start(",initcoef[6],")*x9 + start(0)*x1 + start(0)*x2 + start(0)*x3+start(0)*x4+start(0)*x5+start(0)*x6 \n
               f4=~x10 + start(",initcoef[7],")*x11 + start(",initcoef[8],")*x12+start(0)*x13+start(0)*x14+start(0)*x15+start(0)*x16+start(0)*x17+start(0)*x18 \n
               f5=~x13 + start(",initcoef[9],")*x14 + start(",initcoef[10],")*x15+start(0)*x10+start(0)*x11+start(0)*x12+start(0)*x16+start(0)*x17+start(0)*x18 \n
               f6=~x16 + start(",initcoef[11],")*x17 + start(",initcoef[12],")*x18+start(0)*x10+start(0)*x11+start(0)*x12+start(0)*x13+start(0)*x14+start(0)*x15 \n
               f4~start(0)*f1+start(0)*f2+start(0)*f3+ start(",initcoef[13],")*f5+start(0)*f6\n
               f5~start(0)*f1+start(",initcoef[14],")*f2+start(0)*f3+start(0)*f4+ start(0)*f6\n
               f6~start(0)*f1+start(0)*f2+start(0)*f3+ start(0)*f4+start(0)*f5\n
               x1~~start(",covcoef[1],")*x1 \n
               x2~~start(",covcoef[2],")*x2 \n
               x3~~start(",covcoef[3],")*x3 \n
               x4~~start(",covcoef[4],")*x4 \n
               x5~~start(",covcoef[5],")*x5 \n
               x6~~start(",covcoef[6],")*x6 \n
               x7~~start(",covcoef[7],")*x7 \n
               x8~~start(",covcoef[8],")*x8 \n
               x9~~start(",covcoef[9],")*x9 \n
               x10~~start(",covcoef[10],")*x10 \n
               x11~~start(",covcoef[11],")*x11 \n
               x12~~start(",covcoef[12],")*x12 \n
               x13~~start(",covcoef[13],")*x13 \n
               x14~~start(",covcoef[14],")*x14 \n
               x15~~start(",covcoef[15],")*x15 \n
               x16~~start(",covcoef[16],")*x16 \n
               x17~~start(",covcoef[17],")*x17 \n
               x18~~start(",covcoef[18],")*x18 \n
               f1~~start(",covcoef[19],")*f1 \n
               f2~~start(",covcoef[20],")*f2 \n
               f3~~start(",covcoef[21],")*f3 \n
               f4~~start(",covcoef[22],")*f4 \n
               f5~~start(",covcoef[23],")*f5 \n
               f6~~start(",covcoef[24],")*f6 \n
               f1~~start(0)*f2 \n
               f2~~start(0)*f3 \n
               f1~~start(0)*f3 \n")




L1[[7]]=paste0("f1=~x1+start(",randcoef[1],")*x2+start(",randcoef[2],")*x3+start(0)*x4+start(0)*x5+start(0)*x6+start(0)*x7+start(0)*x8+start(0)*x9 \n
               f2=~x4+start(",randcoef[3],")*x5+start(",randcoef[4],")*x6 + start(0)*x1 + start(0)*x2 + start(0)*x3+start(0)*x7+start(0)*x8+start(0)*x9 \n
               f3=~x7 + start(",randcoef[5],")*x8 + start(",randcoef[6],")*x9 + start(0)*x1 + start(0)*x2 + start(0)*x3+start(0)*x4+start(0)*x5+start(0)*x6 \n
               f4=~x10 + start(",randcoef[7],")*x11 + start(",randcoef[8],")*x12+start(0)*x13+start(0)*x14+start(0)*x15+start(0)*x16+start(0)*x17+start(0)*x18 \n
               f5=~x13 + start(",randcoef[9],")*x14 + start(",randcoef[10],")*x15+start(0)*x10+start(0)*x11+start(0)*x12+start(0)*x16+start(0)*x17+start(0)*x18 \n
               f6=~x16 + start(",randcoef[11],")*x17 + start(",randcoef[12],")*x18+start(0)*x10+start(0)*x11+start(0)*x12+start(0)*x13+start(0)*x14+start(0)*x15 \n
               f4~start(",randcoef[13],")*f1+start(0)*f2+start(0)*f3+ start(0)*f5+start(0)*f6\n
               f5~start(0)*f1+start(",randcoef[14],")*f2+start(0)*f3+start(0)*f4+ start(0)*f6\n
               f6~start(0)*f1+start(0)*f2+start(",randcoef[15],")*f3+ start(0)*f4+start(0)*f5\n
               x1~~start(",VFAC,")*x1 \n
               x2~~start(",VFAC,")*x2 \n
               x3~~start(",VFAC,")*x3 \n
               x4~~start(",VFAC,")*x4 \n
               x5~~start(",VFAC,")*x5 \n
               x6~~start(",VFAC,")*x6 \n
               x7~~start(",VFAC,")*x7 \n
               x8~~start(",VFAC,")*x8 \n
               x9~~start(",VFAC,")*x9 \n
               x10~~start(",VFAC,")*x10 \n
               x11~~start(",VFAC,")*x11 \n
               x12~~start(",VFAC,")*x12 \n
               x13~~start(",VFAC,")*x13 \n
               x14~~start(",VFAC,")*x14 \n
               x15~~start(",VFAC,")*x15 \n
               x16~~start(",VFAC,")*x16 \n
               x17~~start(",VFAC,")*x17 \n
               x18~~start(",VFAC,")*x18 \n
               f1~~start(",VFAC,")*f1 \n
               f2~~start(",VFAC,")*f2 \n
               f3~~start(",VFAC,")*f3 \n
               f4~~start(",VFAC,")*f4 \n
               f5~~start(",VFAC,")*f5 \n
               f6~~start(",VFAC,")*f6 \n
               f1~~start(0.5*",VFAC,")*f2 \n
               f2~~start(-0.5*",VFAC,")*f3 \n
               f1~~start(0)*f3 \n")

L2[[7]]=paste0("f1=~x1+",randcoef[1],"*x2+",randcoef[2],"*x3
               f2=~x4+",randcoef[3],"*x5+",randcoef[4],"*x6
               f3=~x7 + ",randcoef[5],"*x8 + ",randcoef[6],"*x9
               f4=~x10 + ",randcoef[7],"*x11 + ",randcoef[8],"*x12
               f5=~x13 + ",randcoef[9],"*x14 + ",randcoef[10],"*x15
               f6=~x16 + ",randcoef[11],"*x17 + ",randcoef[12],"*x18
               f4~",randcoef[13],"*f1
               f5~",randcoef[14],"*f2
               f6~",randcoef[15],"*f3
               ","x1~~",VFAC,"*x1 \n
               x2~~",VFAC,"*x2 \n
               x3~~",VFAC,"*x3 \n
               x4~~",VFAC,"*x4 \n
               x5~~",VFAC,"*x5 \n
               x6~~",VFAC,"*x6 \n
               x7~~",VFAC,"*x7 \n
               x8~~",VFAC,"*x8 \n
               x9~~",VFAC,"*x9 \n
               x10~~ ",VFAC," *x10 \n
               x11~~ ",VFAC," *x11 \n
               x12~~ ",VFAC," *x12 \n
               x13~~ ",VFAC," *x13 \n
               x14~~ ",VFAC," *x14 \n
               x15~~ ",VFAC," *x15 \n
               x16~~ ",VFAC," *x16 \n
               x17~~ ",VFAC," *x17 \n
               x18~~ ",VFAC," *x18 \n

               f1~~ 1 *f1 \n
               f2~~ 1 *f2 \n
               f3~~ 1 *f3 \n
               f4~~ 1 *f4 \n
               f5~~ 1 *f5 \n
               f6~~ 1 *f6 \n
               f1~~ 0.5 *",VFAC,"*f2 \n
               f2~~ -0.5 *",VFAC,"*f3 \n
               f1~~ 0 *f3 \n"
)
L3[[7]]=paste0("f1=~x1+start(0)*x2+start(0)*x3+start(0)*x4+start(0)*x5+start(0)*x6+start(0)*x7+start(0)*x8+start(0)*x9 \n
               f2=~x4+start(0)*x5+start(0)*x6 + start(0)*x1 + start(0)*x2 + start(0)*x3+start(0)*x7+start(0)*x8+start(0)*x9 \n
               f3=~x7 + start(0)*x8 + start(0)*x9 + start(0)*x1 + start(0)*x2 + start(0)*x3+start(0)*x4+start(0)*x5+start(0)*x6 \n
               f4=~x10 + start(0)*x11 + start(0)*x12+start(0)*x13+start(0)*x14+start(0)*x15+start(0)*x16+start(0)*x17+start(0)*x18 \n
               f5=~x13 + start(0)*x14 + start(0)*x15+start(0)*x10+start(0)*x11+start(0)*x12+start(0)*x16+start(0)*x17+start(0)*x18 \n
               f6=~x16 + start(0)*x17 + start(0)*x18+start(0)*x10+start(0)*x11+start(0)*x12+start(0)*x13+start(0)*x14+start(0)*x15 \n
               f4~start(0)*f1+start(0)*f2+start(0)*f3+ start(0)*f5+start(0)*f6\n
               f5~start(0)*f1+start(0)*f2+start(0)*f3+start(0)*f4+ start(0)*f6\n
               f6~start(0)*f1+start(0)*f2+start(0)*f3+ start(0)*f4+start(0)*f5\n
               x1~~start(",covcoef[1],")*x1 \n
               x2~~start(",covcoef[2],")*x2 \n
               x3~~start(",covcoef[3],")*x3 \n
               x4~~start(",covcoef[4],")*x4 \n
               x5~~start(",covcoef[5],")*x5 \n
               x6~~start(",covcoef[6],")*x6 \n
               x7~~start(",covcoef[7],")*x7 \n
               x8~~start(",covcoef[8],")*x8 \n
               x9~~start(",covcoef[9],")*x9 \n
               x10~~start(",covcoef[10],")*x10 \n
               x11~~start(",covcoef[11],")*x11 \n
               x12~~start(",covcoef[12],")*x12 \n
               x13~~start(",covcoef[13],")*x13 \n
               x14~~start(",covcoef[14],")*x14 \n
               x15~~start(",covcoef[15],")*x15 \n
               x16~~start(",covcoef[16],")*x16 \n
               x17~~start(",covcoef[17],")*x17 \n
               x18~~start(",covcoef[18],")*x18 \n
               f1~~start(",covcoef[19],")*f1 \n
               f2~~start(",covcoef[20],")*f2 \n
               f3~~start(",covcoef[21],")*f3 \n
               f4~~start(",covcoef[22],")*f4 \n
               f5~~start(",covcoef[23],")*f5 \n
               f6~~start(",covcoef[24],")*f6 \n
               f1~~start(0)*f2 \n
               f2~~start(0)*f3 \n
               f1~~start(0)*f3 \n")

L1[[8]]=paste0("f1=~x1+start(",randcoef[1],")*x2+start(",randcoef[2],")*x3+start(0)*x4+start(0)*x5+start(0)*x6+start(0)*x7+start(0)*x8+start(0)*x9 \n
               f2=~x4+start(",randcoef[3],")*x5+start(",randcoef[4],")*x6 + start(0)*x1 + start(0)*x2 + start(0)*x3+start(0)*x7+start(0)*x8+start(0)*x9 \n
               f3=~x7 + start(",randcoef[5],")*x8 + start(",randcoef[6],")*x9 + start(0)*x1 + start(0)*x2 + start(0)*x3+start(0)*x4+start(0)*x5+start(0)*x6 \n
               f4=~x10 + start(",randcoef[7],")*x11 + start(",randcoef[8],")*x12+start(0)*x13+start(0)*x14+start(0)*x15+start(0)*x16+start(0)*x17+start(0)*x18 \n
               f5=~x13 + start(",randcoef[9],")*x14 + start(",randcoef[10],")*x15+start(0)*x10+start(0)*x11+start(0)*x12+start(0)*x16+start(0)*x17+start(0)*x18 \n
               f6=~x16 + start(",randcoef[11],")*x17 + start(",randcoef[12],")*x18+start(0)*x10+start(0)*x11+start(0)*x12+start(0)*x13+start(0)*x14+start(0)*x15 \n
               f4~start(",randcoef[13],")*f1+start(0)*f2+start(0)*f3+ start(0)*f5+start(0)*f6\n
               f5~start(0)*f1+start(",randcoef[14],")*f2+start(0)*f3+start(0)*f4+ start(0)*f6\n
               f6~start(0)*f1+start(0)*f2+start(",randcoef[15],")*f3+ start(0)*f4+start(0)*f5\n
               x1~~start(",VFAC,")*x1 \n
               x2~~start(",VFAC,")*x2 \n
               x3~~start(",VFAC,")*x3 \n
               x4~~start(",VFAC,")*x4 \n
               x5~~start(",VFAC,")*x5 \n
               x6~~start(",VFAC,")*x6 \n
               x7~~start(",VFAC,")*x7 \n
               x8~~start(",VFAC,")*x8 \n
               x9~~start(",VFAC,")*x9 \n
               x10~~start(",VFAC,")*x10 \n
               x11~~start(",VFAC,")*x11 \n
               x12~~start(",VFAC,")*x12 \n
               x13~~start(",VFAC,")*x13 \n
               x14~~start(",VFAC,")*x14 \n
               x15~~start(",VFAC,")*x15 \n
               x16~~start(",VFAC,")*x16 \n
               x17~~start(",VFAC,")*x17 \n
               x18~~start(",VFAC,")*x18 \n
               f1~~start(",VFAC,")*f1 \n
               f2~~start(",VFAC,")*f2 \n
               f3~~start(",VFAC,")*f3 \n
               f4~~start(",VFAC,")*f4 \n
               f5~~start(",VFAC,")*f5 \n
               f6~~start(",VFAC,")*f6 \n
               f1~~start(0.5*",VFAC,")*f2 \n
               f2~~start(-0.5*",VFAC,")*f3 \n
               f1~~start(0)*f3 \n")

L2[[8]]=paste0("f1=~x1+",randcoef[1],"*x2+",randcoef[2],"*x3
               f2=~x4+",randcoef[3],"*x5+",randcoef[4],"*x6
               f3=~x7 + ",randcoef[5],"*x8 + ",randcoef[6],"*x9
               f4=~x10 + ",randcoef[7],"*x11 + ",randcoef[8],"*x12
               f5=~x13 + ",randcoef[9],"*x14 + ",randcoef[10],"*x15
               f6=~x16 + ",randcoef[11],"*x17 + ",randcoef[12],"*x18
               f4~",randcoef[13],"*f1
               f5~",randcoef[14],"*f2
               f6~",randcoef[15],"*f3
               ","x1~~",VFAC,"*x1 \n
               x2~~",VFAC,"*x2 \n
               x3~~",VFAC,"*x3 \n
               x4~~",VFAC,"*x4 \n
               x5~~",VFAC,"*x5 \n
               x6~~",VFAC,"*x6 \n
               x7~~",VFAC,"*x7 \n
               x8~~",VFAC,"*x8 \n
               x9~~",VFAC,"*x9 \n
               x10~~ ",VFAC," *x10 \n
               x11~~ ",VFAC," *x11 \n
               x12~~ ",VFAC," *x12 \n
               x13~~ ",VFAC," *x13 \n
               x14~~ ",VFAC," *x14 \n
               x15~~ ",VFAC," *x15 \n
               x16~~ ",VFAC," *x16 \n
               x17~~ ",VFAC," *x17 \n
               x18~~ ",VFAC," *x18 \n

               f1~~ ",VFAC," *f1 \n
               f2~~ ",VFAC," *f2 \n
               f3~~ ",VFAC," *f3 \n
               f4~~ ",VFAC," *f4 \n
               f5~~ ",VFAC," *f5 \n
               f6~~ ",VFAC," *f6 \n
               f1~~ 0.5*",VFAC," *f2 \n
               f2~~ -0.5*",VFAC," *f3 \n
               f1~~ 0 *f3 \n"
)
L3[[8]]=paste0("f1=~x1+start(",initcoef[1],")*x2+start(",initcoef[2],")*x3+start(0)*x4+start(0)*x5+start(0)*x6+start(0)*x7+start(0)*x8+start(0)*x9 \n
               f2=~x4+start(",initcoef[3],")*x5+start(",initcoef[4],")*x6 + start(0)*x1 + start(0)*x2 + start(0)*x3+start(0)*x7+start(0)*x8+start(0)*x9 \n
               f3=~x7 + start(",initcoef[5],")*x8 + start(",initcoef[6],")*x9 + start(0)*x1 + start(0)*x2 + start(0)*x3+start(0)*x4+start(0)*x5+start(0)*x6 \n
               f4=~x10 + start(",initcoef[7],")*x11 + start(",initcoef[8],")*x12+start(0)*x13+start(0)*x14+start(0)*x15+start(0)*x16+start(0)*x17+start(0)*x18 \n
               f5=~x13 + start(",initcoef[9],")*x14 + start(",initcoef[10],")*x15+start(0)*x10+start(0)*x11+start(0)*x12+start(0)*x16+start(0)*x17+start(0)*x18 \n
               f6=~x16 + start(",initcoef[11],")*x17 + start(",initcoef[12],")*x18+start(0)*x10+start(0)*x11+start(0)*x12+start(0)*x13+start(0)*x14+start(0)*x15 \n
               f4~start(0)*f1+start(0)*f2+start(0)*f3+ start(0)*f5+start(0)*f6\n
               f5~start(0)*f1+start(",initcoef[13],")*f2+start(0)*f3+start(0)*f4+ start(0)*f6\n
               f6~start(0)*f1+start(0)*f2+start(0)*f3+ start(0)*f4+start(0)*f5\n
               x1~~start(",covcoef[1],")*x1 \n
               x2~~start(",covcoef[2],")*x2 \n
               x3~~start(",covcoef[3],")*x3 \n
               x4~~start(",covcoef[4],")*x4 \n
               x5~~start(",covcoef[5],")*x5 \n
               x6~~start(",covcoef[6],")*x6 \n
               x7~~start(",covcoef[7],")*x7 \n
               x8~~start(",covcoef[8],")*x8 \n
               x9~~start(",covcoef[9],")*x9 \n
               x10~~start(",covcoef[10],")*x10 \n
               x11~~start(",covcoef[11],")*x11 \n
               x12~~start(",covcoef[12],")*x12 \n
               x13~~start(",covcoef[13],")*x13 \n
               x14~~start(",covcoef[14],")*x14 \n
               x15~~start(",covcoef[15],")*x15 \n
               x16~~start(",covcoef[16],")*x16 \n
               x17~~start(",covcoef[17],")*x17 \n
               x18~~start(",covcoef[18],")*x18 \n
               f1~~start(",covcoef[19],")*f1 \n
               f2~~start(",covcoef[20],")*f2 \n
               f3~~start(",covcoef[21],")*f3 \n
               f4~~start(",covcoef[22],")*f4 \n
               f5~~start(",covcoef[23],")*f5 \n
               f6~~start(",covcoef[24],")*f6 \n
               f1~~start(0)*f2 \n
               f2~~start(0)*f3 \n
               f1~~start(0)*f3 \n")



L1[[9]]=paste0("f1=~x1+start(",randcoef[1],")*x2+start(",randcoef[2],")*x3+start(0)*x4+start(0)*x5+start(0)*x6+start(0)*x7+start(0)*x8+start(0)*x9 \n
               f2=~x4+start(",randcoef[3],")*x5+start(",randcoef[4],")*x6 + start(0)*x1 + start(0)*x2 + start(0)*x3+start(0)*x7+start(0)*x8+start(0)*x9 \n
               f3=~x7 + start(",randcoef[5],")*x8 + start(",randcoef[6],")*x9 + start(0)*x1 + start(0)*x2 + start(0)*x3+start(0)*x4+start(0)*x5+start(0)*x6 \n
               f4=~x10 + start(",randcoef[7],")*x11 + start(",randcoef[8],")*x12+start(0)*x13+start(0)*x14+start(0)*x15+start(0)*x16+start(0)*x17+start(0)*x18 \n
               f5=~x13 + start(",randcoef[9],")*x14 + start(",randcoef[10],")*x15+start(0)*x10+start(0)*x11+start(0)*x12+start(0)*x16+start(0)*x17+start(0)*x18 \n
               f6=~x16 + start(",randcoef[11],")*x17 + start(",randcoef[12],")*x18+start(0)*x10+start(0)*x11+start(0)*x12+start(0)*x13+start(0)*x14+start(0)*x15 \n
               f4~start(",randcoef[13],")*f1+start(0)*f2+start(0)*f3+ start(",randcoef[14],")*f5+start(0)*f6
               + start(",randcoef[19],")*x19 + start(0)*x20+ start(0)*x21\n
               f5~start(0)*f1+start(",randcoef[15],")*f2+start(0)*f3+start(0)*f4+ start(",randcoef[16],")*f6
               + start(0)*x19 + start(",randcoef[20],")*x20+ start(0)*x21\n
               f6~start(0)*f1+start(0)*f2+start(",randcoef[17],")*f3+ start(",randcoef[18],")*f4+start(0)*f5
               + start(0)*x19 + start(0)*x20+ start(",randcoef[21],")*x21\n
               x1~~start(",VFAC,")*x1 \n
               x2~~start(",VFAC,")*x2 \n
               x3~~start(",VFAC,")*x3 \n
               x4~~start(",VFAC,")*x4 \n
               x5~~start(",VFAC,")*x5 \n
               x6~~start(",VFAC,")*x6 \n
               x7~~start(",VFAC,")*x7 \n
               x8~~start(",VFAC,")*x8 \n
               x9~~start(",VFAC,")*x9 \n
               x10~~start(",VFAC,")*x10 \n
               x11~~start(",VFAC,")*x11 \n
               x12~~start(",VFAC,")*x12 \n
               x13~~start(",VFAC,")*x13 \n
               x14~~start(",VFAC,")*x14 \n
               x15~~start(",VFAC,")*x15 \n
               x16~~start(",VFAC,")*x16 \n
               x17~~start(",VFAC,")*x17 \n
               x18~~start(",VFAC,")*x18 \n

               f1~~start(",VFAC,")*f1 \n
               f2~~start(",VFAC,")*f2 \n
               f3~~start(",VFAC,")*f3 \n
               f4~~start(",VFAC,")*f4 \n
               f5~~start(",VFAC,")*f5 \n
               f6~~start(",VFAC,")*f6 \n
               f1~~start(0.5*",VFAC,")*f2 \n
               f2~~start(-0.5*",VFAC,")*f3 \n
               f1~~start(0)*f3 \n
               x19~~start(",VFAC,")*x19 \n
               x20~~start(",VFAC,")*x20 \n
               x21~~start(",VFAC,")*x21 \n
               x19~~start(0) * x20 \n
               x20~~start(0)*x21 \n
               x19~~start(0)*x21 \n")

L2[[9]]=paste0("f1=~x1+",randcoef[1],"*x2+",randcoef[2],"*x3
               f2=~x4+",randcoef[3],"*x5+",randcoef[4],"*x6
               f3=~x7 + ",randcoef[5],"*x8 + ",randcoef[6],"*x9
               f4=~x10 + ",randcoef[7],"*x11 + ",randcoef[8],"*x12
               f5=~x13 + ",randcoef[9],"*x14 + ",randcoef[10],"*x15
               f6=~x16 + ",randcoef[11],"*x17 + ",randcoef[12],"*x18
               f4~",randcoef[13],"*f1 + ",randcoef[14],"*f5 + ",randcoef[19],"*x19
               f5~",randcoef[15],"*f2 + ",randcoef[16],"*f6 + ",randcoef[20],"*x20
               f6~",randcoef[17],"*f3 + ",randcoef[18],"*f4 + ",randcoef[21],"*x21
               ","x1~~",VFAC,"*x1 \n
               x2~~",VFAC,"*x2 \n
               x3~~",VFAC,"*x3 \n
               x4~~",VFAC,"*x4 \n
               x5~~",VFAC,"*x5 \n
               x6~~",VFAC,"*x6 \n
               x7~~",VFAC,"*x7 \n
               x8~~",VFAC,"*x8 \n
               x9~~",VFAC,"*x9 \n
               x10~~ ",VFAC," *x10 \n
               x11~~ ",VFAC," *x11 \n
               x12~~ ",VFAC," *x12 \n
               x13~~ ",VFAC," *x13 \n
               x14~~ ",VFAC," *x14 \n
               x15~~ ",VFAC," *x15 \n
               x16~~ ",VFAC," *x16 \n
               x17~~ ",VFAC," *x17 \n
               x18~~ ",VFAC," *x18 \n
               x19~~ ",VFAC," *x19 \n
               x20~~ ",VFAC," *x20 \n
               x21~~ ",VFAC," *x21 \n

               f1~~ ",VFAC," *f1 \n
               f2~~ ",VFAC," *f2 \n
               f3~~ ",VFAC," *f3 \n
               f4~~ ",VFAC," *f4 \n
               f5~~ ",VFAC," *f5 \n
               f6~~ ",VFAC," *f6 \n
               f1~~ 0.5 *",VFAC,"*f2 \n
               f2~~ -0.5 *",VFAC,"*f3 \n
               f1~~ 0 *f3 \n

               x19~~0 * x20 \n
               x20~~0*x21 \n
               x19~~0*x21 \n"
)
L3[[9]]=paste0(
  "f1=~x1+start(",initcoef[1],")*x2+start(",initcoef[2],")*x3+start(0)*x4+start(0)*x5+start(0)*x6+start(0)*x7+start(0)*x8+start(0)*x9 \n
  f2=~x4+start(",initcoef[3],")*x5+start(",initcoef[4],")*x6 + start(0)*x1 + start(0)*x2 + start(0)*x3+start(0)*x7+start(0)*x8+start(0)*x9 \n
  f3=~x7 + start(",initcoef[5],")*x8 + start(",initcoef[6],")*x9 + start(0)*x1 + start(0)*x2 + start(0)*x3+start(0)*x4+start(0)*x5+start(0)*x6 \n
  f4=~x10 + start(",initcoef[7],")*x11 + start(",initcoef[8],")*x12+start(0)*x13+start(0)*x14+start(0)*x15+start(0)*x16+start(0)*x17+start(0)*x18 \n
  f5=~x13 + start(",initcoef[9],")*x14 + start(",initcoef[10],")*x15+start(0)*x10+start(0)*x11+start(0)*x12+start(0)*x16+start(0)*x17+start(0)*x18 \n
  f6=~x16 + start(",initcoef[11],")*x17 + start(",initcoef[12],")*x18+start(0)*x10+start(0)*x11+start(0)*x12+start(0)*x13+start(0)*x14+start(0)*x15 \n
  f4~start(0)*f1+start(0)*f2+start(0)*f3+ start(",initcoef[12],")*f5+start(0)*f6
  + start(0)*x19 + start(0)*x20+ start(0)*x21\n
  f5~start(0)*f1+start(",initcoef[13],")*f2+start(0)*f3+start(0)*f4+ start(0)*f6
  + start(0)*x19 + start(0)*x20+ start(0)*x21\n
  f6~start(0)*f1+start(0)*f2+start(0)*f3+ start(0)*f4+start(0)*f5
  + start(0)*x19 + start(0)*x20+ start(0)*x21\n
  x1~~start(",covcoef[1],")*x1 \n
  x2~~start(",covcoef[2],")*x2 \n
  x3~~start(",covcoef[3],")*x3 \n
  x4~~start(",covcoef[4],")*x4 \n
  x5~~start(",covcoef[5],")*x5 \n
  x6~~start(",covcoef[6],")*x6 \n
  x7~~start(",covcoef[7],")*x7 \n
  x8~~start(",covcoef[8],")*x8 \n
  x9~~start(",covcoef[9],")*x9 \n
  x10~~start(",covcoef[10],")*x10 \n
  x11~~start(",covcoef[11],")*x11 \n
  x12~~start(",covcoef[12],")*x12 \n
  x13~~start(",covcoef[13],")*x13 \n
  x14~~start(",covcoef[14],")*x14 \n
  x15~~start(",covcoef[15],")*x15 \n
  x16~~start(",covcoef[16],")*x16 \n
  x17~~start(",covcoef[17],")*x17 \n
  x18~~start(",covcoef[18],")*x18 \n
  f1~~start(",covcoef[19],")*f1 \n
  f2~~start(",covcoef[20],")*f2 \n
  f3~~start(",covcoef[21],")*f3 \n
  f4~~start(",covcoef[22],")*f4 \n
  f5~~start(",covcoef[23],")*f5 \n
  f6~~start(",covcoef[24],")*f6 \n
  f1~~start(0)*f2 \n
  f2~~start(0)*f3 \n
  f1~~start(0)*f3 \n
  x19~~start(",covcoef[25],")*x19 \n
  x20~~start(",covcoef[26],")*x20 \n
  x21~~start(",covcoef[27],")*x21 \n
  x19~~start(0)*x20 \n
  x20~~start(0)*x21 \n
  x19~~start(0)*x21 \n")







L1[[10]]=paste0("f1=~x1+start(",randcoef[1],")*x2+start(",randcoef[2],")*x3+start(0)*x4+start(0)*x5+start(0)*x6+start(0)*x7+start(0)*x8+start(0)*x9 \n
                f2=~x4+start(",randcoef[3],")*x5+start(",randcoef[4],")*x6 + start(0)*x1 + start(0)*x2 + start(0)*x3+start(0)*x7+start(0)*x8+start(0)*x9 \n
                f3=~x7 + start(",randcoef[5],")*x8 + start(",randcoef[6],")*x9 + start(0)*x1 + start(0)*x2 + start(0)*x3+start(0)*x4+start(0)*x5+start(0)*x6 \n
                f4=~x10 + start(",randcoef[7],")*x11 + start(",randcoef[8],")*x12+start(0)*x13+start(0)*x14+start(0)*x15+start(0)*x16+start(0)*x17+start(0)*x18 \n
                f5=~x13 + start(",randcoef[9],")*x14 + start(",randcoef[10],")*x15+start(0)*x10+start(0)*x11+start(0)*x12+start(0)*x16+start(0)*x17+start(0)*x18 \n
                f6=~x16 + start(",randcoef[11],")*x17 + start(",randcoef[12],")*x18+start(0)*x10+start(0)*x11+start(0)*x12+start(0)*x13+start(0)*x14+start(0)*x15 \n
                f4~start(",randcoef[13],")*f1+start(0)*f2+start(0)*f3+ start(",randcoef[14],")*f5+start(0)*f6
                + start(",randcoef[19],")*x19 + start(0)*x20+ start(0)*x21\n
                f5~start(0)*f1+start(",randcoef[15],")*f2+start(0)*f3+start(0)*f4+ start(",randcoef[16],")*f6
                + start(0)*x19 + start(",randcoef[20],")*x20+ start(0)*x21\n
                f6~start(0)*f1+start(0)*f2+start(",randcoef[17],")*f3+ start(",randcoef[18],")*f4+start(0)*f5
                + start(0)*x19 + start(0)*x20+ start(",randcoef[21],")*x21\n
                x1~~start(",VFAC,")*x1 \n
                x2~~start(",VFAC,")*x2 \n
                x3~~start(",VFAC,")*x3 \n
                x4~~start(",VFAC,")*x4 \n
                x5~~start(",VFAC,")*x5 \n
                x6~~start(",VFAC,")*x6 \n
                x7~~start(",VFAC,")*x7 \n
                x8~~start(",VFAC,")*x8 \n
                x9~~start(",VFAC,")*x9 \n
                x10~~start(",VFAC,")*x10 \n
                x11~~start(",VFAC,")*x11 \n
                x12~~start(",VFAC,")*x12 \n
                x13~~start(",VFAC,")*x13 \n
                x14~~start(",VFAC,")*x14 \n
                x15~~start(",VFAC,")*x15 \n
                x16~~start(",VFAC,")*x16 \n
                x17~~start(",VFAC,")*x17 \n
                x18~~start(",VFAC,")*x18 \n
                x19~~start(",VFAC,")*x19 \n
                x20~~start(",VFAC,")*x20 \n
                x21~~start(",VFAC,")*x21 \n

                f1~~start(",VFAC,")*f1 \n
                f2~~start(",VFAC,")*f2 \n
                f3~~start(",VFAC,")*f3 \n
                f4~~start(",VFAC,")*f4 \n
                f5~~start(",VFAC,")*f5 \n
                f6~~start(",VFAC,")*f6 \n
                f1~~start(0.5*",VFAC,")*f2 \n
                f2~~start(-0.5*",VFAC,")*f3 \n
                f1~~start(0)*f3 \n

                x19~~start(0) * x20 \n
                x20~~start(0) * x21 \n
                x19~~start(0) * x21 \n")

L2[[10]]=paste0("f1=~x1+",randcoef[1],"*x2+",randcoef[2],"*x3
                f2=~x4+",randcoef[3],"*x5+",randcoef[4],"*x6
                f3=~x7 + ",randcoef[5],"*x8 + ",randcoef[6],"*x9
                f4=~x10 + ",randcoef[7],"*x11 + ",randcoef[8],"*x12
                f5=~x13 + ",randcoef[9],"*x14 + ",randcoef[10],"*x15
                f6=~x16 + ",randcoef[11],"*x17 + ",randcoef[12],"*x18
                f4~",randcoef[13],"*f1 + ",randcoef[14],"*f5 + ",randcoef[19],"*x19
                f5~",randcoef[15],"*f2 + ",randcoef[16],"*f6 + ",randcoef[20],"*x20
                f6~",randcoef[17],"*f3 + ",randcoef[18],"*f4 + ",randcoef[21],"*x21
                ","x1~~",VFAC,"*x1 \n
                x2~~",VFAC,"*x2 \n
                x3~~",VFAC,"*x3 \n
                x4~~",VFAC,"*x4 \n
                x5~~",VFAC,"*x5 \n
                x6~~",VFAC,"*x6 \n
                x7~~",VFAC,"*x7 \n
                x8~~",VFAC,"*x8 \n
                x9~~",VFAC,"*x9 \n
                x10~~ ",VFAC," *x10 \n
                x11~~ ",VFAC," *x11 \n
                x12~~ ",VFAC," *x12 \n
                x13~~ ",VFAC," *x13 \n
                x14~~ ",VFAC," *x14 \n
                x15~~ ",VFAC," *x15 \n
                x16~~ ",VFAC," *x16 \n
                x17~~ ",VFAC," *x17 \n
                x18~~ ",VFAC," *x18 \n
                x19~~ ",VFAC," *x19 \n
                x20~~ ",VFAC," *x20 \n
                x21~~ ",VFAC," *x21 \n

                f1~~ ",VFAC," *f1 \n
                f2~~ ",VFAC," *f2 \n
                f3~~ ",VFAC," *f3 \n
                f4~~ ",VFAC," *f4 \n
                f5~~ ",VFAC," *f5 \n
                f6~~ ",VFAC," *f6 \n
                f1~~ 0.5  *",VFAC,"*f2 \n
                f2~~ -0.5 *",VFAC,"*f3 \n
                f1~~ 0 *f3 \n

                x19~~0 * x20 \n
                x20~~0*x21 \n
                x19~~0*x21 \n"
)
L3[[10]]=paste0(
  "f1=~x1+start(0)*x2+start(0)*x3+start(0)*x4+start(0)*x5+start(0)*x6+start(0)*x7+start(0)*x8+start(0)*x9 \n
  f2=~x4+start(0)*x5+start(0)*x6 + start(0)*x1 + start(0)*x2 + start(0)*x3+start(0)*x7+start(0)*x8+start(0)*x9 \n
  f3=~x7 + start(0)*x8 + start(0)*x9 + start(0)*x1 + start(0)*x2 + start(0)*x3+start(0)*x4+start(0)*x5+start(0)*x6 \n
  f4=~x10 + start(0)*x11 + start(0)*x12+start(0)*x13+start(0)*x14+start(0)*x15+start(0)*x16+start(0)*x17+start(0)*x18 \n
  f5=~x13 + start(0)*x14 + start(0)*x15+start(0)*x10+start(0)*x11+start(0)*x12+start(0)*x16+start(0)*x17+start(0)*x18 \n
  f6=~x16 + start(0)*x17 + start(0)*x18+start(0)*x10+start(0)*x11+start(0)*x12+start(0)*x13+start(0)*x14+start(0)*x15 \n
  f4~start(0)*f1+start(0)*f2+start(0)*f3+ start(0)*f5+start(0)*f6
  + start(0)*x19 + start(0)*x20+ start(0)*x21\n
  f5~start(0)*f1+start(0)*f2+start(0)*f3+start(0)*f4+ start(0)*f6
  + start(0)*x19 + start(0)*x20+ start(0)*x21\n
  f6~start(0)*f1+start(0)*f2+start(0)*f3+ start(0)*f4+start(0)*f5
  + start(0)*x19 + start(0)*x20+ start(0)*x21\n
  x1~~start(",covcoef[1],")*x1 \n
  x2~~start(",covcoef[2],")*x2 \n
  x3~~start(",covcoef[3],")*x3 \n
  x4~~start(",covcoef[4],")*x4 \n
  x5~~start(",covcoef[5],")*x5 \n
  x6~~start(",covcoef[6],")*x6 \n
  x7~~start(",covcoef[7],")*x7 \n
  x8~~start(",covcoef[8],")*x8 \n
  x9~~start(",covcoef[9],")*x9 \n
  x10~~start(",covcoef[10],")*x10 \n
  x11~~start(",covcoef[11],")*x11 \n
  x12~~start(",covcoef[12],")*x12 \n
  x13~~start(",covcoef[13],")*x13 \n
  x14~~start(",covcoef[14],")*x14 \n
  x15~~start(",covcoef[15],")*x15 \n
  x16~~start(",covcoef[16],")*x16 \n
  x17~~start(",covcoef[17],")*x17 \n
  x18~~start(",covcoef[18],")*x18 \n
  f1~~start(",covcoef[19],")*f1 \n
  f2~~start(",covcoef[20],")*f2 \n
  f3~~start(",covcoef[21],")*f3 \n
  f4~~start(",covcoef[22],")*f4 \n
  f5~~start(",covcoef[23],")*f5 \n
  f6~~start(",covcoef[24],")*f6 \n
  f1~~start(0)*f2 \n
  f2~~start(0)*f3 \n
  f1~~start(0)*f3 \n
  x19~~start(",covcoef[25],")*x19 \n
  x20~~start(",covcoef[26],")*x20 \n
  x21~~start(",covcoef[27],")*x21 \n
  x19~~start(0) * x20 \n
  x20~~start(0)*x21 \n
  x19~~start(0)*x21 \n")



#############################end mod script

regsemt<-function(method){
  t0=proc.time()
  if(simumode==1){
    simu_gamma=3.7
  }
  if(simumode==2){
    simu_gamma=3.7
  }else{
    simu_gamma=3.7
  }
  if(method=="regsem1"){
    r1=tryCatch(
      regsem(lmod,type=Typechar,lambda=l,pars_pen=pars_pen_mod,full=FALSE,optMethod="coord_desc",Start=initx(),line.search=FALSE,momentum=FALSE,gamma=simu_gamma,max.iter = maxit)
      ,error=function(error){
        print("regsem error l=")
        print(l)
        print(error)
        return (NA)
      }
    )
  }
  if(method=="regsem2"){
    r1=tryCatch(
      regsem(lmod,type=Typechar,lambda=l,pars_pen=pars_pen_mod,full=FALSE,optMethod="coord_desc",Start=initx(),line.search=TRUE,momentum=FALSE,gamma=simu_gamma,max.iter = maxit)
      ,error=function(error){
        print("regsem error l=")
        print(l)
        print(error)
        return (NA)
      }
    )
  }
  if(method=="regsem3"){
    r1=tryCatch(
      regsem(lmod,type=Typechar,lambda=l,pars_pen=pars_pen_mod,optMethod="coord_desc",full=TRUE,line.search=FALSE,momentum=FALSE,Start=initx(),gamma=simu_gamma,max.iter = maxit)
      ,error=function(error){
        print("regsem error l=")
        print(l)
        print(error)
        return (NA)
      }
    )
  }
  if(method=="regsem4"){
    r1=tryCatch(
      regsem(lmod,type=Typechar,lambda=l,pars_pen=pars_pen_mod,optMethod="rsolnp",full=FALSE,line.search=FALSE,Start=initx(),gamma=simu_gamma,max.iter = maxit)
      ,error=function(error){
        print("regsem error l=")
        print(l)
        print(error)
        return (NA)
      }
    )
  }
  if(method=="regsem5"){
    r1=tryCatch(
      regsem(lmod,type=Typechar,lambda=l,pars_pen=pars_pen_mod,optMethod="coord_desc",full=FALSE,line.search=FALSE,momentum=TRUE,Start=initx(),gamma=simu_gamma,max.iter = maxit)
      ,error=function(error){
        print("regsem error l=")
        print(l)
        print(error)
        return (NA)
      }
    )
  }

  t0=proc.time()-t0



  status=0
  if(is.na(r1[1])){
    bic=Inf
    coef=numeric(0)
    fit=Inf
    status=-1
  }else{
    if(is.na(r1$fit)){
      bic=Inf
      coef=numeric(0)
      fit=Inf
      status=-2
    }else{
      coef=r1$coefficients


      if(r1$fit<0|r1$fit==Inf){
        status=-2
      }else{
        if(max(abs(coef))>1000){
          status=-3
        }
      }

      bic=computebic(fit=r1$fit,x=coef)

      fit=r1$fit

    }
  }

  list(method=method,coef=coef,fit=fit,bic=bic,l=l,t0=t0[1],r=status,tx=truex,pen=pars_pen_mod,nx=nx)

}


getbic<-function(method,lambda){
  l<<-lambda

  if(method=="regsem1"){#proximal gradient, line search=FALSE
    r=regsemt(method)

    return (r)
  }
  if(method=="regsem2"){#proximal grad, line search=TRUE
    r=regsemt(method)

    return (r)
  }
  if(method=="regsem3"){#palm, line search=FALSE
    r=regsemt(method)

    return (r)
  }

  if(method=="regsem4"){#PALM, line search=TRUE
    r=regsemt(method)

    return (r)
  }
  if(method=="regsem5"){#PALM, line search=TRUE
    r=regsemt(method)

    return (r)
  }
  if(method=="gradsamp"){
    if(simumode==1){
      simu_gamma=3.7
    }
    if(simumode==2){
      simu_gamma=3.7
    }else{
      simu_gamma=3.7
    }
    t0=proc.time()
    #suppressWarnings()
    r1=tryCatch(
      gscpp(initx(),A,S,S_fixed,A_fixed,A_est,S_est,F_b,I_b,dataS_m,pars_pen_mod,l,simumode,const_para_gradsamp_ra,maxit,simu_gamma,TRUE)
      ,error=function(error){
        print("gs error l=")
        print(l)
        print(error)
        return (NA)
      })
    t0=proc.time()-t0
    status=0
    if(is.na(r1[1])){
      bic=Inf
      coef=numeric(0)
      fit=Inf
      status=-1
    }else{
      coef=r1$coefficients
      fit=r1$fit
      if(fit<0|fit==Inf){
        status=-2
      }else{
        if(max(abs(coef))>1000){
          status=-3
        }
      }

      bic=computebic(fit=fit,x=coef)

    }
  }

  if(method=="lbfgs"){
    t0=proc.time()

    r1=tryCatch(
      lbfgs::lbfgs(fl,gl,vars=rev_shuffle(initx()),invisible=1,max_iterations=maxit,
                   orthantwise_c = l,orthantwise_start = 0 ,orthantwise_end = length(pars_pen_mod))
      ,error=function(error){
        print("lbfgs error l=")
        print(l)
        print(error)
        return (NA)
      })
    t0=proc.time()-t0
    status=0
    if(is.na(r1[1])){
      status=-1
      fit=Inf
      coef=numeric(0)
      bic=Inf
    }else{
      coef=shuffle( r1$par )
      fit=fobjml_0pen(coef)



      if(fit==Inf|fit<0){
        status=-2
      }else{
        if(max(abs(coef))>1000){
          status=-3
        }
      }
    }

    bic=computebic(fit=fit,x=coef)


  }

  if(method=="lsl"){#PALM, line search=TRUE

#For each pattern matrix, the element can be set as 1, 0, or NA. 1 means that the corresponding
#    parameter should be freely estimated, 0 indicates that the parameter sould be fixed, and NA
#    makes the parameter to be estimated with penalization
#the bic calculated by lsl should be multiplied by NDATA then it is the bic we use
    rc_sem <- lslSEM()
    rc_sem$input(raw = dataorgS)
    if(mod_num<9){
        lambda_v=inspect(lmod,'start')$lambda
    #lambda=inspect(lmod)$lambda
    #1 freely estimated ,NA penalized,0 fixed
        lambda <- matrix(0, 18, 6)
        lambda[c(1:9),c(1:3)]=1
        if(mod_num==3|mod_num==7){
          lambda[c(2:9), 1]<-NA
          lambda[c(1:3,5:9), 2]<-NA
          lambda[c(1:6,8:9), 3]<-NA
        }else{
          lambda[c(4:9), 1]<-NA
          lambda[c(1:3,7:9), 2]<-NA
          lambda[c(1:6), 3]<-NA
        }
        lambda[c(10:18),c(4:6)]=1
        lambda[1,1]<-lambda[4,2]<-lambda[7,3]<-0
        if(mod_num==3|mod_num==7){
            lambda[c(11:18), 4]<-NA
            lambda[c(10:12,14:18), 5]<-NA
            lambda[c(10:15,17:18), 6]<-NA
          }else{
            lambda[c(13:18), 4]<-NA
            lambda[c(10:12,16:18), 5]<-NA
            lambda[c(10:15), 6]<-NA
          }
  
        lambda[10,4]<-lambda[13,5]<-lambda[16,6]<-0



        #P=18
        #M=6



        beta=matrix(0,6,6)

        if(mod_num==3){

          beta[c(4,5,6),]<-NA
        }
         
        beta[4,4]<-beta[5,5]<-beta[6,6]<-0

        phi=matrix(0,6,6)
        phi[1:3,1:3]=1
        phi[4,4]=phi[5,5]=phi[6,6]=1
        psi=diag(18)


         
        nu=matrix(1,18,1)


        #value

        nu_v=matrix(0,18,1)
        tm=lav_model_set_parameters(lmod@Model,initx())
        tm=lmod@Model
        lambda_v=tm@GLIST$lambda
        beta_v=tm@GLIST$beta
        psi_v=tm@GLIST$theta
        phi_v = tm@GLIST$psi
      } 
    rc_sem$input(raw=dataorgS)
    rc_sem$specify(pattern = list(lambda = lambda,psi=psi,beta=beta,phi=phi,nu=nu),
                  value=list(lambda=lambda_v,beta=beta_v,psi=psi_v,phi=phi_v,nu=nu_v),
                  scale=FALSE)
    #rc_sem$specify(pattern = list(lambda = lambda,psi=psi,beta=beta,phi=phi,nu=nu),value=list(lambda=lambda_v),
#scale=FALSE)
    penalty_text=c('l1','scad','mcp')
    t0=proc.time()
    r1=tryCatch(
    c(1,rc_sem$learn(penalty = list(type = penalty_text[simumode+1], gamma = l),control=list(maxit=maxit))
    ),error=function(error){
      print("lsl error l=")
      print(l)
      print(error)
      return (NA)
    } )
    t0=proc.time()-t0
    if(is.na(r1[1])){
      status=-1
      fit=Inf
      coef=numeric(0)
      bic=Inf
    }else{

      lsl_r=rc_sem$summarize(type = "overall")
      coef_r=rc_sem$summarize(type = "individual")
      #order_coef=c(1:48,c(67,70,73,78,80),c(68,71,74,76,81),c(69,72,75,77,79),49:66,c(82,85,87,88,89,90,83,84,86))
      coef=coef_r[ ,2]
      fit=0.5*lsl_r$`bic optimal`[4]
      bic=NDATA*lsl_r$`bic optimal`[18]
      status=0
    }

  }

  if(method=='lslx'){
    if(mod_num==3|mod_num==7){


      lslxmod=paste0('fix(1) *x1  <=:f1
        x2 + x3 + x4 + x5 + x6+ x7 + x8+ x9 <~: f1
        fix(1) *x4 <=: f2
        x5 + x6 + x1+ x2 + x3 + x7+ x8 + x9 <~: f2
        fix(1) *x7 <=: f3
        x8 + x9 + x1 + x2 + x3+ x4 + x5 + x6 <~: f3
        fix(1) *x10  <=: f4
        x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 <~: f4
        fix(1) *x13  <=: f5
        x14 + x15 + x10 + x11 + x12 + x16 + x17 + x18 <~: f5
        fix(1) *x16  <=: f6
        x17 + x18 + x10 + x11 + x12 + x13 + x14 + x15 <~: f6
        f4<~ f1+f2+f3+f5+f6
        f5<~ f1+f2+f3+f4+f6
        f6<~ f1+f2+f3+f4+f5
        ')



    }
    if(mod_num==4){
      lslxmod=paste0('fix(1) *x1 + free()*x2 + free()*x3 <=: f1
      fix(1) *x4 + free()*x5 + free()*x6 <=: f2
      fix(1) *x7 + free()*x8 + free()*x9 <=: f3
       x4 + x5 + x6 + x7 + x8 + x9 <~: f1
       x1 + x2 + x3 + x7 + x8 + x9 <~: f2
       x4 + x5 + x6 + x7 + x8 + x9 <~: f3
      fix(1) *x10 + free()*x11 + free()*x12 <=: f4
      fix(1) *x13 + free()*x14 + free()*x15 <=: f5
      fix(1) *x16 + free()*x17 + free()*x18 <=: f6
       x13 + x14 + x15 + x16 + x17+ x18 <~: f4
       x10 + x11 + x12 + x16 + x17 + x18 <~: f5
       x10 + x11 + x12 + x13 + x14 + x15 <~: f6
      f4<= f1 +   f2 +   f3 + free()*f5 +  f6
      f5<= f1+free()*f2+ f3+ f4+ f6
      f6<~ f1+ f2+ f3+ f4+ f5

      ')
    }
      if(mod_num==8){
        lslxmod=paste0('fix(1) *x1 + free()*x2 + free()*x3 <=: f1
      fix(1) *x4 + free()*x5 + free()*x6 <=: f2
       fix(1) *x7 + free()*x8 + free()*x9 <=: f3
       x4 + x5 + x6 + x7 + x8 + x9 <~: f1
       x1 + x2 + x3 + x7 + x8 + x9 <~: f2
       x4 + x5 + x6 + x7 + x8 + x9 <~: f3
       fix(1) *x10 + free()*x11 + free()*x12 <=: f4
       fix(1) *x13 + free()*x14 + free()*x15 <=: f5
       fix(1) *x16 + free()*x17 + free()*x18 <=: f6
       x13 + x14 + x15 + x16 + x17+ x18 <~: f4
       x10 + x11 + x12 + x16 + x17 + x18 <~: f5
       x10 + x11 + x12 + x13 + x14 + x15 <~: f6
       f4<= f1 +   f2 +   f3 + f5 +  f6
       f5<= f1+free()*f2+ f3+ f4+ f6
       f6<~ f1+ f2+ f3+ f4+ f5
       ')
      }

    if(mod_num==9){
      lslxmod=paste0('fix(1) *x1 + free(0)*x2 + free(0)*x3 <=: f1
      fix(1) *x4 + free()*x5 + free()*x6 <=: f2
      fix(1) *x7 + free()*x8 + free()*x9 <=: f3
      x4 + x5 + x6 + x7 + x8 + x9 <~: f1
      x1 + x2 + x3 + x7 + x8 + x9 <~: f2
      x4 + x5 + x6 + x7 + x8 + x9 <~: f3
      fix(1) *x10 + free()*x11 + free()*x12 <=: f4
      fix(1) *x13 + free()*x14 + free()*x15 <=: f5
      fix(1) *x16 + free()*x17 + free()*x18 <=: f6
      x13 + x14 + x15 + x16 + x17+ x18 <~: f4
      x10 + x11 + x12 + x16 + x17 + x18 <~: f5
      x10 + x11 + x12 + x13 + x14 + x15 <~: f6
      f4<= f1 +   f2 +   f3 + free()*f5 +  f6 + x19 + x20 + x21
      f5<= f1+free()*f2+ f3+ f4+ f6 + x19 + x20 + x21
      f6<~ f1+ f2+ f3+ f4+ f5 + x19 + x20 + x21
      x19~~0*x20
      x20~~0*x21
      x19~~0*x21
      x19+x20+x21~~0*f1+0*f2+0*f3
      ')

    }

    if(mod_num==10){
      lslxmod=paste0('fix(1) *x1  <=:f1
      x2 + x3 + x4 + x5 + x6+ x7 + x8+ x9   <~: f1
       fix(1) *x4 <=: f2
       x5 + x6 + x1+ x2 + x3 + x7+ x8 + x9 <~: f2
       fix(1) *x7 <=: f3
       x8 + x9 + x1 + x2 + x3+ x4 + x5 + x6 <~: f3
       fix(1) *x10  <=: f4
       x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 <~: f4
       fix(1) *x13  <=: f5
       x14 + x15 + x10 + x11 + x12 + x16 + x17 + x18<~: f5
       fix(1) *x16  <=: f6
       x17 + x18 + x10 + x11 + x12 + x13 + x14 + x15<~: f6
      f4<= f1 +   f2 +   f3 + free()*f5 +  f6 + x19 + x20 + x21
      f5<= f1+free()*f2+ f3+ f4+ f6 + x19 + x20 + x21
                     f6<~ f1+ f2+ f3+ f4+ f5 + x19 + x20 + x21
                     x19~~0*x20
                     x20~~0*x21
                     x19~~0*x21
                     x19+x20+x21~~0*f1+0*f2+0*f3
                     ')

    }
    t0=proc.time()
    lslx_reg <- lslx$new(model = lslxmod,
                         data = dataorgS)
    #lslx_reg <- lslx$new(model = lslxmod,
    #                     sample_cov = dataS,
    #                     sample_size= nrow(dataorgS) )


    r1=tryCatch(
      c(1,lslx_reg$fit(penalty_method = Typechar,
                   lambda_grid = c(l),iter_out_max=maxit)),error=function(error){
        print("lslx error l=")
        print(l)
        print(error)
        return (NA)
      } )
    t0 = proc.time()-t0
    if(is.na(r1[1])){
      status=-1
      fit=Inf
      coef=numeric(0)
      bic=Inf
    }else{

      #lslx_reg$extract_coefficient_matrix(selector = "bic",include_faulty = TRUE)
      #lslx_reg$extract_coefficient_matrix(block='f<-1',selector = "bic",include_faulty = TRUE)
      #lslx_reg$extract_coefficient_matrix(block='y<-1',selector = "bic",include_faulty = TRUE)
      #slx_reg$extract_coefficient_matrix(block='f<-f',selector = "bic",include_faulty = TRUE)
      #lslx_reg$extract_coefficient_matrix(block='f<-y',selector = "bic",include_faulty = TRUE)
      #lslx_reg$extract_coefficient_matrix(block='y<-f',selector = "bic",include_faulty = TRUE)
      #lslx_reg$extract_coefficient_matrix(block='y<-y',selector = "bic",include_faulty = TRUE)
      #lslx_reg$extract_coefficient_matrix(block='f<->f',selector = "bic",include_faulty = TRUE)
      #lslx_reg$extract_coefficient_matrix(block='f<->y',selector = "bic",include_faulty = TRUE)
      #lslx_reg$extract_coefficient_matrix(block='y<->f',selector = "bic",include_faulty = TRUE)
      #lslx_reg$extract_coefficient_matrix(block='y<->y',selector = "bic",include_faulty = TRUE)

      #lslx_reg$get_fitting()$fitted_result
      coef=lslx_reg$extract_coefficient(selector = "bic",include_faulty = TRUE)
      fit=as.numeric(lslx_reg$get_fitting()$fitted_result$numerical_condition[[1]][7])
      bic=computebic(fit,coef,TRUE)
      if(Typechar=='lasso'){
        status=-1+as.numeric(lslx_reg$get_fitting()$fitted_result$is_convergent)
      }else{
        status=-1+as.numeric(lslx_reg$get_fitting()$fitted_result$is_convergent[2])
      }
    }
  }
  return(
    list(method=method,coef=as.double(coef),fit=fit,bic=bic,l=l,t0=t0[1],r=status,tx=truex,pen=pars_pen_mod,nx=nx)
  )

}







bicseq<-function(lseq=SEQU){
  method_text=c()
  if(Typechar=='lasso'){
    method_text[1]="regsem1"
    method_text[2]="regsem2"
    method_text[3]="regsem3"
    method_text[4]="regsem4"
    method_text[5]="regsem5"
    method_text[6]="gradsamp"
    method_text[7]="lbfgs"
    method_text[8]='lslx'
    method_text[9]='lsl'
  }
  if(Typechar=='scad'){
    method_text[1]="regsem1"
    method_text[2]="regsem2"
    method_text[3]="regsem3"
    method_text[4]="regsem4"
    method_text[5]="regsem5"
    method_text[6]="gradsamp"
    method_text[7]='lsl'
     
  }
  if(Typechar=='mcp'){
    method_text[1]="regsem1"
    method_text[2]="regsem2"
    method_text[3]="regsem3"
    method_text[4]="regsem4"
    method_text[5]="regsem5"
    method_text[6]="gradsamp"
    method_text[7]='lslx'
    method_text[8]='lsl'
  }
   
  result=list()
  for(i in 1:length(method_text)){
    result[[i]]=list()
    for(j in 1:length(lseq)){
      result[[i]][[j]]=getbic(method_text[i],lseq[j])
    }
  }
  return (result)
}


#########functions end








#SEQU must go top down


#simumode=1



###########make model

##############mod script


###make model running

updateall<-function(nd=NDATA){
  lmod<<-sem(L3[[mod_num]],sample.nobs=nd,sample.cov=dataS,
             sample.cov.rescale=TRUE,int.ov.free=TRUE,
             int.lv.free=FALSE,fixed.x=FALSE,std.lv=FALSE,
             do.fit=T,meanstructure=TRUE,test='standard')
#so.fit=T for calculating df using fitmeasures(lmod,'df')



  mats=extractMatrices(lmod)
  A <<- mats$A
  A_est <<- mats$A_est
  A_fixed <<- mats$A_fixed
  S <<- mats$S
  S_est <<- mats$S_est
  S_fixed <<- mats$S_fixed
  F_b <<- mats$F
  I_b <<- diag(nrow(A))
  #A_pen = A != 0
  if(mod_num==3|mod_num==7|mod_num==10){
    rind=which(rownames(A)=='1')
    A_pen = A[-rind,-rind]
    pars_pen_mod <<- A_pen[A_pen != 0]
  }else{
    if(mod_num!=8){
      rind=which(rownames(A)=='1')
      A_pen = A[-rind,-rind]
      pars_pen_mod <<- setdiff(A_pen[A_pen != 0],c(A_pen['x2','f1'],A_pen['x3','f1'],A_pen['x5','f2'],A_pen['x6','f2'],
      A_pen['x8','f3'],A_pen['x9','f3'],A_pen['x11','f4'],A_pen['x12','f4'],A_pen['x14','f5'],A_pen['x15','f5'],
      A_pen['x17','f6'],A_pen['x18','f6'],A_pen['f4','f5'],A_pen['f5','f2']) )
    }else{
      rind=which(rownames(A)=='1')
      A_pen = A[-rind,-rind]
      pars_pen_mod <<- setdiff(A_pen[A_pen != 0],c(A_pen['x2','f1'],A_pen['x3','f1'],A_pen['x5','f2'],A_pen['x6','f2'],
      A_pen['x8','f3'],A_pen['x9','f3'],A_pen['x11','f4'],A_pen['x12','f4'],A_pen['x14','f5'],A_pen['x15','f5'],
      A_pen['x17','f6'],A_pen['x18','f6'],A_pen['f5','f2']) )
    }

  }
  #trying to modify dataS if meanstructure==TRUE
  if(mats$mean == TRUE){
    mm = mats$A[,"1"]

    SampMean <- lmod@SampleStats@mean[][[1]]
    ss = match(names(mm[mm > 0]),lmod@Data@ov$name)
    ss <- ss[!is.na(ss)] # NAs may arise from fixed means
    SampMean[-c(ss)] = 0
    SampCov <- lmod@SampleStats@cov[][[1]]

    SampCov2 <- SampCov + SampMean%*%t(SampMean)
    # try changing size of SampCov
    SampCov3 = cbind(SampCov2,SampMean)
    dataS_m <<- rbind(SampCov3,append(SampMean,1))
  }else{
    dataS_m <<- dataS
    colnames(dataS_m)=names(dataorgS)
    rownames(dataS_m)=names(dataorgS)
  }
}

if(mod_num<11){
  dataorgS=simulateData(L2[[mod_num]],sample.nobs=NDATA,meanstructure = TRUE)
  dataS<<-cov(dataorgS)


  colnames(dataS)=names(dataorgS)
  dataS0<<-dataS

  updateall(NDATA)
  #data for cv
  covS=list()
  covrS=list()
  nd=NDATA

  for(i in 1:5){

    covrS[[i]]=cov(dataorgS[-c((nd*(i-1)/5+1):(nd*i/5)),])
    covS[[i]]=cov(dataorgS[c((nd*(i-1)/5+1):(nd*i/5)),])
  }
  l=NULL
  #nx=max(max(A),max(S))
  nx=fitmeasures(lmod,'df')
  tmod=suppressWarnings(
    sem(L1[[mod_num]],sample.nobs=NDATA,sample.cov=dataS,sample.cov.rescale=TRUE,fixed.x=FALSE,
        do.fit=FALSE,control=list(iter.max=0),meanstructure=TRUE))
  truex=coef(tmod)
}else{
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
    lmod<<-cfa(mod[1],data=dat[[1]],fixed.x=FALSE,meanstructure=TRUE)
  }else{
    lmod<<-sem(mod[2], data = pisa.be.2003, std.lv = FALSE,
               fixed.x=TRUE,conditional.x=TRUE,
               int.ov.free=TRUE,int.lv.free=FALSE,meanstructure = TRUE)
  }
  mats=extractMatrices(lmod)
  A = mats$A
  A_est = mats$A_est
  A_fixed = mats$A_fixed
  S = mats$S
  S_est = mats$S_est
  S_fixed = mats$S_fixed
  F_b <<- mats$F
  I_b <<- diag(nrow(A))
  #A_pen = A != 0
  rind=which(rownames(A)=='1')
  A_pen = A[-rind,-rind]
  pars_pen_mod <<- A_pen[A_pen != 0]
  
  NDATA=nrow(dat[[mod_num-10]])
  dataS=NDATA/(NDATA-1)*cov(dat[[mod_num-10]])
  nx=fitmeasures(lmod,'df')
  covS=list()
  covrS=list()
  nd=NDATA
  dataorgS=dat[[mod_num-10]]
  for(i in 1:5){
    covrS[[i]]=cov(dataorgS[-c((nd*(i-1)/5+1):(nd*i/5)),])
    covS[[i]]=cov(dataorgS[c((nd*(i-1)/5+1):(nd*i/5)),])
  }
  truex=NULL
}
####init globval





##################################load my functions


result1=bicseq()
 
save(result1,truex,pars_pen_mod
     ,file=paste0('paper_simu/',MYSEED,"_",NDATA,"_",maxit,"_",100*VFAC,"_",simumode,"_",mod_num,".RData"))




