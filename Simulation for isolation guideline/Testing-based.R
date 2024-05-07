install.packages("deSolve")
install.packages("writexl")

library(ggplot2)
library(dplyr)
library(deSolve)
library(writexl)
library(readxl)


setwd("~/Desktop/")


#################### Setting ####################
Tmin <- 0
Tmax <- 100
step_size <- 1
mtime <- seq(Tmin,Tmax,step_size)

N <- 1000 ## # of samples

First <- 0 ## First testing

pop <- read.csv("populationParameters_mpox.txt", row.names = 1)

##################################### Viral dynamics model #####################################
Covfun<-function(pars){
  
  delta <- as.numeric(pars[1])
  V0 <- as.numeric(pars[2])
  
  derivs<-function(time,y,pars){
    with(as.list(c(pars,y)),{
      
      dV <- -delta*V
      
      return(list(c(dV)))
    })
  }
  y<-c(V=10^V0)
  
  times<-c(seq(Tmin,Tmax,step_size))
  out<-lsoda(y=y,parms=pars,times=times,func=derivs,rtol=0.00004,atol=0)
  out2<-cbind(time=out[,1],aV=((log10(out[,2]))))
  as.data.frame(out2)
}




######################################################################################
######################################################################################
##################################### Simulation #####################################
######################################################################################
######################################################################################

#### Risk of prematurely ending isolation (Probability => P)
Prob11 <- data.frame(); Prob21 <- data.frame(); Prob31 <- data.frame()
Prob12 <- data.frame(); Prob22 <- data.frame(); Prob32 <- data.frame()
Prob13 <- data.frame(); Prob23 <- data.frame(); Prob33 <- data.frame()
Prob14 <- data.frame(); Prob24 <- data.frame(); Prob34 <- data.frame()
Prob15 <- data.frame(); Prob25 <- data.frame(); Prob35 <- data.frame()

#### Infectious period after ending isolation (Length => L)
Leng11 <- data.frame(); Leng21 <- data.frame(); Leng31 <- data.frame()
Leng12 <- data.frame(); Leng22 <- data.frame(); Leng32 <- data.frame()
Leng13 <- data.frame(); Leng23 <- data.frame(); Leng33 <- data.frame()
Leng14 <- data.frame(); Leng24 <- data.frame(); Leng34 <- data.frame()
Leng15 <- data.frame(); Leng25 <- data.frame(); Leng35 <- data.frame()

#### Unnecessarily prolonged isolation period (Burden => B)
Burd11 <- data.frame(); Burd21 <- data.frame(); Burd31 <- data.frame()
Burd12 <- data.frame(); Burd22 <- data.frame(); Burd32 <- data.frame()
Burd13 <- data.frame(); Burd23 <- data.frame(); Burd33 <- data.frame()
Burd14 <- data.frame(); Burd24 <- data.frame(); Burd34 <- data.frame()
Burd15 <- data.frame(); Burd25 <- data.frame(); Burd35 <- data.frame()

#### Corresponding isolation period (Isolation => I)
Isol11 <- data.frame(); Isol21 <- data.frame(); Isol31 <- data.frame()
Isol12 <- data.frame(); Isol22 <- data.frame(); Isol32 <- data.frame()
Isol13 <- data.frame(); Isol23 <- data.frame(); Isol33 <- data.frame()
Isol14 <- data.frame(); Isol24 <- data.frame(); Isol34 <- data.frame()
Isol15 <- data.frame(); Isol25 <- data.frame(); Isol35 <- data.frame()


simulation <- 1


for (s in 1:simulation) {
  
  
  while (TRUE) {
  #################### Parameters ####################
  
  popdelta <- rlnorm(N, log(pop["delta_pop", "value"]), pop["omega_delta", "value"])
  popv0    <- rlnorm(N, log(pop["V0_pop", "value"]), pop["omega_V0", "value"])
  fitt     <- data.frame(delta=popdelta,V0=popv0)
  
  
  #################### Sampling 1000 COVID-19 patients ####################
  IF <- log10(10^5.5)
  
  df_total = data.frame()
  
  for ( i in 1:N ) {
    
    pars <- c(delta=fitt$delta[i],V_0=fitt$V0[i])
    out  <- Covfun(pars)
    
    ppp<-matrix(NA,nrow=length(mtime),ncol=1)
    for(jj in 1:length(mtime)){
      ii<-jj
      a<-mtime[ii]
      dd<-out[out$time==a,]
      ppp[ii]<-dd$aV
    }
    pp<-data.frame(ID=i, time=mtime, V=ppp)
    df_total<-rbind(df_total,pp)
    
  }
  
  
  if ( length(which(df_total$V < IF & df_total$time == 0))==0 & length(which(df_total$V > IF & df_total$time == Tmax))==0 )
    break
  
  }

  
  #################### Including error for viral dynamics of COVID-19 patients ####################
  T  <- (Tmax+1)*N   
  TM <- T+Tmax+1
  error <- rnorm(T, mean=0, sd=pop$value[5])
  
  df_totalM  <- df_total
  df_totalM2 <- df_totalM$V+error
  df_totalM$obs <- df_totalM2 ## Data frame of total patients from viral dynamics model
  
  nomeaning <- data.frame(ID=rep(N+1,times=Tmax+1),time=mtime,V=rep(10,times=Tmax+1),obs=rep(10,times=Tmax+1))
  df_totalM <- rbind(df_totalM,nomeaning)
  
  
  
  
  
  
  
  
  
  
  #########################################################################################################################
  #########################################################################################################################
  #########################################################################################################################
  ###################################### Infectiousness threshold = 10^5.5 copies/ml ######################################
  #########################################################################################################################
  #########################################################################################################################
  #########################################################################################################################
  
  IF <- log10(10^5.5)
  DL <- IF
  
  ################### Time point where the patient loses infectiousness ###################
  cri <- c()
  
  for ( i in 1:TM ) {
    if (df_totalM$V[i]<IF) {
      cri[i] <- df_totalM$time[i]
    }else {
      cri[i] <- Tmax+1
    }
  }
  
  df_totalM   <-cbind(df_totalM,cri)
  df_totalM$V <- as.numeric(df_totalM$V)
  df_totalM$obs <- as.numeric(df_totalM$obs)
  
  
  

  ##############################################################################################################
  ######################################## 1 Consecutive negative result #######################################
  ##############################################################################################################
  
  K <- 1
  
  ############################################################ Interval of tests = 1 day ############################################################
  P <- c(); PP <- c(); 
  L <- c(); LL <- c(); 
  B <- c(); BB <- c(); 
  I <- c(); II <- c(); 
  
  k <- 1
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        # & i%%k==1
        #& df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL 
        #& df_totalM$obs[i+3*k]<DL 
        #& df_totalM$obs[i+4*k]<DL 
        #& df_totalM$time[i]<df_totalM$time[i+1*k] 
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k] 
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 2 day ############################################################
  P <- c(); PP <- c(); 
  L <- c(); LL <- c(); 
  B <- c(); BB <- c(); 
  I <- c(); II <- c(); 
  
  k <- 2
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & i%%k==1
        #& df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL 
        #& df_totalM$obs[i+3*k]<DL 
        #& df_totalM$obs[i+4*k]<DL 
        #& df_totalM$time[i]<df_totalM$time[i+1*k] 
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k] 
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 3 day ############################################################
  P <- c(); PP <- c(); 
  L <- c(); LL <- c(); 
  B <- c(); BB <- c(); 
  I <- c(); II <- c(); 
  
  k <- 3
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & i%%k==1
        #& df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL 
        #& df_totalM$obs[i+3*k]<DL 
        #& df_totalM$obs[i+4*k]<DL 
        #& df_totalM$time[i]<df_totalM$time[i+1*k] 
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k] 
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 4 day ############################################################
  P <- c(); PP <- c(); 
  L <- c(); LL <- c(); 
  B <- c(); BB <- c(); 
  I <- c(); II <- c(); 
  
  k <- 4
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & i%%k==1
        #& df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL 
        #& df_totalM$obs[i+3*k]<DL 
        #& df_totalM$obs[i+4*k]<DL 
        #& df_totalM$time[i]<df_totalM$time[i+1*k] 
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k] 
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 5 day ############################################################
  P <- c(); PP <- c(); 
  L <- c(); LL <- c(); 
  B <- c(); BB <- c(); 
  I <- c(); II <- c(); 
  
  k <- 5
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & i%%k==1
        #& df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL 
        #& df_totalM$obs[i+3*k]<DL 
        #& df_totalM$obs[i+4*k]<DL 
        #& df_totalM$time[i]<df_totalM$time[i+1*k] 
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k] 
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  #################### Result ####################
  
  Prob11 <- rbind(Prob11,data.frame(P1=sum(P1)*100/N,P2=sum(P2)*100/N,P3=sum(P3)*100/N,P4=sum(P4)*100/N,P5=sum(P5)*100/N))
  Leng11 <- rbind(Leng11,data.frame(L1=mean(L1),L2=mean(L2),L3=mean(L3),L4=mean(L4),L5=mean(L5)))
  Burd11 <- rbind(Burd11,data.frame(B1=mean(B1),B2=mean(B2),B3=mean(B3),B4=mean(B4),B5=mean(B5)))
  Isol11 <- rbind(Isol11,data.frame(I1=mean(I1),I2=mean(I2),I3=mean(I3),I4=mean(I4),I5=mean(I5)))
  
  
  
  
  ##############################################################################################################
  ######################################## 2 Consecutive negative result #######################################
  ##############################################################################################################
  
  K <- 2
  
  ############################################################ Interval of tests = 1 day ############################################################
  P <- c(); PP <- c(); 
  L <- c(); LL <- c(); 
  B <- c(); BB <- c(); 
  I <- c(); II <- c(); 
  
  k <- 1
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL 
        #& df_totalM$obs[i+3*k]<DL 
        #& df_totalM$obs[i+4*k]<DL 
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k] 
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 2 day ############################################################
  P <- c(); PP <- c(); 
  L <- c(); LL <- c(); 
  B <- c(); BB <- c(); 
  I <- c(); II <- c(); 
  
  k <- 2
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL 
        #& df_totalM$obs[i+3*k]<DL 
        #& df_totalM$obs[i+4*k]<DL 
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k] 
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 3 day ############################################################
  P <- c(); PP <- c(); 
  L <- c(); LL <- c(); 
  B <- c(); BB <- c(); 
  I <- c(); II <- c(); 
  
  k <- 3
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL 
        #& df_totalM$obs[i+3*k]<DL 
        #& df_totalM$obs[i+4*k]<DL 
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k] 
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 4 day ############################################################
  P <- c(); PP <- c(); 
  L <- c(); LL <- c(); 
  B <- c(); BB <- c(); 
  I <- c(); II <- c(); 
  
  k <- 4
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL 
        #& df_totalM$obs[i+3*k]<DL 
        #& df_totalM$obs[i+4*k]<DL 
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k] 
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 5 day ############################################################
  P <- c(); PP <- c(); 
  L <- c(); LL <- c(); 
  B <- c(); BB <- c(); 
  I <- c(); II <- c(); 
  
  k <- 5
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL 
        #& df_totalM$obs[i+3*k]<DL 
        #& df_totalM$obs[i+4*k]<DL 
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k] 
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  #################### Result ####################
  
  Prob12 <- rbind(Prob12,data.frame(P1=sum(P1)*100/N,P2=sum(P2)*100/N,P3=sum(P3)*100/N,P4=sum(P4)*100/N,P5=sum(P5)*100/N))
  Leng12 <- rbind(Leng12,data.frame(L1=mean(L1),L2=mean(L2),L3=mean(L3),L4=mean(L4),L5=mean(L5)))
  Burd12 <- rbind(Burd12,data.frame(B1=mean(B1),B2=mean(B2),B3=mean(B3),B4=mean(B4),B5=mean(B5)))
  Isol12 <- rbind(Isol12,data.frame(I1=mean(I1),I2=mean(I2),I3=mean(I3),I4=mean(I4),I5=mean(I5)))
  
  
  
  
  ##############################################################################################################
  ######################################## 3 Consecutive negative result #######################################
  ##############################################################################################################
  
  K <- 3
  
  ############################################################ Interval of tests = 1 day ############################################################
  P <- c(); PP <- c(); 
  L <- c(); LL <- c(); 
  B <- c(); BB <- c(); 
  I <- c(); II <- c(); 
  
  k <- 1
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL 
        #& df_totalM$obs[i+4*k]<DL 
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 2 day ############################################################
  P <- c(); PP <- c(); 
  L <- c(); LL <- c(); 
  B <- c(); BB <- c(); 
  I <- c(); II <- c(); 
  
  k <- 2
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL 
        #& df_totalM$obs[i+4*k]<DL 
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 3 day ############################################################
  P <- c(); PP <- c(); 
  L <- c(); LL <- c(); 
  B <- c(); BB <- c(); 
  I <- c(); II <- c(); 
  
  k <- 3
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL 
        #& df_totalM$obs[i+4*k]<DL 
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 4 day ############################################################
  P <- c(); PP <- c(); 
  L <- c(); LL <- c(); 
  B <- c(); BB <- c(); 
  I <- c(); II <- c(); 
  
  k <- 4
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL 
        #& df_totalM$obs[i+4*k]<DL 
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 5 day ############################################################
  P <- c(); PP <- c(); 
  L <- c(); LL <- c(); 
  B <- c(); BB <- c(); 
  I <- c(); II <- c(); 
  
  k <- 5
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL 
        #& df_totalM$obs[i+4*k]<DL 
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  #################### Result ####################
  
  Prob13 <- rbind(Prob13,data.frame(P1=sum(P1)*100/N,P2=sum(P2)*100/N,P3=sum(P3)*100/N,P4=sum(P4)*100/N,P5=sum(P5)*100/N))
  Leng13 <- rbind(Leng13,data.frame(L1=mean(L1),L2=mean(L2),L3=mean(L3),L4=mean(L4),L5=mean(L5)))
  Burd13 <- rbind(Burd13,data.frame(B1=mean(B1),B2=mean(B2),B3=mean(B3),B4=mean(B4),B5=mean(B5)))
  Isol13 <- rbind(Isol13,data.frame(I1=mean(I1),I2=mean(I2),I3=mean(I3),I4=mean(I4),I5=mean(I5)))
  
  
  
  
  ##############################################################################################################
  ######################################## 4 Consecutive negative result #######################################
  ##############################################################################################################
  
  K <- 4
  
  ############################################################ Interval of tests = 1 day ############################################################
  P <- c(); PP <- c(); 
  L <- c(); LL <- c(); 
  B <- c(); BB <- c(); 
  I <- c(); II <- c(); 
  
  k <- 1
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL 
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 2 day ############################################################
  P <- c(); PP <- c(); 
  L <- c(); LL <- c(); 
  B <- c(); BB <- c(); 
  I <- c(); II <- c(); 
  
  k <- 2
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL 
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 3 day ############################################################
  P <- c(); PP <- c(); 
  L <- c(); LL <- c(); 
  B <- c(); BB <- c(); 
  I <- c(); II <- c(); 
  
  k <- 3
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL 
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 4 day ############################################################
  P <- c(); PP <- c(); 
  L <- c(); LL <- c(); 
  B <- c(); BB <- c(); 
  I <- c(); II <- c(); 
  
  k <- 4
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL 
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 5 day ############################################################
  P <- c(); PP <- c(); 
  L <- c(); LL <- c(); 
  B <- c(); BB <- c(); 
  I <- c(); II <- c(); 
  
  k <- 5
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL 
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  #################### Result ####################
  
  Prob14 <- rbind(Prob14,data.frame(P1=sum(P1)*100/N,P2=sum(P2)*100/N,P3=sum(P3)*100/N,P4=sum(P4)*100/N,P5=sum(P5)*100/N))
  Leng14 <- rbind(Leng14,data.frame(L1=mean(L1),L2=mean(L2),L3=mean(L3),L4=mean(L4),L5=mean(L5)))
  Burd14 <- rbind(Burd14,data.frame(B1=mean(B1),B2=mean(B2),B3=mean(B3),B4=mean(B4),B5=mean(B5)))
  Isol14 <- rbind(Isol14,data.frame(I1=mean(I1),I2=mean(I2),I3=mean(I3),I4=mean(I4),I5=mean(I5)))
  
  
  
  
  ##############################################################################################################
  ######################################## 5 Consecutive negative result #######################################
  ##############################################################################################################
  
  K <- 5
  
  ############################################################ Interval of tests = 1 day ############################################################
  P <- c(); PP <- c(); 
  L <- c(); LL <- c(); 
  B <- c(); BB <- c(); 
  I <- c(); II <- c(); 
  
  k <- 1
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        & df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        & df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 2 day ############################################################
  P <- c(); PP <- c(); 
  L <- c(); LL <- c(); 
  B <- c(); BB <- c(); 
  I <- c(); II <- c(); 
  
  k <- 2
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        & df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        & df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 3 day ############################################################
  P <- c(); PP <- c(); 
  L <- c(); LL <- c(); 
  B <- c(); BB <- c(); 
  I <- c(); II <- c(); 
  
  k <- 3
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        & df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        & df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 4 day ############################################################
  P <- c(); PP <- c(); 
  L <- c(); LL <- c(); 
  B <- c(); BB <- c(); 
  I <- c(); II <- c(); 
  
  k <- 4
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        & df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        & df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 5 day ############################################################
  P <- c(); PP <- c(); 
  L <- c(); LL <- c(); 
  B <- c(); BB <- c(); 
  I <- c(); II <- c(); 
  
  k <- 5
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        & df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        & df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  #################### Result ####################
  
  Prob15 <- rbind(Prob15,data.frame(P1=sum(P1)*100/N,P2=sum(P2)*100/N,P3=sum(P3)*100/N,P4=sum(P4)*100/N,P5=sum(P5)*100/N))
  Leng15 <- rbind(Leng15,data.frame(L1=mean(L1),L2=mean(L2),L3=mean(L3),L4=mean(L4),L5=mean(L5)))
  Burd15 <- rbind(Burd15,data.frame(B1=mean(B1),B2=mean(B2),B3=mean(B3),B4=mean(B4),B5=mean(B5)))
  Isol15 <- rbind(Isol15,data.frame(I1=mean(I1),I2=mean(I2),I3=mean(I3),I4=mean(I4),I5=mean(I5)))
  
  
  df_totalM <- subset(df_totalM, select=-cri)
  
  
  
  
  
  
  
  
  
  #########################################################################################################################
  #########################################################################################################################
  #########################################################################################################################
  ###################################### Infectiousness threshold = 10^6.0 copies/ml ######################################
  #########################################################################################################################
  #########################################################################################################################
  #########################################################################################################################
  
  IF <- log10(10^6.0)
  DL <- IF
  
  ################### Time point where the patient loses infectiousness ###################
  cri <- c()
  
  for ( i in 1:TM ) {
    if (df_totalM$V[i]<IF) {
      cri[i] <- df_totalM$time[i]
    }else {
      cri[i] <- Tmax+1
    }
  }
  
  df_totalM   <-cbind(df_totalM,cri)
  df_totalM$V <- as.numeric(df_totalM$V)
  df_totalM$obs <- as.numeric(df_totalM$obs)
  
  
  
  
  ##############################################################################################################
  ######################################## 1 Consecutive negative result #######################################
  ##############################################################################################################
  
  K <- 1
  
  ############################################################ Interval of tests = 1 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 1
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        # & i%%k==1
        #& df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        #& df_totalM$time[i]<df_totalM$time[i+1*k]
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 2 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 2
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & i%%k==1
        #& df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        #& df_totalM$time[i]<df_totalM$time[i+1*k]
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 3 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 3
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & i%%k==1
        #& df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        #& df_totalM$time[i]<df_totalM$time[i+1*k]
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 4 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 4
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & i%%k==1
        #& df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        #& df_totalM$time[i]<df_totalM$time[i+1*k]
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 5 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 5
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & i%%k==1
        #& df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        #& df_totalM$time[i]<df_totalM$time[i+1*k]
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  #################### Result ####################
  
  Prob21 <- rbind(Prob21,data.frame(P1=sum(P1)*100/N,P2=sum(P2)*100/N,P3=sum(P3)*100/N,P4=sum(P4)*100/N,P5=sum(P5)*100/N))
  Leng21 <- rbind(Leng21,data.frame(L1=mean(L1),L2=mean(L2),L3=mean(L3),L4=mean(L4),L5=mean(L5)))
  Burd21 <- rbind(Burd21,data.frame(B1=mean(B1),B2=mean(B2),B3=mean(B3),B4=mean(B4),B5=mean(B5)))
  Isol21 <- rbind(Isol21,data.frame(I1=mean(I1),I2=mean(I2),I3=mean(I3),I4=mean(I4),I5=mean(I5)))
  
  
  
  
  ##############################################################################################################
  ######################################## 2 Consecutive negative result #######################################
  ##############################################################################################################
  
  K <- 2
  
  ############################################################ Interval of tests = 1 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 1
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 2 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 2
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 3 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 3
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 4 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 4
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 5 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 5
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  #################### Result ####################
  
  Prob22 <- rbind(Prob22,data.frame(P1=sum(P1)*100/N,P2=sum(P2)*100/N,P3=sum(P3)*100/N,P4=sum(P4)*100/N,P5=sum(P5)*100/N))
  Leng22 <- rbind(Leng22,data.frame(L1=mean(L1),L2=mean(L2),L3=mean(L3),L4=mean(L4),L5=mean(L5)))
  Burd22 <- rbind(Burd22,data.frame(B1=mean(B1),B2=mean(B2),B3=mean(B3),B4=mean(B4),B5=mean(B5)))
  Isol22 <- rbind(Isol22,data.frame(I1=mean(I1),I2=mean(I2),I3=mean(I3),I4=mean(I4),I5=mean(I5)))
  
  
  
  
  ##############################################################################################################
  ######################################## 3 Consecutive negative result #######################################
  ##############################################################################################################
  
  K <- 3
  
  ############################################################ Interval of tests = 1 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 1
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 2 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 2
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 3 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 3
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 4 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 4
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 5 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 5
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  #################### Result ####################
  
  Prob23 <- rbind(Prob23,data.frame(P1=sum(P1)*100/N,P2=sum(P2)*100/N,P3=sum(P3)*100/N,P4=sum(P4)*100/N,P5=sum(P5)*100/N))
  Leng23 <- rbind(Leng23,data.frame(L1=mean(L1),L2=mean(L2),L3=mean(L3),L4=mean(L4),L5=mean(L5)))
  Burd23 <- rbind(Burd23,data.frame(B1=mean(B1),B2=mean(B2),B3=mean(B3),B4=mean(B4),B5=mean(B5)))
  Isol23 <- rbind(Isol23,data.frame(I1=mean(I1),I2=mean(I2),I3=mean(I3),I4=mean(I4),I5=mean(I5)))
  
  
  
  
  ##############################################################################################################
  ######################################## 4 Consecutive negative result #######################################
  ##############################################################################################################
  
  K <- 4
  
  ############################################################ Interval of tests = 1 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 1
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 2 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 2
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 3 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 3
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 4 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 4
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 5 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 5
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  #################### Result ####################
  
  Prob24 <- rbind(Prob24,data.frame(P1=sum(P1)*100/N,P2=sum(P2)*100/N,P3=sum(P3)*100/N,P4=sum(P4)*100/N,P5=sum(P5)*100/N))
  Leng24 <- rbind(Leng24,data.frame(L1=mean(L1),L2=mean(L2),L3=mean(L3),L4=mean(L4),L5=mean(L5)))
  Burd24 <- rbind(Burd24,data.frame(B1=mean(B1),B2=mean(B2),B3=mean(B3),B4=mean(B4),B5=mean(B5)))
  Isol24 <- rbind(Isol24,data.frame(I1=mean(I1),I2=mean(I2),I3=mean(I3),I4=mean(I4),I5=mean(I5)))
  
  
  
  
  ##############################################################################################################
  ######################################## 5 Consecutive negative result #######################################
  ##############################################################################################################
  
  K <- 5
  
  ############################################################ Interval of tests = 1 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 1
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        & df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        & df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 2 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 2
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        & df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        & df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 3 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 3
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        & df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        & df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 4 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 4
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        & df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        & df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 5 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 5
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        & df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        & df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  #################### Result ####################
  
  Prob25 <- rbind(Prob25,data.frame(P1=sum(P1)*100/N,P2=sum(P2)*100/N,P3=sum(P3)*100/N,P4=sum(P4)*100/N,P5=sum(P5)*100/N))
  Leng25 <- rbind(Leng25,data.frame(L1=mean(L1),L2=mean(L2),L3=mean(L3),L4=mean(L4),L5=mean(L5)))
  Burd25 <- rbind(Burd25,data.frame(B1=mean(B1),B2=mean(B2),B3=mean(B3),B4=mean(B4),B5=mean(B5)))
  Isol25 <- rbind(Isol25,data.frame(I1=mean(I1),I2=mean(I2),I3=mean(I3),I4=mean(I4),I5=mean(I5)))
  
  
  df_totalM <- subset(df_totalM, select=-cri)
  
  
  
  
  
  
  
  
  
  
  #########################################################################################################################
  #########################################################################################################################
  #########################################################################################################################
  ###################################### Infectiousness threshold = 10^6.5 copies/ml ######################################
  #########################################################################################################################
  #########################################################################################################################
  #########################################################################################################################
  
  IF <- log10(10^6.5)
  DL <- IF
  
  ################### Time point where the patient loses infectiousness ###################
  cri <- c()
  
  for ( i in 1:TM ) {
    if (df_totalM$V[i]<IF) {
      cri[i] <- df_totalM$time[i]
    }else {
      cri[i] <- Tmax+1
    }
  }
  
  df_totalM   <-cbind(df_totalM,cri)
  df_totalM$V <- as.numeric(df_totalM$V)
  df_totalM$obs <- as.numeric(df_totalM$obs)
  
  
  
  
  ##############################################################################################################
  ######################################## 1 Consecutive negative result #######################################
  ##############################################################################################################
  
  K <- 1
  
  ############################################################ Interval of tests = 1 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 1
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        # & i%%k==1
        #& df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        #& df_totalM$time[i]<df_totalM$time[i+1*k]
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 2 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 2
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & i%%k==1
        #& df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        #& df_totalM$time[i]<df_totalM$time[i+1*k]
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 3 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 3
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & i%%k==1
        #& df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        #& df_totalM$time[i]<df_totalM$time[i+1*k]
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 4 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 4
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & i%%k==1
        #& df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        #& df_totalM$time[i]<df_totalM$time[i+1*k]
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 5 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 5
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & i%%k==1
        #& df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        #& df_totalM$time[i]<df_totalM$time[i+1*k]
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  #################### Result ####################
  
  Prob31 <- rbind(Prob31,data.frame(P1=sum(P1)*100/N,P2=sum(P2)*100/N,P3=sum(P3)*100/N,P4=sum(P4)*100/N,P5=sum(P5)*100/N))
  Leng31 <- rbind(Leng31,data.frame(L1=mean(L1),L2=mean(L2),L3=mean(L3),L4=mean(L4),L5=mean(L5)))
  Burd31 <- rbind(Burd31,data.frame(B1=mean(B1),B2=mean(B2),B3=mean(B3),B4=mean(B4),B5=mean(B5)))
  Isol31 <- rbind(Isol31,data.frame(I1=mean(I1),I2=mean(I2),I3=mean(I3),I4=mean(I4),I5=mean(I5)))
  
  
  
  
  ##############################################################################################################
  ######################################## 2 Consecutive negative result #######################################
  ##############################################################################################################
  
  K <- 2
  
  ############################################################ Interval of tests = 1 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 1
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 2 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 2
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 3 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 3
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 4 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 4
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 5 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 5
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        #& df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        #& df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  #################### Result ####################
  
  Prob32 <- rbind(Prob32,data.frame(P1=sum(P1)*100/N,P2=sum(P2)*100/N,P3=sum(P3)*100/N,P4=sum(P4)*100/N,P5=sum(P5)*100/N))
  Leng32 <- rbind(Leng32,data.frame(L1=mean(L1),L2=mean(L2),L3=mean(L3),L4=mean(L4),L5=mean(L5)))
  Burd32 <- rbind(Burd32,data.frame(B1=mean(B1),B2=mean(B2),B3=mean(B3),B4=mean(B4),B5=mean(B5)))
  Isol32 <- rbind(Isol32,data.frame(I1=mean(I1),I2=mean(I2),I3=mean(I3),I4=mean(I4),I5=mean(I5)))
  
  
  
  
  ##############################################################################################################
  ######################################## 3 Consecutive negative result #######################################
  ##############################################################################################################
  
  K <- 3
  
  ############################################################ Interval of tests = 1 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 1
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 2 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 2
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 3 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 3
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 4 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 4
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 5 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 5
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        #& df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        #& df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  #################### Result ####################
  
  Prob33 <- rbind(Prob33,data.frame(P1=sum(P1)*100/N,P2=sum(P2)*100/N,P3=sum(P3)*100/N,P4=sum(P4)*100/N,P5=sum(P5)*100/N))
  Leng33 <- rbind(Leng33,data.frame(L1=mean(L1),L2=mean(L2),L3=mean(L3),L4=mean(L4),L5=mean(L5)))
  Burd33 <- rbind(Burd33,data.frame(B1=mean(B1),B2=mean(B2),B3=mean(B3),B4=mean(B4),B5=mean(B5)))
  Isol33 <- rbind(Isol33,data.frame(I1=mean(I1),I2=mean(I2),I3=mean(I3),I4=mean(I4),I5=mean(I5)))
  
  
  
  
  ##############################################################################################################
  ######################################## 4 Consecutive negative result #######################################
  ##############################################################################################################
  
  K <- 4
  
  ############################################################ Interval of tests = 1 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 1
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 2 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 2
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 3 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 3
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 4 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 4
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 5 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 5
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        #& df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        #& df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  #################### Result ####################
  
  Prob34 <- rbind(Prob34,data.frame(P1=sum(P1)*100/N,P2=sum(P2)*100/N,P3=sum(P3)*100/N,P4=sum(P4)*100/N,P5=sum(P5)*100/N))
  Leng34 <- rbind(Leng34,data.frame(L1=mean(L1),L2=mean(L2),L3=mean(L3),L4=mean(L4),L5=mean(L5)))
  Burd34 <- rbind(Burd34,data.frame(B1=mean(B1),B2=mean(B2),B3=mean(B3),B4=mean(B4),B5=mean(B5)))
  Isol34 <- rbind(Isol34,data.frame(I1=mean(I1),I2=mean(I2),I3=mean(I3),I4=mean(I4),I5=mean(I5)))
  
  
  
  
  ##############################################################################################################
  ######################################## 5 Consecutive negative result #######################################
  ##############################################################################################################
  
  K <- 5
  
  ############################################################ Interval of tests = 1 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 1
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        & df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        & df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 2 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 2
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        & df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        & df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 3 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 3
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        & df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        & df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 4 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 4
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        & df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        & df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  ############################################################ Interval of tests = 5 day ############################################################
  P <- c(); PP <- c();
  L <- c(); LL <- c();
  B <- c(); BB <- c();
  I <- c(); II <- c();
  
  k <- 5
  
  #### Isolation & Burden
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1*k]<DL
        & df_totalM$obs[i+2*k]<DL
        & df_totalM$obs[i+3*k]<DL
        & df_totalM$obs[i+4*k]<DL
        & df_totalM$time[i]<df_totalM$time[i+1*k]
        & df_totalM$time[i+1*k]<df_totalM$time[i+2*k]
        & df_totalM$time[i+2*k]<df_totalM$time[i+3*k]
        & df_totalM$time[i+3*k]<df_totalM$time[i+4*k]
    ){
      I[i]<-df_totalM$time[i+k*(K-1)]
    }else {
      I[i]<-Tmax+1 ## meaningless
    }
  }
  
  I = data.frame(I)
  df_totalM<-cbind(df_totalM,I)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    II[j]<-min(df_totalM3$I)
    BB[j]<-min(df_totalM3$I) - min(df_totalM3$cri)
    
  }
  
  
  #### Length
  L <- -BB
  for( j in 1:N ) {
    
    if (L[j]>=0){
      LL[j]<-L[j]
    }else {
      LL[j]<-0
    }
    
  }
  
  
  #### Probability
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (df_totalM3$V[which(df_totalM3$time==II[j])] > IF){
      PP[j]<-1
    }else {
      PP[j]<-0
    }
    
  }
  
  assign(paste0('P',k), PP)
  assign(paste0('L',k), LL)
  assign(paste0('B',k), BB)
  assign(paste0('I',k), II)
  
  df_totalM <- subset(df_totalM, select=-I)
  
  
  #################### Result ####################
  
  Prob35 <- rbind(Prob35,data.frame(P1=sum(P1)*100/N,P2=sum(P2)*100/N,P3=sum(P3)*100/N,P4=sum(P4)*100/N,P5=sum(P5)*100/N))
  Leng35 <- rbind(Leng35,data.frame(L1=mean(L1),L2=mean(L2),L3=mean(L3),L4=mean(L4),L5=mean(L5)))
  Burd35 <- rbind(Burd35,data.frame(B1=mean(B1),B2=mean(B2),B3=mean(B3),B4=mean(B4),B5=mean(B5)))
  Isol35 <- rbind(Isol35,data.frame(I1=mean(I1),I2=mean(I2),I3=mean(I3),I4=mean(I4),I5=mean(I5)))
  
  
  df_totalM <- subset(df_totalM, select=-cri)
  
}

#################### Results ####################

setwd("~/capsule/results")

## Infectiousness threshold = 10^6 copies/ml (Default)
## For example, choose a criterion of 5-day interval between tests for ending isolation
## Column name "number" in the below means the number of consecutive negative results for ending isolation

P2 <- data.frame(number=c(1,2,3,4,5),value=t(Prob25))
L2 <- data.frame(number=c(1,2,3,4,5),value=t(Leng25))
B2 <- data.frame(number=c(1,2,3,4,5),value=t(Burd25))
I2 <- data.frame(number=c(1,2,3,4,5),value=t(Isol25))


write_xlsx(P2,"Test_P.xlsx") 
write_xlsx(L2,"Test_L.xlsx") 
write_xlsx(B2,"Test_B.xlsx") 
write_xlsx(I2,"Test_I.xlsx") 

