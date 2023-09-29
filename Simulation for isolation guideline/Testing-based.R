library(ggplot2)
library(plyr)
library(deSolve)
library(writexl)
library(readxl)
library(xlsx)
library(pheatmap)
library(RColorBrewer)
library(patchwork)
library(desplot)

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

Prob11 <- data.frame(); Prob21 <- data.frame(); Prob31 <- data.frame()
Prob12 <- data.frame(); Prob22 <- data.frame(); Prob32 <- data.frame()
Prob13 <- data.frame(); Prob23 <- data.frame(); Prob33 <- data.frame()
Prob14 <- data.frame(); Prob24 <- data.frame(); Prob34 <- data.frame()
Prob15 <- data.frame(); Prob25 <- data.frame(); Prob35 <- data.frame()

Leng11 <- data.frame(); Leng21 <- data.frame(); Leng31 <- data.frame()
Leng12 <- data.frame(); Leng22 <- data.frame(); Leng32 <- data.frame()
Leng13 <- data.frame(); Leng23 <- data.frame(); Leng33 <- data.frame()
Leng14 <- data.frame(); Leng24 <- data.frame(); Leng34 <- data.frame()
Leng15 <- data.frame(); Leng25 <- data.frame(); Leng35 <- data.frame()

Burd11 <- data.frame(); Burd21 <- data.frame(); Burd31 <- data.frame()
Burd12 <- data.frame(); Burd22 <- data.frame(); Burd32 <- data.frame()
Burd13 <- data.frame(); Burd23 <- data.frame(); Burd33 <- data.frame()
Burd14 <- data.frame(); Burd24 <- data.frame(); Burd34 <- data.frame()
Burd15 <- data.frame(); Burd25 <- data.frame(); Burd35 <- data.frame()

Isol11 <- data.frame(); Isol21 <- data.frame(); Isol31 <- data.frame()
Isol12 <- data.frame(); Isol22 <- data.frame(); Isol32 <- data.frame()
Isol13 <- data.frame(); Isol23 <- data.frame(); Isol33 <- data.frame()
Isol14 <- data.frame(); Isol24 <- data.frame(); Isol34 <- data.frame()
Isol15 <- data.frame(); Isol25 <- data.frame(); Isol35 <- data.frame()


simulation <- 2


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


write_xlsx(Prob11,"Prob11.xlsx"); write_xlsx(Prob21,"Prob21.xlsx"); write_xlsx(Prob31,"Prob31.xlsx")
write_xlsx(Prob12,"Prob12.xlsx"); write_xlsx(Prob22,"Prob22.xlsx"); write_xlsx(Prob32,"Prob32.xlsx")
write_xlsx(Prob13,"Prob13.xlsx"); write_xlsx(Prob23,"Prob23.xlsx"); write_xlsx(Prob33,"Prob33.xlsx")
write_xlsx(Prob14,"Prob14.xlsx"); write_xlsx(Prob24,"Prob24.xlsx"); write_xlsx(Prob34,"Prob34.xlsx")
write_xlsx(Prob15,"Prob15.xlsx"); write_xlsx(Prob25,"Prob25.xlsx"); write_xlsx(Prob35,"Prob35.xlsx")

write_xlsx(Leng11,"Leng11.xlsx"); write_xlsx(Leng21,"Leng21.xlsx"); write_xlsx(Leng31,"Leng31.xlsx")
write_xlsx(Leng12,"Leng12.xlsx"); write_xlsx(Leng22,"Leng22.xlsx"); write_xlsx(Leng32,"Leng32.xlsx")
write_xlsx(Leng13,"Leng13.xlsx"); write_xlsx(Leng23,"Leng23.xlsx"); write_xlsx(Leng33,"Leng33.xlsx")
write_xlsx(Leng14,"Leng14.xlsx"); write_xlsx(Leng24,"Leng24.xlsx"); write_xlsx(Leng34,"Leng34.xlsx")
write_xlsx(Leng15,"Leng15.xlsx"); write_xlsx(Leng25,"Leng25.xlsx"); write_xlsx(Leng35,"Leng35.xlsx")

write_xlsx(Burd11,"Burd11.xlsx"); write_xlsx(Burd21,"Burd21.xlsx"); write_xlsx(Burd31,"Burd31.xlsx")
write_xlsx(Burd12,"Burd12.xlsx"); write_xlsx(Burd22,"Burd22.xlsx"); write_xlsx(Burd32,"Burd32.xlsx")
write_xlsx(Burd13,"Burd13.xlsx"); write_xlsx(Burd23,"Burd23.xlsx"); write_xlsx(Burd33,"Burd33.xlsx")
write_xlsx(Burd14,"Burd14.xlsx"); write_xlsx(Burd24,"Burd24.xlsx"); write_xlsx(Burd34,"Burd34.xlsx")
write_xlsx(Burd15,"Burd15.xlsx"); write_xlsx(Burd25,"Burd25.xlsx"); write_xlsx(Burd35,"Burd35.xlsx")

write_xlsx(Isol11,"Isol11.xlsx"); write_xlsx(Isol21,"Isol21.xlsx"); write_xlsx(Isol31,"Isol31.xlsx")
write_xlsx(Isol12,"Isol12.xlsx"); write_xlsx(Isol22,"Isol22.xlsx"); write_xlsx(Isol32,"Isol32.xlsx")
write_xlsx(Isol13,"Isol13.xlsx"); write_xlsx(Isol23,"Isol23.xlsx"); write_xlsx(Isol33,"Isol33.xlsx")
write_xlsx(Isol14,"Isol14.xlsx"); write_xlsx(Isol24,"Isol24.xlsx"); write_xlsx(Isol34,"Isol34.xlsx")
write_xlsx(Isol15,"Isol15.xlsx"); write_xlsx(Isol25,"Isol25.xlsx"); write_xlsx(Isol35,"Isol35.xlsx")










#############################################################################################################################################################
#############################################################################################################################################################
########################################################################## Figures ##########################################################################
#############################################################################################################################################################
#############################################################################################################################################################


Prob11 <- read.xlsx("Prob11.xlsx", sheetIndex = 1); Prob21 <- read.xlsx("Prob21.xlsx", sheetIndex = 1); Prob31 <- read.xlsx("Prob31.xlsx", sheetIndex = 1)
Prob12 <- read.xlsx("Prob12.xlsx", sheetIndex = 1); Prob22 <- read.xlsx("Prob22.xlsx", sheetIndex = 1); Prob32 <- read.xlsx("Prob32.xlsx", sheetIndex = 1)
Prob13 <- read.xlsx("Prob13.xlsx", sheetIndex = 1); Prob23 <- read.xlsx("Prob23.xlsx", sheetIndex = 1); Prob33 <- read.xlsx("Prob33.xlsx", sheetIndex = 1)
Prob14 <- read.xlsx("Prob14.xlsx", sheetIndex = 1); Prob24 <- read.xlsx("Prob24.xlsx", sheetIndex = 1); Prob34 <- read.xlsx("Prob34.xlsx", sheetIndex = 1)
Prob15 <- read.xlsx("Prob15.xlsx", sheetIndex = 1); Prob25 <- read.xlsx("Prob25.xlsx", sheetIndex = 1); Prob35 <- read.xlsx("Prob35.xlsx", sheetIndex = 1)

Leng11 <- read.xlsx("Leng11.xlsx", sheetIndex = 1); Leng21 <- read.xlsx("Leng21.xlsx", sheetIndex = 1); Leng31 <- read.xlsx("Leng31.xlsx", sheetIndex = 1)
Leng12 <- read.xlsx("Leng12.xlsx", sheetIndex = 1); Leng22 <- read.xlsx("Leng22.xlsx", sheetIndex = 1); Leng32 <- read.xlsx("Leng32.xlsx", sheetIndex = 1)
Leng13 <- read.xlsx("Leng13.xlsx", sheetIndex = 1); Leng23 <- read.xlsx("Leng23.xlsx", sheetIndex = 1); Leng33 <- read.xlsx("Leng33.xlsx", sheetIndex = 1)
Leng14 <- read.xlsx("Leng14.xlsx", sheetIndex = 1); Leng24 <- read.xlsx("Leng24.xlsx", sheetIndex = 1); Leng34 <- read.xlsx("Leng34.xlsx", sheetIndex = 1)
Leng15 <- read.xlsx("Leng15.xlsx", sheetIndex = 1); Leng25 <- read.xlsx("Leng25.xlsx", sheetIndex = 1); Leng35 <- read.xlsx("Leng35.xlsx", sheetIndex = 1)

Burd11 <- read.xlsx("Burd11.xlsx", sheetIndex = 1); Burd21 <- read.xlsx("Burd21.xlsx", sheetIndex = 1); Burd31 <- read.xlsx("Burd31.xlsx", sheetIndex = 1)
Burd12 <- read.xlsx("Burd12.xlsx", sheetIndex = 1); Burd22 <- read.xlsx("Burd22.xlsx", sheetIndex = 1); Burd32 <- read.xlsx("Burd32.xlsx", sheetIndex = 1)
Burd13 <- read.xlsx("Burd13.xlsx", sheetIndex = 1); Burd23 <- read.xlsx("Burd23.xlsx", sheetIndex = 1); Burd33 <- read.xlsx("Burd33.xlsx", sheetIndex = 1)
Burd14 <- read.xlsx("Burd14.xlsx", sheetIndex = 1); Burd24 <- read.xlsx("Burd24.xlsx", sheetIndex = 1); Burd34 <- read.xlsx("Burd34.xlsx", sheetIndex = 1)
Burd15 <- read.xlsx("Burd15.xlsx", sheetIndex = 1); Burd25 <- read.xlsx("Burd25.xlsx", sheetIndex = 1); Burd35 <- read.xlsx("Burd35.xlsx", sheetIndex = 1)

Isol11 <- read.xlsx("Isol11.xlsx", sheetIndex = 1); Isol21 <- read.xlsx("Isol21.xlsx", sheetIndex = 1); Isol31 <- read.xlsx("Isol31.xlsx", sheetIndex = 1)
Isol12 <- read.xlsx("Isol12.xlsx", sheetIndex = 1); Isol22 <- read.xlsx("Isol22.xlsx", sheetIndex = 1); Isol32 <- read.xlsx("Isol32.xlsx", sheetIndex = 1)
Isol13 <- read.xlsx("Isol13.xlsx", sheetIndex = 1); Isol23 <- read.xlsx("Isol23.xlsx", sheetIndex = 1); Isol33 <- read.xlsx("Isol33.xlsx", sheetIndex = 1)
Isol14 <- read.xlsx("Isol14.xlsx", sheetIndex = 1); Isol24 <- read.xlsx("Isol24.xlsx", sheetIndex = 1); Isol34 <- read.xlsx("Isol34.xlsx", sheetIndex = 1)
Isol15 <- read.xlsx("Isol15.xlsx", sheetIndex = 1); Isol25 <- read.xlsx("Isol25.xlsx", sheetIndex = 1); Isol35 <- read.xlsx("Isol35.xlsx", sheetIndex = 1)


####################################### Infectiousness1
Probability1 <- data.frame(results=rep(c(1,2,3,4,5),each=5), 
                           interval=rep(c(1,2,3,4,5),times=5),
                           value=c(mean(Prob11$P1),mean(Prob11$P2),mean(Prob11$P3),mean(Prob11$P4),mean(Prob11$P5),
                                   mean(Prob12$P1),mean(Prob12$P2),mean(Prob12$P3),mean(Prob12$P4),mean(Prob12$P5),
                                   mean(Prob13$P1),mean(Prob13$P2),mean(Prob13$P3),mean(Prob13$P4),mean(Prob13$P5),
                                   mean(Prob14$P1),mean(Prob14$P2),mean(Prob14$P3),mean(Prob14$P4),mean(Prob14$P5),
                                   mean(Prob15$P1),mean(Prob15$P2),mean(Prob15$P3),mean(Prob15$P4),mean(Prob15$P5)))

Length1      <- data.frame(results=rep(c(1,2,3,4,5),each=5), 
                           interval=rep(c(1,2,3,4,5),times=5),
                           value=c(mean(Leng11$L1),mean(Leng11$L2),mean(Leng11$L3),mean(Leng11$L4),mean(Leng11$L5),
                                   mean(Leng12$L1),mean(Leng12$L2),mean(Leng12$L3),mean(Leng12$L4),mean(Leng12$L5),
                                   mean(Leng13$L1),mean(Leng13$L2),mean(Leng13$L3),mean(Leng13$L4),mean(Leng13$L5),
                                   mean(Leng14$L1),mean(Leng14$L2),mean(Leng14$L3),mean(Leng14$L4),mean(Leng14$L5),
                                   mean(Leng15$L1),mean(Leng15$L2),mean(Leng15$L3),mean(Leng15$L4),mean(Leng15$L5)))

Burden1      <- data.frame(results=rep(c(1,2,3,4,5),each=5), 
                           interval=rep(c(1,2,3,4,5),times=5),
                           value=c(mean(Burd11$B1),mean(Burd11$B2),mean(Burd11$B3),mean(Burd11$B4),mean(Burd11$B5),
                                   mean(Burd12$B1),mean(Burd12$B2),mean(Burd12$B3),mean(Burd12$B4),mean(Burd12$B5),
                                   mean(Burd13$B1),mean(Burd13$B2),mean(Burd13$B3),mean(Burd13$B4),mean(Burd13$B5),
                                   mean(Burd14$B1),mean(Burd14$B2),mean(Burd14$B3),mean(Burd14$B4),mean(Burd14$B5),
                                   mean(Burd15$B1),mean(Burd15$B2),mean(Burd15$B3),mean(Burd15$B4),mean(Burd15$B5)))

Isolation1   <- data.frame(results=rep(c(1,2,3,4,5),each=5), 
                           interval=rep(c(1,2,3,4,5),times=5),
                           value=c(mean(Isol11$I1),mean(Isol11$I2),mean(Isol11$I3),mean(Isol11$I4),mean(Isol11$I5),
                                   mean(Isol12$I1),mean(Isol12$I2),mean(Isol12$I3),mean(Isol12$I4),mean(Isol12$I5),
                                   mean(Isol13$I1),mean(Isol13$I2),mean(Isol13$I3),mean(Isol13$I4),mean(Isol13$I5),
                                   mean(Isol14$I1),mean(Isol14$I2),mean(Isol14$I3),mean(Isol14$I4),mean(Isol14$I5),
                                   mean(Isol15$I1),mean(Isol15$I2),mean(Isol15$I3),mean(Isol15$I4),mean(Isol15$I5)))


####################################### Infectiousness2
Probability2 <- data.frame(results=rep(c(1,2,3,4,5),each=5), 
                           interval=rep(c(1,2,3,4,5),times=5),
                           value=c(mean(Prob21$P1),mean(Prob21$P2),mean(Prob21$P3),mean(Prob21$P4),mean(Prob21$P5),
                                   mean(Prob22$P1),mean(Prob22$P2),mean(Prob22$P3),mean(Prob22$P4),mean(Prob22$P5),
                                   mean(Prob23$P1),mean(Prob23$P2),mean(Prob23$P3),mean(Prob23$P4),mean(Prob23$P5),
                                   mean(Prob24$P1),mean(Prob24$P2),mean(Prob24$P3),mean(Prob24$P4),mean(Prob24$P5),
                                   mean(Prob25$P1),mean(Prob25$P2),mean(Prob25$P3),mean(Prob25$P4),mean(Prob25$P5)))

Length2      <- data.frame(results=rep(c(1,2,3,4,5),each=5), 
                           interval=rep(c(1,2,3,4,5),times=5),
                           value=c(mean(Leng21$L1),mean(Leng21$L2),mean(Leng21$L3),mean(Leng21$L4),mean(Leng21$L5),
                                   mean(Leng22$L1),mean(Leng22$L2),mean(Leng22$L3),mean(Leng22$L4),mean(Leng22$L5),
                                   mean(Leng23$L1),mean(Leng23$L2),mean(Leng23$L3),mean(Leng23$L4),mean(Leng23$L5),
                                   mean(Leng24$L1),mean(Leng24$L2),mean(Leng24$L3),mean(Leng24$L4),mean(Leng24$L5),
                                   mean(Leng25$L1),mean(Leng25$L2),mean(Leng25$L3),mean(Leng25$L4),mean(Leng25$L5)))

Burden2      <- data.frame(results=rep(c(1,2,3,4,5),each=5), 
                           interval=rep(c(1,2,3,4,5),times=5),
                           value=c(mean(Burd21$B1),mean(Burd21$B2),mean(Burd21$B3),mean(Burd21$B4),mean(Burd21$B5),
                                   mean(Burd22$B1),mean(Burd22$B2),mean(Burd22$B3),mean(Burd22$B4),mean(Burd22$B5),
                                   mean(Burd23$B1),mean(Burd23$B2),mean(Burd23$B3),mean(Burd23$B4),mean(Burd23$B5),
                                   mean(Burd24$B1),mean(Burd24$B2),mean(Burd24$B3),mean(Burd24$B4),mean(Burd24$B5),
                                   mean(Burd25$B1),mean(Burd25$B2),mean(Burd25$B3),mean(Burd25$B4),mean(Burd25$B5)))

Isolation2   <- data.frame(results=rep(c(1,2,3,4,5),each=5), 
                           interval=rep(c(1,2,3,4,5),times=5),
                           value=c(mean(Isol21$I1),mean(Isol21$I2),mean(Isol21$I3),mean(Isol21$I4),mean(Isol21$I5),
                                   mean(Isol22$I1),mean(Isol22$I2),mean(Isol22$I3),mean(Isol22$I4),mean(Isol22$I5),
                                   mean(Isol23$I1),mean(Isol23$I2),mean(Isol23$I3),mean(Isol23$I4),mean(Isol23$I5),
                                   mean(Isol24$I1),mean(Isol24$I2),mean(Isol24$I3),mean(Isol24$I4),mean(Isol24$I5),
                                   mean(Isol25$I1),mean(Isol25$I2),mean(Isol25$I3),mean(Isol25$I4),mean(Isol25$I5)))


####################################### Infectiousness3
Probability3 <- data.frame(results=rep(c(1,2,3,4,5),each=5), 
                           interval=rep(c(1,2,3,4,5),times=5),
                           value=c(mean(Prob31$P1),mean(Prob31$P2),mean(Prob31$P3),mean(Prob31$P4),mean(Prob31$P5),
                                   mean(Prob32$P1),mean(Prob32$P2),mean(Prob32$P3),mean(Prob32$P4),mean(Prob32$P5),
                                   mean(Prob33$P1),mean(Prob33$P2),mean(Prob33$P3),mean(Prob33$P4),mean(Prob33$P5),
                                   mean(Prob34$P1),mean(Prob34$P2),mean(Prob34$P3),mean(Prob34$P4),mean(Prob34$P5),
                                   mean(Prob35$P1),mean(Prob35$P2),mean(Prob35$P3),mean(Prob35$P4),mean(Prob35$P5)))

Length3      <- data.frame(results=rep(c(1,2,3,4,5),each=5), 
                           interval=rep(c(1,2,3,4,5),times=5),
                           value=c(mean(Leng31$L1),mean(Leng31$L2),mean(Leng31$L3),mean(Leng31$L4),mean(Leng31$L5),
                                   mean(Leng32$L1),mean(Leng32$L2),mean(Leng32$L3),mean(Leng32$L4),mean(Leng32$L5),
                                   mean(Leng33$L1),mean(Leng33$L2),mean(Leng33$L3),mean(Leng33$L4),mean(Leng33$L5),
                                   mean(Leng34$L1),mean(Leng34$L2),mean(Leng34$L3),mean(Leng34$L4),mean(Leng34$L5),
                                   mean(Leng35$L1),mean(Leng35$L2),mean(Leng35$L3),mean(Leng35$L4),mean(Leng35$L5)))

Burden3      <- data.frame(results=rep(c(1,2,3,4,5),each=5), 
                           interval=rep(c(1,2,3,4,5),times=5),
                           value=c(mean(Burd31$B1),mean(Burd31$B2),mean(Burd31$B3),mean(Burd31$B4),mean(Burd31$B5),
                                   mean(Burd32$B1),mean(Burd32$B2),mean(Burd32$B3),mean(Burd32$B4),mean(Burd32$B5),
                                   mean(Burd33$B1),mean(Burd33$B2),mean(Burd33$B3),mean(Burd33$B4),mean(Burd33$B5),
                                   mean(Burd34$B1),mean(Burd34$B2),mean(Burd34$B3),mean(Burd34$B4),mean(Burd34$B5),
                                   mean(Burd35$B1),mean(Burd35$B2),mean(Burd35$B3),mean(Burd35$B4),mean(Burd35$B5)))

Isolation3   <- data.frame(results=rep(c(1,2,3,4,5),each=5), 
                           interval=rep(c(1,2,3,4,5),times=5),
                           value=c(mean(Isol31$I1),mean(Isol31$I2),mean(Isol31$I3),mean(Isol31$I4),mean(Isol31$I5),
                                   mean(Isol32$I1),mean(Isol32$I2),mean(Isol32$I3),mean(Isol32$I4),mean(Isol32$I5),
                                   mean(Isol33$I1),mean(Isol33$I2),mean(Isol33$I3),mean(Isol33$I4),mean(Isol33$I5),
                                   mean(Isol34$I1),mean(Isol34$I2),mean(Isol34$I3),mean(Isol34$I4),mean(Isol34$I5),
                                   mean(Isol35$I1),mean(Isol35$I2),mean(Isol35$I3),mean(Isol35$I4),mean(Isol35$I5)))



t_col <- function(color, percent = 50, name = NULL) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  invisible(t.col)
}

col1 <- "#1C115C"
col2 <- "#5D1A66"
col3 <- "#DF67DA"

col11 <- t_col(col1)
col22 <- t_col(col2)
col33 <- t_col(col3)

xlabel <- "Interval between tests (days)"
ylabel <- "Consecutive negative results"
plt1<-ggplot(Probability1, aes(x=interval, y=results, fill=value)) +
  geom_tile() +
  xlab(xlabel) +
  ylab(ylabel) +
  scale_fill_gradient2(name = "(%)", low = "gray97", high = col1, mid = col11, midpoint = 50, breaks = c(0,25,50,75,100), limits = c(0,100))+
  theme(plot.title = element_text(size = 18, hjust = 0.5, vjust = -2),
        legend.title = element_text(size = 15,vjust = 4),
        legend.text = element_text(size = 15),
        legend.key.height = unit(1.72,'cm'),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 16, vjust = 2),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18, vjust = 1),
        axis.title.y = element_text(size = 18, vjust = 3.5)) +
  coord_fixed()


plt2<-ggplot(Probability2, aes(x=interval, y=results, fill=value)) +
  geom_tile() +
  xlab(xlabel) +
  ylab(ylabel) +
  scale_fill_gradient2(name = "(%)", low = "gray97", high = col2, mid = col22, midpoint = 50, breaks = c(0,25,50,75,100), limits = c(0,100))+
  theme(plot.title = element_text(size = 18, hjust = 0.5, vjust = -2),
        legend.title = element_text(size = 15,vjust = 4),
        legend.text = element_text(size = 15),
        legend.key.height = unit(1.72,'cm'),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 16, vjust = 2),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18, vjust = 1),
        axis.title.y = element_text(size = 18, vjust = 3.5)) +
  coord_fixed()


plt3<-ggplot(Probability3, aes(x=interval, y=results, fill=value)) +
  geom_tile() +
  xlab(xlabel) +
  ylab(ylabel) +
  scale_fill_gradient2(name = "(%)", low = "gray97", high = col3, mid = col33, midpoint = 50, breaks = c(0,25,50,75,100), limits = c(0,100))+
  theme(plot.title = element_text(size = 18, hjust = 0.5, vjust = -2),
        legend.title = element_text(size = 15,vjust = 4),
        legend.text = element_text(size = 15),
        legend.key.height = unit(1.72,'cm'),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 16, vjust = 2),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18, vjust = 1),
        axis.title.y = element_text(size = 18, vjust = 3.5)) +
  coord_fixed()


plt4<-ggplot(Length1, aes(x=interval, y=results, fill=value)) +
  geom_tile() +
  xlab(xlabel) +
  ylab(ylabel) +
  scale_fill_gradient2(name = "(days)", low = "gray97", high = col1, mid = col11, midpoint = 5, breaks = c(0,2,4,6,8,10), limits = c(0,10))+
  theme(plot.title = element_text(size = 18, hjust = 0.5, vjust = -2),
        legend.title = element_text(size = 15,vjust = 4),
        legend.text = element_text(size = 15),
        legend.key.height = unit(1.72,'cm'),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 16, vjust = 2),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18, vjust = 1),
        axis.title.y = element_text(size = 18, vjust = 3.5)) +
  coord_fixed()


plt5<-ggplot(Length2, aes(x=interval, y=results, fill=value)) +
  geom_tile() +
  xlab(xlabel) +
  ylab(ylabel) +
  scale_fill_gradient2(name = "(days)", low = "gray97", high = col2, mid = col22, midpoint = 5, breaks = c(0,2,4,6,8,10), limits = c(0,10))+
  theme(plot.title = element_text(size = 18, hjust = 0.5, vjust = -2),
        legend.title = element_text(size = 15,vjust = 4),
        legend.text = element_text(size = 15),
        legend.key.height = unit(1.72,'cm'),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 16, vjust = 2),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18, vjust = 1),
        axis.title.y = element_text(size = 18, vjust = 3.5)) +
  coord_fixed()


plt6<-ggplot(Length3, aes(x=interval, y=results, fill=value)) +
  geom_tile() +
  xlab(xlabel) +
  ylab(ylabel) +
  scale_fill_gradient2(name = "(days)", low = "gray97", high = col3, mid = col33, midpoint = 5, breaks = c(0,2,4,6,8,10), limits = c(0,10))+
  theme(plot.title = element_text(size = 18, hjust = 0.5, vjust = -2),
        legend.title = element_text(size = 15,vjust = 4),
        legend.text = element_text(size = 15),
        legend.key.height = unit(1.72,'cm'),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 16, vjust = 2),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18, vjust = 1),
        axis.title.y = element_text(size = 18, vjust = 3.5)) +
  coord_fixed()


plt7<-ggplot(Burden1, aes(x=interval, y=results, fill=value)) +
  geom_tile() +
  xlab(xlabel) +
  ylab(ylabel) +
  scale_fill_gradient2(name = "(days)", low = "gray45", high = col1, mid = "gray97", midpoint = 0, limits = c(-10,20))+
  theme(plot.title = element_text(size = 18, hjust = 0.5, vjust = -2),
        legend.title = element_text(size = 15,vjust = 4),
        legend.text = element_text(size = 15),
        legend.key.height = unit(1.72,'cm'),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 16, vjust = 2),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18, vjust = 1),
        axis.title.y = element_text(size = 18, vjust = 3.5)) +
  coord_fixed()


plt8<-ggplot(Burden2, aes(x=interval, y=results, fill=value)) +
  geom_tile() +
  xlab(xlabel) +
  ylab(ylabel) +
  scale_fill_gradient2(name = "(days)", low = "gray45", high = col2, mid = "gray97", midpoint = 0, limits = c(-10,20))+
  theme(plot.title = element_text(size = 18, hjust = 0.5, vjust = -2),
        legend.title = element_text(size = 15,vjust = 4),
        legend.text = element_text(size = 15),
        legend.key.height = unit(1.72,'cm'),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 16, vjust = 2),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18, vjust = 1),
        axis.title.y = element_text(size = 18, vjust = 3.5)) +
  coord_fixed()


plt9<-ggplot(Burden3, aes(x=interval, y=results, fill=value)) +
  geom_tile() +
  xlab(xlabel) +
  ylab(ylabel) +
  scale_fill_gradient2(name = "(days)", low = "gray45", high = col3, mid = "gray97", midpoint = 0, limits = c(-10,20))+
  theme(plot.title = element_text(size = 18, hjust = 0.5, vjust = -2),
        legend.title = element_text(size = 15,vjust = 4),
        legend.text = element_text(size = 15),
        legend.key.height = unit(1.72,'cm'),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 16, vjust = 2),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18, vjust = 1),
        axis.title.y = element_text(size = 18, vjust = 3.5)) +
  coord_fixed()


plt10<-ggplot(Isolation1, aes(x=interval, y=results, fill=value)) +
  geom_tile() +
  xlab(xlabel) +
  ylab(ylabel) +
  scale_fill_gradient2(name = "(days)", low = "gray97", high = col1, mid = col11, midpoint = 20, limits = c(0,40))+
  theme(plot.title = element_text(size = 18, hjust = 0.5, vjust = -2),
        legend.title = element_text(size = 15,vjust = 4),
        legend.text = element_text(size = 15),
        legend.key.height = unit(1.72,'cm'),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 16, vjust = 2),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18, vjust = 1),
        axis.title.y = element_text(size = 18, vjust = 3.5)) +
  coord_fixed()


plt11<-ggplot(Isolation2, aes(x=interval, y=results, fill=value)) +
  geom_tile() +
  xlab(xlabel) +
  ylab(ylabel) +
  scale_fill_gradient2(name = "(days)", low = "gray97", high = col2, mid = col22, midpoint = 20, limits = c(0,40))+
  theme(plot.title = element_text(size = 18, hjust = 0.5, vjust = -2),
        legend.title = element_text(size = 15,vjust = 4),
        legend.text = element_text(size = 15),
        legend.key.height = unit(1.72,'cm'),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 16, vjust = 2),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18, vjust = 1),
        axis.title.y = element_text(size = 18, vjust = 3.5)) +
  coord_fixed()


plt12<-ggplot(Isolation3, aes(x=interval, y=results, fill=value)) +
  geom_tile() +
  # geom_tileborder(aes(group=1, grp=Optimal)) +
  xlab(xlabel) +
  ylab(ylabel) +
  scale_fill_gradient2(name = "(days)", low = "gray97", high = col3, mid = col33, midpoint = 20, limits = c(0,40))+
  theme(plot.title = element_text(size = 18, hjust = 0.5, vjust = -2),
        legend.title = element_text(size = 15,vjust = 4),
        legend.text = element_text(size = 15),
        legend.key.height = unit(1.72,'cm'),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 16, vjust = 2),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18, vjust = 1),
        axis.title.y = element_text(size = 18, vjust = 3.5)) +
  coord_fixed()


sfn <- "Personalized.pdf"
# pdf(sfn,height = (6)*3, width = (10)*3)
pdf(sfn,height = (6)*4, width = (10)*3)
theme_set(theme_classic(base_size = 32, base_family = "Helvetica"))
plt <- wrap_elements(plt1)  + wrap_elements(plt2)  + wrap_elements(plt3)  + 
  wrap_elements(plt4)  + wrap_elements(plt5)  + wrap_elements(plt6)  +
  wrap_elements(plt7)  + wrap_elements(plt8)  + wrap_elements(plt9)  +
  wrap_elements(plt10) + wrap_elements(plt11) + wrap_elements(plt12) +
  # plot_layout(nrow=3, heights = c(1,1,1))
  plot_layout(nrow=4, heights = c(1,1,1,1))
print( plt )
dev.off()


