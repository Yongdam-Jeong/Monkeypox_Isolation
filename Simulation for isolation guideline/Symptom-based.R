library(ggplot2)
library(plyr)
library(deSolve)
library(writexl)
library(readxl)
library(xlsx)


setwd("~/Desktop/")


######################################################################################
######################################################################################
######################################## Model #######################################
######################################################################################
######################################################################################

#################### Setting ####################
Tmin <- 0
Tmax <- 100
step_size <- 1
mtime <- seq(Tmin,Tmax,step_size)

N <- 1000 ## # of samples

IF1 <- log10(10^5.5) ## Infectiousness threshold
IF2 <- log10(10^6.0)
IF3 <- log10(10^6.5)

Lmin <- ceiling(qlnorm(0.025,meanlog=3.14,sdlog=0.41)) ## Duration of lesion presence
Lmax <- floor(qlnorm(0.975,meanlog=3.14,sdlog=0.41))


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


##################################### Simulation #####################################
simulation <- 2

P1 <- matrix(NA, nrow = 21, ncol = simulation)
P2 <- matrix(NA, nrow = 21, ncol = simulation)
P3 <- matrix(NA, nrow = 21, ncol = simulation)

L1 <- matrix(NA, nrow = 21, ncol = simulation)
L2 <- matrix(NA, nrow = 21, ncol = simulation)
L3 <- matrix(NA, nrow = 21, ncol = simulation)

B1 <- matrix(NA, nrow = 21, ncol = simulation)
B2 <- matrix(NA, nrow = 21, ncol = simulation)
B3 <- matrix(NA, nrow = 21, ncol = simulation)


for (s in 1:simulation) {
  
  
  while (TRUE) {
    #################### Parameters ####################
    pop <- read.csv("populationParameters_mpox.txt", row.names = 1)
    
    popdelta <- rlnorm(N, log(pop["delta_pop", "value"]), pop["omega_delta", "value"])
    popv0    <- rlnorm(N, log(pop["V0_pop", "value"]), pop["omega_V0", "value"])
    fitt     <- data.frame(delta=popdelta,V0=popv0)
    
    
    #################### Sampling 1000 COVID-19 patients ####################
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
    
    
    if ( length(which(df_total$V < IF1 & df_total$time == 0))==0 & length(which(df_total$V > IF1 & df_total$time == Tmax))==0 )
      break
    
  }
  

  #################### Sampling timing of lesion clearance (Lesion period)
  LP <- c()
  j <- 1
  while (TRUE) {
    
    # lp <- rgamma(1,shape=Fitt$estimate[1],scale=Fitt$estimate[2])
    lp <- rlnorm(1,meanlog=3.14,sdlog=0.41)
    
    if (lp >= Lmin & lp <= Lmax) {
      LP[j] <- round(lp, digits=0)
      j<-j+1
      if (j>N)
        break
    }
  }
  
  
  ##################################################################################################
  ##################################################################################################
  ######################################## Computing metrics #######################################
  ##################################################################################################
  ##################################################################################################
  
  df_totalM <- df_total
  df_totalM$time <- round(df_totalM$time,digits = 0)
  
  
  #################### Combining data frame with timing of infectiousness loss
  IF11 <- c()
  IF22 <- c()
  IF33 <- c()
  
  for (k in 1:(length(mtime)*N)) {
    
    if (df_totalM$V[k] < IF1) {
      IF11[k] <- 1
    } else {
      IF11[k] <- 0
    }
    
    if (df_totalM$V[k] < IF2) {
      IF22[k] <- 1
    } else {
      IF22[k] <- 0
    }
    
    if (df_totalM$V[k] < IF3) {
      IF33[k] <- 1
    } else {
      IF33[k] <- 0
    }
    
  }
  
  df_totalM$IF1 <- IF11
  df_totalM$IF2 <- IF22
  df_totalM$IF3 <- IF33
  
  
  for (k in 1:21) {
    
    Risk1 <- c()
    Risk2 <- c()
    Risk3 <- c()
    
    Length1 <- c()
    Length2 <- c()
    Length3 <- c()
    
    Burden1 <- c()
    Burden2 <- c()
    Burden3 <- c()
    
    for (i in 1:N) {
      
      df_ind <- subset(df_totalM,ID==i)
      
      lp <- LP[i] -5 + (k-1)
      
      
      ####################################################################################
      ######################################## Probability of prematurely ending isolation
      ####################################################################################
      
      if (df_ind$V[which(df_ind$time==lp)] >= IF1) {
        Risk1[i] <- 1
      } else {
        Risk1[i] <- 0
      }
      
      if (df_ind$V[which(df_ind$time==lp)] >= IF2) {
        Risk2[i] <- 1
      } else {
        Risk2[i] <- 0
      }
      
      if (df_ind$V[which(df_ind$time==lp)] >= IF3) {
        Risk3[i] <- 1
      } else {
        Risk3[i] <- 0
      }
     
      
      #############################################################################################
      ######################################## Length of infectious period after ending isolation
      #############################################################################################
      
      lost1 <- min(df_ind$time[which(df_ind$IF1==1)])
      lost2 <- min(df_ind$time[which(df_ind$IF2==1)])
      lost3 <- min(df_ind$time[which(df_ind$IF3==1)])
      
      if ((lost1-lp) >= 0) {
        Length1[i] <- lost1-lp
      } else {
        Length1[i] <- 0
      }
      
      if ((lost2-lp) >= 0) {
        Length2[i] <- lost2-lp
      } else {
        Length2[i] <- 0
      }
      
      if ((lost3-lp) >= 0) {
        Length3[i] <- lost3-lp
      } else {
        Length3[i] <- 0
      }
      
      
      ####################################################################################
      ######################################## Burden of unnecessarily prolonged isolation
      ####################################################################################
      
      Burden1[i] <- lp - lost1
      Burden2[i] <- lp - lost2
      Burden3[i] <- lp - lost3
      
    }
    
    P1[k,s] <- sum(Risk1)/N*100
    P2[k,s] <- sum(Risk2)/N*100
    P3[k,s] <- sum(Risk3)/N*100
    
    L1[k,s] <- sum(Length1)/N
    L2[k,s] <- sum(Length2)/N
    L3[k,s] <- sum(Length3)/N
    
    B1[k,s] <- sum(Burden1)/N
    B2[k,s] <- sum(Burden2)/N
    B3[k,s] <- sum(Burden3)/N
    
  }
  
}
