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

IF1 <- log10(10^5.5) ## Infectiousness threshold
IF2 <- log10(10^6.0) ## Default value
IF3 <- log10(10^6.5)


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
simulation <- 1


#### Risk of prematurely ending isolation (Probability => P)
P1 <- matrix(NA, nrow = length(mtime), ncol = simulation)
P2 <- matrix(NA, nrow = length(mtime), ncol = simulation)
P3 <- matrix(NA, nrow = length(mtime), ncol = simulation)

#### Infectious period after ending isolation (Length => L)
L1 <- matrix(NA, nrow = length(mtime), ncol = simulation)
L2 <- matrix(NA, nrow = length(mtime), ncol = simulation)
L3 <- matrix(NA, nrow = length(mtime), ncol = simulation)

#### Unnecessarily prolonged isolation period (Burden => B)
B1 <- matrix(NA, nrow = length(mtime), ncol = simulation)
B2 <- matrix(NA, nrow = length(mtime), ncol = simulation)
B3 <- matrix(NA, nrow = length(mtime), ncol = simulation)


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


##################################################################################################
##################################################################################################
######################################## Computing metrics #######################################
##################################################################################################
##################################################################################################

df_totalM <- df_total
df_totalM$time <- round(df_totalM$time,digits = 0)

####################################################################################
######################################## Probability of prematurely ending isolation
####################################################################################

for (k in 1:length(mtime)) {
  
  P1[k,s] <- length(which(df_totalM$V > IF1 & df_totalM$time == (k-1)))/N*100
  P2[k,s] <- length(which(df_totalM$V > IF2 & df_totalM$time == (k-1)))/N*100
  P3[k,s] <- length(which(df_totalM$V > IF3 & df_totalM$time == (k-1)))/N*100
  
}

#############################################################################################
######################################## Length of infectious period after ending isolation
#############################################################################################

cri1 <- c()
cri2 <- c()
cri3 <- c()

for ( i in 1:(length(mtime)*N) ) {
  
  if (df_totalM$V[i] < IF1){
    cri1[i] <- df_totalM$time[i]
  }else {
    cri1[i] <- Inf ## meaningless
  }
  
  if (df_totalM$V[i] < IF2){
    cri2[i] <- df_totalM$time[i]
  }else {
    cri2[i] <- Inf ## meaningless
  }
  
  if (df_totalM$V[i] < IF3){
    cri3[i] <- df_totalM$time[i]
  }else {
    cri3[i] <- Inf ## meaningless
  }

}

df_totalM <- cbind(df_totalM,cri1,cri2,cri3)
# df_totalM <- subset(df_totalM,select=-c(cri1,cri2,cri3))

cri11 <- c()
cri22 <- c()
cri33 <- c()

for ( i in 1:N ) {
  
  df_totalM1 <- subset(df_totalM,df_totalM$ID==i)
  
  cri11[i]<- min(df_totalM1$cri1)
  cri22[i]<- min(df_totalM1$cri2)
  cri33[i]<- min(df_totalM1$cri3)
  
}

for (k in 1:length(mtime)) {
  
  cri111 <- cri11 - (k-1)
  cri222 <- cri22 - (k-1)
  cri333 <- cri33 - (k-1)
  
  cri111[cri111<0] <- 0
  cri222[cri222<0] <- 0
  cri333[cri333<0] <- 0
  
  L1[k,s] <- mean(cri111)
  L2[k,s] <- mean(cri222)
  L3[k,s] <- mean(cri333)
  
}


####################################################################################
######################################## Burden of unnecessarily prolonged isolation
####################################################################################

for (k in 1:length(mtime)) {
  
  B1[k,s] <- (k-1) - mean(cri11)
  B2[k,s] <- (k-1) - mean(cri22)
  B3[k,s] <- (k-1) - mean(cri33)
  
}

}




#################### Results ####################

setwd("~/capsule/results")

## Infectiousness threshold = 10^6 copies/ml (Default)

P2 <- data.frame(time=mtime,value=P2)
L2 <- data.frame(time=mtime,value=L2)
B2 <- data.frame(time=mtime,value=B2)

write_xlsx(P2,"Fix_P.xlsx") 
write_xlsx(L2,"Fix_L.xlsx") 
write_xlsx(B2,"Fix_B.xlsx") 


