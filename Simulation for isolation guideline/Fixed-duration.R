library(ggplot2)
library(plyr)
library(deSolve)
library(writexl)
library(readxl)
library(xlsx)


setwd("~/Desktop/")


#################### Setting ####################
Tmin <- 0
Tmax <- 100
step_size <- 1
mtime <- seq(Tmin,Tmax,step_size)

N <- 1000 ## # of samples

IF1 <- log10(10^5.5) ## Infectiousness threshold
IF2 <- log10(10^6.0)
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
simulation <- 2

P1 <- matrix(NA, nrow = length(mtime), ncol = simulation)
P2 <- matrix(NA, nrow = length(mtime), ncol = simulation)
P3 <- matrix(NA, nrow = length(mtime), ncol = simulation)

L1 <- matrix(NA, nrow = length(mtime), ncol = simulation)
L2 <- matrix(NA, nrow = length(mtime), ncol = simulation)
L3 <- matrix(NA, nrow = length(mtime), ncol = simulation)

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


P1 <- data.frame(P1); write.xlsx(P1,"One_P1.xlsx"); 
P2 <- data.frame(P2); write.xlsx(P2,"One_P2.xlsx"); 
P3 <- data.frame(P3); write.xlsx(P3,"One_P3.xlsx"); 

L1 <- data.frame(L1); write.xlsx(L1,"One_L1.xlsx"); 
L2 <- data.frame(L2); write.xlsx(L2,"One_L2.xlsx"); 
L3 <- data.frame(L3); write.xlsx(L3,"One_L3.xlsx"); 

B1 <- data.frame(B1); write.xlsx(B1,"One_B1.xlsx"); 
B2 <- data.frame(B2); write.xlsx(B2,"One_B2.xlsx"); 
B3 <- data.frame(B3); write.xlsx(B3,"One_B3.xlsx"); 


#######################################################################################
#######################################################################################
######################################## Figure #######################################
#######################################################################################
#######################################################################################

P1 <- as.matrix(read.xlsx("One_P1.xlsx", sheetIndex = 1)); P1 <- P1[,-1]
P2 <- as.matrix(read.xlsx("One_P2.xlsx", sheetIndex = 1)); P2 <- P2[,-1]
P3 <- as.matrix(read.xlsx("One_P3.xlsx", sheetIndex = 1)); P3 <- P3[,-1]

L1 <- as.matrix(read.xlsx("One_L1.xlsx", sheetIndex = 1)); L1 <- L1[,-1]
L2 <- as.matrix(read.xlsx("One_L2.xlsx", sheetIndex = 1)); L2 <- L2[,-1]
L3 <- as.matrix(read.xlsx("One_L3.xlsx", sheetIndex = 1)); L3 <- L3[,-1]

B1 <- as.matrix(read.xlsx("One_B1.xlsx", sheetIndex = 1)); B1 <- B1[,-1]
B2 <- as.matrix(read.xlsx("One_B2.xlsx", sheetIndex = 1)); B2 <- B2[,-1]
B3 <- as.matrix(read.xlsx("One_B3.xlsx", sheetIndex = 1)); B3 <- B3[,-1]


P11 <- matrix(NA, nrow = length(mtime), ncol = 5)
P22 <- matrix(NA, nrow = length(mtime), ncol = 5)
P33 <- matrix(NA, nrow = length(mtime), ncol = 5)

L11 <- matrix(NA, nrow = length(mtime), ncol = 5)
L22 <- matrix(NA, nrow = length(mtime), ncol = 5)
L33 <- matrix(NA, nrow = length(mtime), ncol = 5)

B11 <- matrix(NA, nrow = length(mtime), ncol = 5)
B22 <- matrix(NA, nrow = length(mtime), ncol = 5)
B33 <- matrix(NA, nrow = length(mtime), ncol = 5)

for (k in 1:length(mtime)) {
  
  P11[k,1] <- k-1; L11[k,1] <- k-1; B11[k,1] <- k-1;
  P22[k,1] <- k-1; L22[k,1] <- k-1; B22[k,1] <- k-1;
  P33[k,1] <- k-1; L33[k,1] <- k-1; B33[k,1] <- k-1;
  
  P11[k,2] <- quantile(as.numeric(P1[k,]),0.025); L11[k,2] <- quantile(as.numeric(L1[k,]),0.025); B11[k,2] <- quantile(as.numeric(B1[k,]),0.025);
  P22[k,2] <- quantile(as.numeric(P2[k,]),0.025); L22[k,2] <- quantile(as.numeric(L2[k,]),0.025); B22[k,2] <- quantile(as.numeric(B2[k,]),0.025);
  P33[k,2] <- quantile(as.numeric(P3[k,]),0.025); L33[k,2] <- quantile(as.numeric(L3[k,]),0.025); B33[k,2] <- quantile(as.numeric(B3[k,]),0.025);
  
  P11[k,3] <- mean(as.numeric(P1[k,])); L11[k,3] <- mean(as.numeric(L1[k,])); B11[k,3] <- mean(as.numeric(B1[k,]));
  P22[k,3] <- mean(as.numeric(P2[k,])); L22[k,3] <- mean(as.numeric(L2[k,])); B22[k,3] <- mean(as.numeric(B2[k,]));
  P33[k,3] <- mean(as.numeric(P3[k,])); L33[k,3] <- mean(as.numeric(L3[k,])); B33[k,3] <- mean(as.numeric(B3[k,]));
  
  P11[k,4] <- quantile(as.numeric(P1[k,]),0.975); L11[k,4] <- quantile(as.numeric(L1[k,]),0.975); B11[k,4] <- quantile(as.numeric(B1[k,]),0.975);
  P22[k,4] <- quantile(as.numeric(P2[k,]),0.975); L22[k,4] <- quantile(as.numeric(L2[k,]),0.975); B22[k,4] <- quantile(as.numeric(B2[k,]),0.975);
  P33[k,4] <- quantile(as.numeric(P3[k,]),0.975); L33[k,4] <- quantile(as.numeric(L3[k,]),0.975); B33[k,4] <- quantile(as.numeric(B3[k,]),0.975);
  
  P11[k,5] <- median(as.numeric(P1[k,])); L11[k,5] <- median(as.numeric(L1[k,])); B11[k,5] <- median(as.numeric(B1[k,]));
  P22[k,5] <- median(as.numeric(P2[k,])); L22[k,5] <- median(as.numeric(L2[k,])); B22[k,5] <- median(as.numeric(B2[k,]));
  P33[k,5] <- median(as.numeric(P3[k,])); L33[k,5] <- median(as.numeric(L3[k,])); B33[k,5] <- median(as.numeric(B3[k,]));
  
}

P11 <- data.frame(P11); colnames(P11) <- c("day","min","mean","max","median")
P22 <- data.frame(P22); colnames(P22) <- c("day","min","mean","max","median")
P33 <- data.frame(P33); colnames(P33) <- c("day","min","mean","max","median")

L11 <- data.frame(L11); colnames(L11) <- c("day","min","mean","max","median")
L22 <- data.frame(L22); colnames(L22) <- c("day","min","mean","max","median")
L33 <- data.frame(L33); colnames(L33) <- c("day","min","mean","max","median")

B11 <- data.frame(B11); colnames(B11) <- c("day","min","mean","max","median")
B22 <- data.frame(B22); colnames(B22) <- c("day","min","mean","max","median")
B33 <- data.frame(B33); colnames(B33) <- c("day","min","mean","max","median")

P111 <- subset(P11,P11$mean<=5); P111 <- min(P111$day) 
L111 <- subset(L11,L11$mean<=1); L111 <- min(L111$day)
OneI1 <- max(P111,L111); OneB1 <- B11$mean[which(B11$day==OneI1)]

P222 <- subset(P22,P22$mean<=5); P222 <- min(P222$day) 
L222 <- subset(L22,L22$mean<=1); L222 <- min(L222$day)
OneI2 <- max(P222,L222); OneB2 <- B22$mean[which(B22$day==OneI2)]

P333 <- subset(P33,P33$mean<=5); P333 <- min(P333$day) 
L333 <- subset(L33,L33$mean<=1); L333 <- min(L333$day)
OneI3 <- max(P333,L333); OneB3 <- B33$mean[which(B33$day==OneI3)]


col1 <- "#1C115C"
col2 <- "#5D1A66"
col3 <- "#DF67DA"


sfn <- "One_risk.pdf"
pdf(sfn,width=14,height=10)
theme_set(theme_classic(base_size = 32, base_family = "Helvetica"))
xlabel <- "Isolation period (days)"
ylabel <- "Probability of prematurely ending isolation (%)"
plt1 <- ggplot() +

  geom_ribbon(data=P11,aes(x=day,ymin=min,ymax=max), fill=col1, alpha=0.2) +
  geom_path(data=P11,aes(x=day,y=mean), color=col1, lwd=1.2) +

  geom_ribbon(data=P22,aes(x=day,ymin=min,ymax=max), fill=col2, alpha=0.2) +
  geom_path(data=P22,aes(x=day,y=mean), color=col2, lwd=1.2) +

  geom_ribbon(data=P33,aes(x=day,ymin=min,ymax=max), fill=col3, alpha=0.2) +
  geom_path(data=P33,aes(x=day,y=mean), color=col3, lwd=1.2) +
  
  xlab(xlabel) +
  ylab(ylabel) +

  scale_x_continuous(breaks=c(seq(0,30,by=5)),limits=c(0,30)) +
  scale_y_continuous(breaks=seq(0,100,by=20),limits=c(0,100)) +
  
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(linetype = "dotted",color = "gray90"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position='none',
        axis.title.y = element_text(family="Helvetica"),
        axis.title.x = element_text(family="Helvetica"))
print( plt1 )
dev.off()


sfn <- "One_length.pdf"
pdf(sfn,width=14,height=10)
theme_set(theme_classic(base_size = 32, base_family = "Helvetica"))
xlabel <- "Isolation period (days)"
ylabel <- "Infectious period after ending isolation (days)"
plt2 <- ggplot() +
  
  geom_ribbon(data=L11,aes(x=day,ymin=min,ymax=max), fill=col1, alpha=0.2) +
  geom_path(data=L11,aes(x=day,y=mean), color=col1, lwd=1.2) +
  
  geom_ribbon(data=L22,aes(x=day,ymin=min,ymax=max), fill=col2, alpha=0.2) +
  geom_path(data=L22,aes(x=day,y=mean), color=col2, lwd=1.2) +
  
  geom_ribbon(data=L33,aes(x=day,ymin=min,ymax=max), fill=col3, alpha=0.2) +
  geom_path(data=L33,aes(x=day,y=mean), color=col3, lwd=1.2) +
  
  xlab(xlabel) +
  ylab(ylabel) +
  
  scale_x_continuous(breaks=c(seq(0,30,by=5)),limits=c(0,30)) +
  scale_y_continuous(breaks=seq(0,20,by=4),limits=c(0,20)) +
  
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(linetype = "dotted",color = "gray90"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position='none',
        axis.title.y = element_text(family="Helvetica"),
        axis.title.x = element_text(family="Helvetica"))
print( plt2 )
dev.off()


sfn <- "One_burden.pdf"
pdf(sfn,width=14,height=10)
theme_set(theme_classic(base_size = 32, base_family = "Helvetica"))
xlabel <- "Isolation period (days)"
ylabel <- "Unnecessarily prolonged isolation period (days)"
plt3 <- ggplot() +
  
  geom_ribbon(data=B11,aes(x=day,ymin=min,ymax=max), fill=col1, alpha=0.2) +
  geom_path(data=B11,aes(x=day,y=mean), color=col1, lwd=1.2) +

  geom_ribbon(data=B22,aes(x=day,ymin=min,ymax=max), fill=col2, alpha=0.2) +
  geom_path(data=B22,aes(x=day,y=mean), color=col2, lwd=1.2) +

  geom_ribbon(data=B33,aes(x=day,ymin=min,ymax=max), fill=col3, alpha=0.2) +
  geom_path(data=B33,aes(x=day,y=mean), color=col3, lwd=1.2) +

  xlab(xlabel) +
  ylab(ylabel) +
  
  scale_x_continuous(breaks=c(seq(0,30,by=5)),limits=c(0,30)) +
  scale_y_continuous(breaks=seq(-20,30,by=10),limits=c(-20,30)) +
  
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(linetype = "dotted",color = "gray90"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position='none',
        axis.title.y = element_text(family="Helvetica"),
        axis.title.x = element_text(family="Helvetica"))
print( plt3 )
dev.off()


sfn <- "One_size_fits_all.pdf"
pdf(sfn,width=15*3,height=11)
theme_set(theme_classic(base_size = 32, base_family = "Helvetica"))
plt <- wrap_elements(plt1) + wrap_elements(plt2) + wrap_elements(plt3) +
  plot_layout(ncol=3, widths = c(1,1,1))
print( plt )
dev.off()
