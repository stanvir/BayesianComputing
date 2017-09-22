#Clear Workspace
rm(list=ls())

#setwd("C:/Users/SHAMS/Google Drive/Bayesian Volume Prediction")

library(splines)
library(rjags)
library(ggplot2)

#Raw Data Import
Traffic<- read.csv("input _3months.csv")
# Setting Up Initial Values
#add rush hour flag and avg. of adjacent tmc travel rate
Traffic$flag<-ifelse(Traffic$hour %in% c(7:10,16:19),1,0) #rush hour flag

Segment <- as.integer(as.numeric(Traffic[,10]))
NSegmenttype <-nlevels(Traffic[,10]) 

Traffic$Cnt<- as.numeric(Traffic[,2])   #station
Cnt <- Traffic$Cnt

Y<-Traffic[,4]/Traffic[,12]  #Normalized volume
Y <- as.integer(Y)
X<-Traffic[,5]  #travel rate
X<-X/max(X)       #normalizing travel rate

X1<-Traffic[,7]  #travel rate @t-1
X1<-X1/max(X1)    #normalizing

Flag<-Traffic$flag  #rush hour flag

Time<-Traffic[,3]/max(Traffic[,3])  #hour of day

Adj<-(Traffic[,11]+Traffic[,12])/2
Adj <- Adj/max(Adj)
#Sgmnt<-Traffic[,11]

AADT<-Traffic[,6]

##Need to set out variables for test data
##Yp, Xp, AADTp, Adjp, Flagp, X1p

##need to define: spline for travel rate
J <- 3
B<-matrix(0,432,3)    #Cubic(3) function for 6 stations * 72 times = 432. Question: should it be 18*24?
for(i in 1:6){
        B[c(((i-1)*72+1):(i*72)),] <- bs(X[c(((i-1)*72+1):(i*72))],J)
}

B

##need to define: spline for hour
T <- 5
C <- bs(Time,T)
dim(C)
matplot(Time,C,type="l",xlab="Time",ylab="Basis function, Bj(X)",cex.lab=1.5,cex.axis=1.5,lwd=2)


n<-length(Y)
n


# segment type <- 4
# time_period <- 72

#Define Model

model <- "model{

#likelihood
for (i in 1:n){
Y[i] ~ dpois(AADT[i]*lam[i])
lam[i] <-  exp(
alpha + 
alpha_segment[Segment[i]]+
#alpha_Time[] # either random effect or spline
inprod((beta2[]+ beta22[Segment[i],]),X[i,])+
(beta3*X1[i])+
(beta_Flag*Flag[i])+
(beta5* Adj[i]) +
inprod(beta6[],C[i,])   #beta6 is the parameter for hour
)

}

# # Prediction
# for(i in 1:np){
# Yp[i]  ~ dpois(AADTp[i]*lam)
# lam <-  exp(
# alpha + 
# alpha_segment[Segmentp[i]]+
# inprod(beta2[],Xp[i,])+
# (beta3*X1p[i])+
# (beta_Flag*Flagp[i])+
# (beta5* Adjp[i]) 
# )
# }

#Prior - Intercept   
alpha ~ dnorm(0,0.01)


sigma_s ~ dunif(0, 100) # standard deviation of random effect (variance between segments)
tau_s <- 1/(sigma_s * sigma_s) #precision for segments
#Prior- Segment Type
for (i in 1:NSegmenttype){
alpha_segment[i] ~ dnorm(0, tau_s)
}

for (i in 1:NSegmenttype){
for (j in 1:J){
beta22[i,j] ~ ddexp(0,inv.var4) 
}
}

#Prior Splines for travel rate
for(j in 1:J){
beta2[j] ~ ddexp(0,inv.var1) 
}	


#Prior (t-1)travel rate 
beta3 ~ ddexp(0,inv.var2)


#Prior Flags
beta_Flag  ~ ddexp(0,inv.varF)


#Prior adjacent
beta5 ~ ddexp(0,inv.var3)


#Prior Splines for hour
for(k in 1:5){
beta6[k] ~ ddexp(0,inv.var5) 
}	

#Hyperparameters
inv.var1 ~ dgamma(0.01, 0.01)
inv.var2 ~ dgamma(0.01, 0.01)
inv.varF ~ dgamma(0.01, 0.01)
inv.var3 ~ dgamma(0.01, 0.01)
inv.var4 ~ dgamma(0.01, 0.01)
inv.var5 ~ dgamma(0.01, 0.01)

#Predictions
for (i in 1:n){
pred[i]<-AADT[i]*lam[i]
}
}"

#########################################################################
# Calculate DIC for the model.
data <- list(AADT=AADT, Y=Y, J=J, X=B, X1=X1, n=n,
             Flag=Flag, Adj=Adj, Segment= Segment,
             NSegmenttype= NSegmenttype,Time=Time,C=C)

model <- jags.model(textConnection(model), n.chains = 3,
                    data= data)


update(model, 6000); # Burnin for 6000 samples
dic <- dic.samples(model, n.iter=10000,
                   variable.names=c("pred"))
dic

#########################################################################

# Run the model, summarize convergence and estimates

params <- c("alpha","beta2","beta3",
            "beta5","beta_Flag","pred","lam")

params2<-c("pred")

samp <- coda.samples(model, n.iter=10000,
                     variable.names= params2)

#summary(samp)  # this is too long output


## e.g. to get the mean and SD of alpha


meanpred<-as.data.frame(summary(samp)[[1]][,1])
plotdata<-as.data.frame(cbind(meanpred[,1],Y))



colnames(plotdata)<-c("predicted","observed")
str(plotdata)
ggplot(plotdata,aes(observed,predicted))+
  geom_point()+
  geom_abline(slope=1,intercept = 0)+
  ggtitle('with time as a variable')

R2 <- 1 - (sum((plotdata$observed-plotdata$predicted )^2)/sum((plotdata$observed-mean(plotdata$observed))^2))
R2


#for beta2- basis function parameters
summary(samp)[[1]][7:24,1:2]

## and to get its 2.5 and 97.5 percentiles
summary(samp)[[2]][1:6,c(1,5)]


# Convergence diagnostics
pdf("project_mod1_20160424_trace.pdf"); plot(samp); dev.off() # Trace plots
pdf("project_mod1_20160424_GR.pdf"); gelman.plot(samp); dev.off() # G-R plots
pdf("project_mod1_20160424_ACF.pdf"); autocorr.plot(samp); dev.off() # ACF plots

sink("project_mod1_20160424_sampsize.txt")
(sampsize <- effectiveSize(samp)) # Effective sample sizes
sink()


#########################################################################

# Check convergence and size of lambda itself.

update(model, 5000); # Burnin for 5000 samples
samp.lda <- coda.samples(model, n.iter=10000,
                         variable.names=c("lam"))

# Convergence diagnostics
pdf("project_mod1_20160424_ldatrace.pdf"); plot(samp.lda); dev.off() # Trace plots
pdf("project_mod1_20160424_ldaGR.pdf"); gelman.plot(samp.lda); dev.off() # G-R plots
pdf("project_mod1_20160424_ldaACF.pdf"); autocorr.plot(samp.lda); dev.off() # ACF plots

# Plot the estimated lambdas vs. time by station.
post.lda <- summary(samp.lda) # Summary of posteriors
dat.mod1 <- data.frame(station=Traffic$itsStationId,
                       hour=Traffic$Hour,
                       lda.mean=post.lda$statistics[,1])
head(dat.mod1)

library(ggplot2)
plot.lda <- ggplot(dat.mod1, aes(x=hour, y=lda.mean)) +
        geom_line(aes(color=station)) #+ geom_point(aes(color=station))
x11(); plot(plot.lda)

dat.mod1a <- cbind(dat.mod1,AADT)
plot.lda.v.aadt <- ggplot(dat.mod1a) +
        geom_point(aes(x=AADT,y=lda.mean, color=station))
x11(); plot(plot.lda.v.aadt)


sum <- summary(samp)

q<-sum$quantiles
st<-sum$statistics

st

# data_list<-list(n=n,Y=Y,Cnt=Cnt,X=B,J=J,X1=X1,Flag=Flag,Vol=Hist_Vol,AADT=AADT,Time=Hist_Time,Confounders=Confounders)
# model <- jags.model(textConnection(model),data_list,n.chains = 3)
# tick <- proc.time()[3]
# Sys.time()
# 
# update(model, 5000, progress.bar="none"); # Burnin for 5000 samples
# tock <- proc.time()[3]
# tock-tick
# Sys.time()
# 
# 
# tick <- proc.time()[3]
# samp <- coda.samples(model, 
#                      variable.names=c("pred"),
#                      n.iter=10000, progress.bar="none")
# 
# tock <- proc.time()[3]
# tock-tick
# Sys.time()
# effectiveSize(samp)
# 
# sum <- summary(samp)
# 
# q<-sum$quantiles
# st<-sum$statistics
# 
# st