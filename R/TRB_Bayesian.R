#Clear Workspace
rm(list=ls())

rm(list = ls())
#Raw Data Import
Traffic<- read.csv("input.csv")
# Setting Up Initial Values
#add rush hour flag
Traffic$flag<-ifelse(Traffic$hour %in% c(7:10,16:19),1,0)

str(Traffic)
Cnt<- Traffic[,3]
Y<-Traffic[,5]
X<-Traffic[,6]
X<-X/max(X)
X1<-Traffic[,8]
Flag<-Traffic[,15]
X1<-X1/max(X1)
Time<-as.factor(Traffic[,4])
Adj<-Traffic[,14]
Sgmnt<-Traffic[,11]
AADT<-Traffic[,7]
##need to define:

library(splines)

J <- 3
B<-matrix(0,144,3)
for(i in 1:6)
{
        B[c(((i-1)*24+1):(i*24)),] <- bs(X[c(((i-1)*24+1):(i*24))],J)
        
}
##Q: What is bs?
B

n<-length(Y)
n

library(rjags)


#Define Model

model <- "model{
#likelihood
for (i in 1:n){
Y[i] ~ dpois(AADT[i]*lam[i])
lam[i] <-  exp(
alpha[Cnt[i]] + 
inprod(beta2[Cnt[i],],X[i,])+
(beta3[Cnt[i]]*X1[i])+
(beta_Flag[Cnt[i]]*Flag[i])+
(beta5[Cnt[i]]*Adj[i]) 
)
}
#Prior - Intercept   
for(i in 1:6){
alpha[i] ~ dnorm(0,0.01)
}
#Prior Splines 
for(i in 1:6){
for(j in 1:J){
beta2[i,j] ~ ddexp(0,inv.var1) 
}	
}
# #Prior (t-1)travel rate 
for(i in 1:6){
beta3[i] ~ ddexp(0,inv.var2)
}
# #Prior Flags
for(i in 1:6){
beta_Flag[i]  ~ ddexp(0,inv.varF)
}
#Prior adjacent
for(i in 1:6){
beta5[i] ~ ddexp(0,inv.var3)
}

#Hyperparameters
inv.var1 ~ dgamma(0.01, 0.01)
inv.var2 ~ dgamma(0.01, 0.01)
inv.varF ~ dgamma(0.01, 0.01)
inv.var3 ~ dgamma(0.01, 0.01)
#Predictions
for (i in 1:n){
pred[i]<-AADT[i]*lam[i]
}
}"
#########################################################################
# Calculate DIC for the model.
model <- jags.model(textConnection(model), n.chains = 3,
                    list(AADT=AADT, Y=Y, J=J, X=B, X1=X1, n=n,
                         Flag=Flag, Adj=Adj,Cnt=Cnt))

###model is showing error after running till this part###

update(model, 6000); # Burnin for 6000 samples
dic <- dic.samples(model, n.iter=10000,
                   variable.names=c("pred"))
dic

#########################################################################

# Run the model, summarize convergence and estimates
model <- jags.model(textConnection(model), n.chains = 3,
                    list(AADT=AADT, Y=Y, J=J, X=B, X1=X1, n=n, 
                         Flag=Flag, Adj=Adj,Cnt=Cnt))
update(model, 6000); # Burnin for 6000 samples
samp <- coda.samples(model, n.iter=10000,
                     variable.names=c("alpha","beta2","beta3",
                                      "beta5","beta_Flag","pred","lam"))

# Convergence diagnostics
pdf("project_mod1_20160424_trace.pdf"); plot(samp); dev.off() # Trace plots
pdf("project_mod1_20160424_GR.pdf"); gelman.plot(samp); dev.off() # G-R plots
pdf("project_mod1_20160424_ACF.pdf"); autocorr.plot(samp); dev.off() # ACF plots

sink("project_mod1_20160424_sampsize.txt")
(sampsize <- effectiveSize(samp)) # Effective sample sizes
sink()


#########################################################################

# Check convergence and size of lambda itself.
model.lda <- jags.model(textConnection(model.string), n.chains = 3,
                        list(AADT=AADT, Y=Y, J=J, X=B, X1=X1, n=n, 
                             Flag=Flag, Adj=Adj,Cnt=Cnt))

update(model.lda, 5000); # Burnin for 5000 samples
samp.lda <- coda.samples(model.lda, n.iter=10000,
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