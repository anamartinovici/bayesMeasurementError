mydata <- read.csv("nhanes.csv", header=TRUE)

library(lme4)
library(rjags)
library(lattice)
library(MASS)
library(survival)

#naive CCA
ccaData <- subset(mydata, (is.na(sbp1)==FALSE) & (is.na(smoke)==FALSE))
n <- dim(ccaData)[1]
naiveweibull <- survreg(Surv(t,d)~sbp1+sex+age+smoke+diabetes,dist="weibull", data=ccaData)
summary(naiveweibull)
#log hazard ratios are
-coef(naiveweibull)/naiveweibull$scale
#for CIs, you could apply delta method to vcov(naiveweibull), or I just used Stata, which uses 
#parametrization for Weibull regression which gives log hazard ratios directly

#equivalent Bayes analysis
#construct vectors to model as censored weibull
set.seed(87893)
nChains <- 5
t <- ccaData$t
t[ccaData$d==0] <- NA
c <- ccaData$t
iscensored <- 1-ccaData$d

jags.data<-list(n = n, t=t, c=c, is.censored=iscensored,
                sbp1 = ccaData$sbp1, sex=ccaData$sex, age=ccaData$age, smoke=ccaData$smoke, diabetes=ccaData$diabetes)

#JAGS treats the event times of those subjects who were censored as missing data
#we therefore need to specify sensible initial values for these
jags.inits <- vector("list",nChains)
tInits <- array(0, dim=c(nChains,n))
logscaleinits <- array(0, dim=c(nChains, 6))
rinits <- array(0, dim=c(nChains))

for (i in 1:nChains) {
  tInits[i,] <- rep(NA,n)
  tInits[i,iscensored==1] <- ccaData$t[iscensored==1]+0.5*runif(sum(iscensored))
  logscaleinits[i,] <- mvrnorm(1, coef(naiveweibull), 2*vcov(naiveweibull)[1:6,1:6])
  newlogscale <- rnorm(1, log(naiveweibull$scale), 2*vcov(naiveweibull)[7,7]^0.5)
  rinits[i] <- 1/exp(newlogscale)
  jags.inits[[i]] <- list(t=tInits[i,], log.scale=logscaleinits[i,], r=rinits[i])
}

jagsmodel1 <- jags.model(data=jags.data, file="weibull_naive_cca.bug",n.chains=nChains, inits=jags.inits)
burnin1 <- coda.samples(jagsmodel1, variable.names=c("beta","r"), n.iter=5000, thin=1)
gelman.diag(burnin1)
summary(burnin1)
mainsample1 <- coda.samples(jagsmodel1, variable.names=c("beta","r"), n.iter=5000, thin=1)
gelman.diag(mainsample1)
summary(mainsample1)

#regression calibration on complete cases
rc <- function(data, index) {
  locdata <- data[index,]
  n <- length(locdata$sbp1)
  locdata$id <- 1:n
  locdata$meas <- 1
  stackeddata <- rbind(locdata, locdata)
  stackeddata$sbp1[(n+1):(2*n)] <- locdata$sbp2
  stackeddata$meas[(n+1):(2*n)] <- 2
  rc_mod <- lmer(sbp1~1+meas+sex+age+smoke+diabetes+(1|id), data=stackeddata)
  locdata$x_blup <- fitted(rc_mod)[1:n]
  rc_out_mod <- survreg(Surv(t,d)~x_blup+sex+age+smoke+diabetes,dist="weibull", data=locdata)
  -coef(rc_out_mod)/rc_out_mod$scale
}

rc(ccaData, 1:n)
library(boot)

set.seed(6723431)
boot <- boot(ccaData, rc, R=2000)
for (i in 1:6) {
  print(boot.ci(boot, index=i, type="perc"))
}


#Bayes analysis, allowing for measurement error but still CCA
set.seed(87893)
gammamean <- rep(0,5)
gammaprec <- 0.0001*diag(5)

#fit measurement error model to get some overdispersed initial values
wideccaData <- ccaData
wideccaData$id <- 1:n
wideccaData$meas <- 1
longData <- rbind(wideccaData, wideccaData)
longData$sbp1[(n+1):(2*n)] <- ccaData$sbp2
longData$meas[(n+1):(2*n)] <- 2
rc_mod <- lmer(sbp1~1+sex+age+smoke+diabetes+factor(meas)+(1|id), data=longData)

gammaInits <- array(0, dim=c(nChains,5))
shiftInits <- array(0, dim=c(nChains))
tauxInits <- array(0, dim=nChains)
tauuInits <- array(0, dim=nChains)
for (i in 1:nChains) {
  draw <- mvrnorm(1, fixef(rc_mod), 2*vcov(rc_mod))
  gammaInits[i,] <- draw[1:5]
  shiftInits[i] <- draw[6]
  tauxInits[i] <- ((6-i)/2.5)/(VarCorr(rc_mod)$id[1]^2)
  tauuInits[i] <- (i/2.5)/(sigma(rc_mod)^2)
  jags.inits[[i]] <- list(t=tInits[i,], log.scale=logscaleinits[i,], r=rinits[i], gamma=gammaInits[i,], shift=shiftInits[i], taux=tauxInits[i], tauu=tauuInits[i])
}

jags.data<-list(n = n, t=t, c=c, is.censored=iscensored,
                sbp1 = ccaData$sbp1, sbp2 = ccaData$sbp2, sex=ccaData$sex, age=ccaData$age, smoke=ccaData$smoke, diabetes=ccaData$diabetes,
                "gammamean"=gammamean, "gammaprec"=gammaprec,"taux_alpha"=0.5, "taux_beta"=0.5,"tauu_alpha"=0.5, "tauu_beta"=0.5)

jagsmodel2 <- jags.model(data=jags.data, file="weibull_adj_cca.bug",n.chains=nChains, inits=jags.inits)
burnin2 <- coda.samples(jagsmodel2, variable.names=c("beta","r","sigmau","sigmax"), n.iter=5000, thin=1)
gelman.diag(burnin2)
summary(burnin2)
mainsample2 <- coda.samples(jagsmodel2, variable.names=c("beta","r","sigmau","sigmax"), n.iter=5000, thin=1)
gelman.diag(mainsample2)
summary(mainsample2)


#now make use of full sample, handling missing data in smoking
set.seed(87893)
alphamean <- rep(0,4)
alphaprec <- 0.0001*diag(4)

n <- dim(mydata)[1]
t <- mydata$t
t[mydata$d==0] <- NA
c <- mydata$t
iscensored <- 1-mydata$d

#set up initial values
jags.inits <- vector("list",nChains)
tInits <- array(0, dim=c(nChains,n))

smokeModel <- glm(smoke~sex+age+diabetes, mydata, family="binomial")

for (i in 1:nChains) {
  tInits[i,] <- rep(NA,n)
  tInits[i,iscensored==1] <- mydata$t[iscensored==1]+0.5*runif(sum(iscensored))
  alphaInit <- mvrnorm(1, mu=coef(smokeModel), Sigma=2*vcov(smokeModel))
  jags.inits[[i]] <- list(t=tInits[i,], log.scale=logscaleinits[i,], r=rinits[i], gamma=gammaInits[i,], shift=shiftInits[i], 
                          taux=tauxInits[i], tauu=tauuInits[i], alpha=alphaInit)
}

jags.data<-list(n = n, t=t, c=c, is.censored=iscensored ,
                sbp1 = mydata$sbp1, sbp2 = mydata$sbp2, sex=mydata$sex, age=mydata$age, smoke=mydata$smoke, diabetes=mydata$diabetes,
                "alphamean"=alphamean, "alphaprec"=alphaprec, "gammamean"=gammamean, "gammaprec"=gammaprec,
                "taux_alpha"=0.5, "taux_beta"=0.5,"tauu_alpha"=0.5, "tauu_beta"=0.5)

jagsmodel3 <- jags.model(data=jags.data, file="weibull_adj_fullsample.bug",n.chains=nChains, inits=jags.inits)
burnin3 <- coda.samples(jagsmodel3, variable.names=c("beta","r","sigmau","sigmax","shift"), n.iter=5000)
gelman.diag(burnin3)
mainsample3 <- coda.samples(jagsmodel3, variable.names=c("beta","r","alpha","sigmau","sigmax","shift"), n.iter=5000)
gelman.diag(mainsample3)
summary(mainsample3)
