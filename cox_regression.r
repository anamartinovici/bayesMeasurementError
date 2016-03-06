library(rjags)
#library(lattice)
library(MASS)

myburnin <- 1000
mysample <- 5000

#small simulation study
n<- 1000 #sample size
n1<- 100 #size of calibration sample
xzcorr <- 0.25
kappa<-2 #weibull shape parameter

gammamean <- c(0,0)
gammaprec <- 0.0001*diag(2)
betamean <- c(0,0)
betaprec <- diag(c(0.72,0.72))

reliability <- 0.5
b1 <- 0.5
lambda<-0.00164 #weibull scale parameter

cov <- mvrnorm(n, mu=c(0,0), Sigma=array(c(1,xzcorr,xzcorr,1),dim=c(2,2)))
x <- cov[,1]
z <- cov[,2]

u<-runif(n,0,1)
t.event<-((-log(u)/lambda)*(kappa^(lambda+1))*exp(-(b1*x+b1*z)))^(1/kappa)
d<-ifelse(t.event<10,1,0)
t<-ifelse(d==1,t.event,10)
mean(d)

sigma_u_sq <- 1/reliability - 1
w1 <- x+rnorm(n, sd=sigma_u_sq^0.5)
w2 <- x+rnorm(n, sd=sigma_u_sq^0.5)
w2[(n1+1):n] <- NA

mydata <- data.frame(t,d,z,w1,w2)

inits1 <- list(beta=c(0,0), gamma=c(0,0), taux=1, tauu=1)
inits2 <- list(beta=c(1,-1), gamma=c(5,1), taux=1.5, tauu=0.5)
inits3 <- list(beta=c(2,2), gamma=c(-5,1), taux=2, tauu=0.2)
inits4 <- list(beta=c(-1,-1.5), gamma=c(5,5), taux=0.5, tauu=2)
inits5 <- list(beta=c(-2,-0.5), gamma=c(-5,-5), taux=0.25, tauu=4)

nevent<-length(mydata$t[mydata$d==1])
tevent<-sort(mydata$t[mydata$d==1])

jags.data<-list(N = n, N1=n1, T = nevent, eps = 1.0E-10,c=0.001,r=0.01,
                obs.t = t,
                fail = d,
                w1 = w1,
                w2=w2,
                z=z,
                t = c(tevent,10),
                "betamean"=betamean, 
                "betaprec"=betaprec, "gammamean"=gammamean, "gammaprec"=gammaprec,
                "tauu_alpha"=0.5, "tauu_beta"=0.5, "taux_alpha"=0.5, "taux_beta"=0.5)

jags.inits <- list(inits1,inits2,inits3,inits4,inits5)

jags.params <- c("gamma", "tauu", "taux", "beta")

jagsmodel <- jags.model(data=jags.data, file="cox_reg.bug",n.chains=5, inits=jags.inits)
burnin <- coda.samples(jagsmodel, variable.names=c("beta"), n.iter=myburnin)
mainsample <- coda.samples(jagsmodel, variable.names=c("beta"), n.iter=mysample)
