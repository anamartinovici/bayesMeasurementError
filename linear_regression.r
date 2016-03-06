library(rjags)
#library(lattice)
library(MASS)

myburnin <- 1000
mysample <- 5000

#small simulation study
n <- 1000
n1 <- 100
xzcorr <- 0.25

nSims <- 10
nParams <- 3

inits1 <- list(beta=c(0,0,0), gamma=c(0,0), taux=1, tauu=1, tauy=0.1)
inits2 <- list(beta=c(0,1,-1), gamma=c(5,1), taux=1.5, tauu=0.5, tauy=5)
inits3 <- list(beta=c(0,2,2), gamma=c(-5,1), taux=2, tauu=0.2, tauy=0.6)
inits4 <- list(beta=c(0,-1,-1.5), gamma=c(5,5), taux=0.5, tauu=2, tauy=0.4)
inits5 <- list(beta=c(0,-2,-0.5), gamma=c(-5,-5), taux=0.25, tauu=4, tauy=1.5)

jags.inits <- list(inits1,inits2,inits3,inits4,inits5)

reliability <- 0.5
r2 <- 0.1
resvar = (2+2*xzcorr-2*r2-2*xzcorr*r2)/r2

cov <- mvrnorm(n, mu=c(0,0), Sigma=array(c(1,xzcorr,xzcorr,1),dim=c(2,2)))
x <- cov[,1]
z <- cov[,2]

y <- x+z+rnorm(n,sd=resvar^0.5)
sigma_u_sq <- 1/reliability - 1
w1 <- x+rnorm(n, sd=sigma_u_sq^0.5)
w2 <- x+rnorm(n, sd=sigma_u_sq^0.5)
w2[(n1+1):n] <- NA

mydata <- data.frame(y,z,w1,w2)

betamean <- c(0,0,0)
betaprec <- diag(c(0.0001,0.0001,0.0001))
gammamean <- c(0,0)
gammaprec <- 0.0001*diag(2)

jags.data <- list("y"=mydata$y, "z"=mydata$z, "w1"=mydata$w1, "w2"=mydata$w2, "n"=n, "n1"=n1, "betamean"=betamean, 
                  "betaprec"=betaprec, "gammamean"=gammamean, "gammaprec"=gammaprec,
                  "tauu_alpha"=0.5, "tauu_beta"=0.5, "taux_alpha"=0.5, "taux_beta"=0.5,
                  "tauy_alpha"=0.5, "tauy_beta"=0.5)
jags.params <- c("gamma", "tauu", "taux", "tauy", "beta")

jagsmodel <- jags.model(data=jags.data, file="lin_reg.bug",n.chains=5, inits=jags.inits)
burnin <- coda.samples(jagsmodel, variable.names=c("beta"), n.iter=myburnin)
mainsample <- coda.samples(jagsmodel, variable.names=c("beta"), n.iter=mysample)
