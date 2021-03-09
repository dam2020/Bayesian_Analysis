library(tidyverse)
library(dplyr)
library(plyr)
library(R2jags)
library(ggplot2)
library(nlme)
library(standardize)
library(MCMCvis)
library(bayesplot)
library(coda)
library(batchmeans)
library(lattice)
library(ggmcmc)
##################################################################

#cadmium <- read.csv("C:\\Users\\NSOH TANIH\\OneDrive\\Documents\\bayesian inference\\level 2_2\\Bayesian2021\\exposure cadmium.csv",
#sep = ';')
cadmium <- read.csv("C:\\Users\\NSOH TANIH\\OneDrive\\Documents\\bayesian inference\\level 2_2\\Bayesian2021\\exposure cadmium.csv",
                    sep = ';',dec = ',')

str(cadmium)
###
cadmium$exmd <- as.numeric(cadmium$exmd)
cadmium$weight <- as.numeric(cadmium$weight)
mean(cadmium$exmd)
####
summary(cadmium$exmd)
summary(cadmium$person)
summary(cadmium$age)
summary(cadmium$weight)
summary(cadmium$sex)
summary(cadmium$ur)

######Histogram
hist(cadmium$exmd,nclass = 5, probability=T)
ggplot(cadmium, aes(x=logy))+geom_histogram()
##ggplot(cadmium, aes(x=exmd))+geom_histogram(aes(y=..density..))+geom_density()

####Density
ggplot(cadmium, aes(x=logy))+geom_density()
ggplot(cadmium, aes(x=age))+geom_density()
ggplot(cadmium, aes(x=weight))+geom_density()

####histograms
ggplot(cadmium, aes(y=logy,x=sex))+geom_boxplot()
ggplot(cadmium, aes(y=exmd,x=ur))+geom_boxplot()

################################################################################
###############Without any covariate###########################################
################################################################################
#######centering all numeric variables
task1 <- cadmium

#task1 <- within(task1, {
#age <- as.numeric(scale(age, scale = FALSE))
# weight <- as.numeric(scale(weight, scale = FALSE))
#})

####################
task1.list <- with(task1, list(y = logy, n = nrow(task1)))
#task1.list

task1 <- "model {
        #Likelihood
        for (i in 1:n) {
          y[i] ~ dnorm(mu[i],tau)
          mu[i] <- beta0
          y.err[i] <- y[i] - mu[i]
        }
        
        #Priors
        beta0 ~ dnorm(0.0,1.0E-6)
        tau <- 1 / (sigma * sigma)
        sigma2 <- 1/tau
        sigma~dgamma(0.001,0.001)
      }"

set.seed(11)
params =c('beta0','sigma2') 
writeLines(task1,con='task1.bug')

task1 <- R2jags::jags(model='task1.bug', data=task1.list, inits=NULL,param= params, 
                      n.chain=3,n.iter=50000, n.thin=5, n.burnin=15000)

stats::update(task1, 1e3)
print(task1,digits=3,intervals=c(0.025,0.75,0.95,0.99))

#MCMC error: Method of batch mean
mcmc_task1 <- batchSE(as.mcmc(task1), batchSize=100)
mcmc_task1 <- as.vector(mcmc_task1 )
formatC(mcmc_task1 ,format="fg")

###
task1_results <- as.data.frame(MCMCsummary(task1 , round = 3))

### Check convergence using classical diagnostics
task1.fit.mcmc <- as.mcmc(task1)
#par(mfrow=c(2,4))
mod1 <- ggs(task1.fit.mcmc)

###trace plots
ggs_traceplot(mod1)

###autocorrelation
ggs_autocorrelation(mod1)

###Gelman and Rubin's plot
gelman.plot(task1.fit.mcmc)

##cross -correlation of parameters
ggs_crosscorrelation(mod1)

##ggs_ppmean(mod1,outcome = 'mu')
###Histograms
ggplot(cadmium, aes(x=logy))+geom_histogram() + geom_vline(xintercept=-0.613,color ='blue') +
  xlab("log Cadmium")