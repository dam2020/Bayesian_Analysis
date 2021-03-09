##################################################################

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
##################################################################

cadmium <- read.csv("C:\\Users\\NSOH TANIH\\OneDrive\\Documents\\bayesian inference\\level 2_2\\Bayesian2021\\exposure cadmium.csv",
                    sep = ';',dec = ',')

###################
#########
##Here i fit regression model with all variables,
###recoding variable  sex and ur;
cadmium$sex <- ifelse(cadmium$sex == 'Male', 1, 0)  ### Male=1, female = 0
cadmium$ur <- ifelse(cadmium$ur == 'Urban', 1, 0) #### Urban=1, rural= 0

###
cadmium$exmd <- as.numeric(cadmium$exmd)
cadmium$weight <- as.numeric(cadmium$weight)
cadmium$age <- as.numeric(cadmium$age)
cadmium$sex <- as.factor(cadmium$sex)
cadmium$ur <- as.factor(cadmium$ur)
cadmium$logy <- log(cadmium$exmd)
####

#######count of unique househould id

(NSUB1<-length(which(plyr::count(cadmium,vars="houseid") == 1))+
    length(which(plyr::count(cadmium,vars="houseid") == 2))+
    length(which(plyr::count(cadmium,vars="houseid") == 3))+
    length(which(plyr::count(cadmium,vars="houseid") == 4))+
    length(which(plyr::count(cadmium,vars="houseid") == 5))+
    length(which(plyr::count(cadmium,vars="houseid") == 6))+
    length(which(plyr::count(cadmium,vars="houseid") == 7)))

#################



#######centering all numeric variables
task2 <- cadmium
#task2 <- within(task2, {
 # age <- as.numeric(scale(age, scale = FALSE))
  #weight <- as.numeric(scale(weight, scale = FALSE))
#})

#################### making list of all variables
task2.list <- with(task2, list(y = logy,age=age, weight=weight, sex=sex, ur=ur, n = nrow(task2)))

###model
task2 <- "model {
        #Likelihood
        for (iobs in 1:n) {
          y[iobs] ~ dnorm(mu[iobs],tau) 
          
          mu[iobs] <- beta0 + beta1*age[iobs] + beta2*weight[iobs] + beta3*sex[iobs] + beta4*ur[iobs]
          
          # ordinary residual
        r[iobs] <- y[iobs] - mu[iobs]
        
     # Studentized residuals 
    t1[iobs] <- r[iobs]/sigma
    #t2[iobs] <- r[iobs]/(sigma*sqrt(1-h[iobs]))
    
    # PPOi
    log.ppo[iobs] <- -0.5*log(tau) - 0.5*tau*pow(y[iobs]-mu[iobs],2)
    ppo[iobs] <- exp(log.ppo[iobs])
    
    # basis for computing CPO
    icpo[iobs]  <- 1/ppo[iobs]
    
              }
        
                  
        #Priors
        beta0 ~ dnorm(0.0,1.0E-6)   
        beta1 ~ dnorm(0.0,1.0E-6)  
        beta2 ~ dnorm(0.0,1.0E-6)  
        beta3 ~ dnorm(0.0,1.0E-6) 
        beta4 ~ dnorm(0.0,1.0E-6) 
        tau ~ dunif(0,1000)
        sigma2 <- 1/tau
        sigma <- 1 / sqrt(tau)
       
        
                                
      }"

set.seed(11)
#params =c('beta0','beta1','beta2','beta3','beta4','tau','r','t1','ppo','icpo') 
params =c('sigma2') 
writeLines(task2,con='task2.bug')

task2 <- R2jags::jags(model='task2.bug', data=task2.list, inits=NULL,param= params, 
                      n.chain=3,n.iter=70000, n.thin=5, n.burnin=15000)
stats::update(task2, 1e3)

task2_ouput <- task2
task2_ouput_results <- as.data.frame(MCMCsummary(task2_ouput, round = 3))
####
print(task2,digits=3,intervals=c(0.025,0.75,0.95,0.99))


#MCMC error: Method of batch mean
mcmc_LR <- batchSE(as.mcmc(task2_ouput), batchSize=100)
mcmc_LR<- as.vector(mcmc_LR)
formatC(mcmc_LR,format="fg")


###########################################################################################################
#######################Convergence for linear model#######################################################
###########################################################################################################

### Check convergence using classical diagnostics
task2.fit.mcmc <- as.mcmc(task2)
#par(mfrow=c(2,4))
LR <- ggs(task2.fit.mcmc)

###trace plots
ggs_traceplot(LR)
#traceplot(task2.fit.mcmc)
#xyplot(task2.fit.mcmc, layout=c(3,4), aspect="fill")
#plot(task2.fit.mcmc)
###Gelman and Rubin's plot
gelman.plot(task2.fit.mcmc)
###autocorrelation
ggs_autocorrelation(LR)
#autocorr.plot(task2.fit.mcmc,auto.layout = FALSE)
##cross -correlation of parameters
ggs_crosscorrelation(LR,absolute_scale = FALSE)


###########################################################################################################
#######################OUTPUT FOR STANDARD RESIDUAL#######################################################
###########################################################################################################
set.seed(11)
params =c('t1') 
writeLines(task2,con='task2.bug')
task2 <- R2jags::jags(model='task2.bug', data=task2.list, inits=NULL,param= params, 
                      n.chain=3,n.iter=70000, n.thin=5, n.burnin=15000)
stats::update(task2, 1e3)
##index
index<- as.data.frame(index <- seq(1, 1006, by=1))
names(index)[names(index) == 'index <- seq(1, 1006, by = 1)'] <- 'index'
####output for standardised residuals(r)
test_stand <- as.data.frame(MCMCsummary(task2, round = 3))
test_stand <- test_stand[-c(1), ]
test_stand <- cbind.data.frame(test_stand,index)
##standardized residuals PLOT
ggplot(test_stand, aes(y=mean,x=index))+geom_point()+
  geom_hline(yintercept=2,color ='blue')+geom_hline(yintercept=-2,color ='blue')+
  ggtitle("Detecting Outlying Standardized Residuals") +
  xlab("index") + ylab("Standardized Residuals")

###########################################################################################################
#######################creating a new dataset without outliers to test for inflence#######################################################
###########################################################################################################
test_stand <- cbind.data.frame(test_stand,cadmium)
##deleting observations >2 and <-2
infle <- subset(test_stand, mean < 2 & mean > -2)

task22.list <- with(infle, list(y = logy,age=age, weight=weight, sex=sex, ur=ur, n = nrow(infle)))
###
set.seed(11)
params =c('beta0','beta1','beta2','beta3','beta4','sigma2') 
writeLines(task2,con='task2.bug')
task22 <- R2jags::jags(model='task2.bug', data=task22.list, inits=NULL,param= params, 
                      n.chain=3,n.iter=70000, n.thin=5, n.burnin=15000)
stats::update(task22, 1e3)

task22_ouput <- task22
#MCMC error: Method of batch mean
mcmc_LR <- batchSE(as.mcmc(task22), batchSize=100)
mcmc_LR<- as.vector(mcmc_LR)
formatC(mcmc_LR,format="fg")




###########################################################################################################
#######################OUTPUT FOR CPO#######################################################
###########################################################################################################
set.seed(11)
params =c('icpo') 
writeLines(task2,con='task2.bug')
task2 <- R2jags::jags(model='task2.bug', data=task2.list, inits=NULL,param= params, 
                      n.chain=3,n.iter=70000, n.thin=5, n.burnin=15000)
stats::update(task2, 1e3)
####output for cpo
test_cpo<- as.data.frame(MCMCsummary(task2, round = 3))
test_cpo <- test_cpo[-c(1), ]
test_cpo<- cbind.data.frame(test_cpo,index)
##CPO
test_cpo$in_cpo <- 1/(test_cpo$mean)
###Influencial observations using cpo
ggplot(test_cpo, aes(y=mean,x=index))+geom_point()+
  ggtitle("Index Plot of 1/CPO") +
  xlab("index") + ylab("1/CPO") 
##+scale_y_continuous(limits = c(0, 250000))

###########################################################################################################
#######################OUTPUT FOR PPO#######################################################
###########################################################################################################
set.seed(11)
params =c('ppo') 
writeLines(task2,con='task2.bug')
task2 <- R2jags::jags(model='task2.bug', data=task2.list, inits=NULL,param= params, 
                      n.chain=3,n.iter=70000, n.thin=5, n.burnin=15000)
stats::update(task2, 1e3)
####output for ppo
test_ppo<- as.data.frame(MCMCsummary(task2, round = 3))
test_ppo <- test_ppo[-c(1), ]
test_ppo<- cbind.data.frame(test_ppo,index)
summary(test_ppo$mean)
test_ppo$std <- ((test_ppo$mean)/(0.3300)) 

####making relative difference of cpo vs ppo plot
test_cpo$mean_cpo <- test_cpo$mean
test_ppo$mean_ppp <- test_ppo$mean
ppo_cpo <- cbind.data.frame(test_ppo,test_cpo)
ppo_cpo$new <- (100*(test_cpo$in_cpo - test_ppo$mean)/(test_ppo$mean))

ppo_cpo2 <- subset.data.frame(ppo_cpo, select = c("mean_ppp","new"))


###Relative difference of CPO vs PPO 
ggplot(ppo_cpo2, aes(y=ppo_cpo$new,x=ppo_cpo$mean_ppp))+geom_point()+
  ggtitle("Relative difference of CPO vs PPO") +
  xlab("PPO") + ylab("(CPO-PPO)/PP0 %") + scale_y_continuous(limits = c(-20, 0))

###Index Plot of ppo 
ggplot(test_ppo, aes(y=std,x=index))+geom_point()+
  ggtitle("Index Plot Of PPO") +
  xlab("Index") + ylab("PP0")



+ scale_y_continuous(limits = c(-20, 0))







