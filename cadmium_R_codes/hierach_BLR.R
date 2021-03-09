##########################################################################################################################
##################taking hierarchiy into consideration###################
########################################################################################################################
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
cadmium$provinceid<-as.factor(cadmium$provinceid)
cadmium$houseid <- as.factor(cadmium$houseid)

#######centering all numeric variables
task3 <- cadmium
#task2 <- within(task2, {
# age <- as.numeric(scale(age, scale = FALSE))
#weight <- as.numeric(scale(weight, scale = FALSE))
#})

#######count of unique househould id

(NSUB1<-length(which(plyr::count(cadmium,vars="houseid") == 1))+
   length(which(plyr::count(cadmium,vars="houseid") == 2))+
   length(which(plyr::count(cadmium,vars="houseid") == 3))+
   length(which(plyr::count(cadmium,vars="houseid") == 4))+
   length(which(plyr::count(cadmium,vars="houseid") == 5))+
   length(which(plyr::count(cadmium,vars="houseid") == 6))+
   length(which(plyr::count(cadmium,vars="houseid") == 7)))

#################

#####################
task3.list <- with(task3, list(y =logy,age=age, weight=weight, sex=sex, ur=ur, SUB=houseid,Nsub=NSUB1,Nobs=nrow(task3)))

##################
task3 <- "model
	{
	for( iobs in 1 : Nobs ) 
	    {	
    	##likelihood
    	
        y[iobs] ~ dnorm(mu[iobs],tau_iobs)
              mu[iobs] <- beta0 +beta1*age[iobs] + beta2*weight[iobs]
                                  +beta3*sex[iobs] + beta4*ur[iobs] + b0[SUB[iobs]]
                                  
      # prediction new observation from current set of observations 
            newresp[iobs] ~ dnorm(mu[iobs], tau_iobs)
            

                                  
	    }
	     
	### Hyperpriors Priors
  	  tau_iobs ~ dgamma(0.001,0.001)  ###this how we specify an inverse gamma in jags
      sigma_iobs <- 1 / sqrt(tau_iobs)
      sigma_iobs2 <- 1/tau_iobs
    
	 #Random Effects
        for( isub in 1 : Nsub) {	
          b0[isub] ~ dnorm(0,taub_b02)  ###can vary prior here using t-distribution or inverser wishart

    
          ####distribution of future b0's
             b0_rep[isub] ~ dnorm(0,taub_b02)
    
	 
        }

    # Hyperpriors for Random Effects:
            sigma_b0 ~ dunif(0,100)
            sigma_b02  <- pow(sigma_b0,2)
            taub_b02 <- pow(sigma_b02,-1)
  
  #Priors Fixed effect
        beta0 ~ dt(0,0.1,100)   
        beta1 ~ dt(0,0.1,100) 
        beta2 ~ dt(0,0.1,100) 
        beta3 ~ dt(0,0.1,100) 
        beta4 ~ dt(0,0.1,100) 
                  

        
        # intra-class correlation
      r <- sigma_b02/(sigma_b02+sigma_iobs2)
      
 # PPCs checking distribution of random intercepts

          # Min and max of b0s
          ##observed bo's
                  b0.min <- min(b0[])
                  b0.max <- max(b0[])
                  
            ### replicated b0's
                   b0.rep.min <- min(b0_rep[])
                   b0.rep.max <- max(b0_rep[])
                   
           #tests
           tmin_test <-step(b0.rep.min - b0.min)
           tmax_test <- step(b0.rep.max - b0.max)
           
      # Sinharay and Stern test for b0s

           nmed <- round(Nsub/2)
           b.sort <- sort(b0)
           b_med <- b.sort[nmed]
           b.rep.sort <- sort(b0_rep)
           b.rep_med <- b.rep.sort[nmed]


          ss <- abs(b0.max - b_med) - abs(b0.min - b_med)
          ss.rep<-abs(b0.rep.max - b.rep_med)-abs(b0.rep.min - b.rep_med)

          ss.test <- step(ss.rep - ss)  	
           
      
      # Checking skewness and kurtosis of the b0's making use of posterior mean

          for (isub in 1:Nsub){
             m3.b0[isub] <- pow((b0[isub])/sigma_b0,3)
             m4.b0[isub] <- pow((b0[isub])/sigma_b0,4)
             m3.b0.rep[isub] <- pow((b0_rep[isub])/sigma_b0,3)
             m4.b0.rep[isub] <- pow((b0_rep[isub])/sigma_b0,4)
          }  

                 m3b0 <- sum(m3.b0[])/Nsub
                 m4b0 <- sum(m4.b0[])/Nsub- 3
                 m3b0.rep <- sum(m3.b0.rep[])/Nsub
                 m4b0.rep <- sum(m4.b0.rep[])/Nsub - 3
              
                 skewness.b0.test  <- step(m3b0.rep - m3b0)
                 kurtosis.b0.test <- step(m4b0.rep-m4b0)
      
      ###tests
          test[1] <- tmin_test
          test[2] <- tmax_test
          test[3] <- skewness.b0.test
          test[4] <- kurtosis.b0.test
          test[5] <- ss.test
          
      
      # PPC measures
          meas[1] 	<- b0.min
          meas.rep[1] 	<- b0.rep.min
          meas[2] 	<- b0.max
          meas.rep[2] 	<- b0.rep.max
          meas[3] 	<- m3b0
          meas.rep[3] 	<- m3b0.rep	
          meas[4] 	<- m4b0
          meas.rep[4] 	<- m4b0.rep	
          meas[5] 	<- ss
          meas.rep[5] 	<- ss.rep	
####end PPC
      

}"

set.seed(11)
params =c('beta0','beta1','beta2','beta3','beta4') 
#params =c('beta0','beta1','beta2','beta3','beta4','sigma_iobs2','sigma_b02','r') 
writeLines(task3,con='task3.bug')

HBLR <- R2jags::jags(model='task3.bug', data=task3.list, inits=NULL,param= params, 
                      n.chain=3,n.iter=70000, n.thin=5, n.burnin=15000)

stats::update(HBLR, 1e3)

HBLR_ouput_jags <- HBLR

HBLR_ouput_results <- as.data.frame(MCMCsummary(HBLR_ouput_jags, round = 3))

print(HBLR_ouput_jags,digits=3,intervals=c(0.025,0.75,0.95,0.99))

#MCMC error: Method of batch mean
mcmc_LR <- batchSE(as.mcmc(HBLR_ouput_jags), batchSize=100)
mcmc_LR<- as.vector(mcmc_LR)
formatC(mcmc_LR,format="fg")
###


### Check convergence using classical diagnostics
task3.fit.mcmc <- as.mcmc(HBLR_ouput_jags)
#par(mfrow=c(2,4))
HLR <- ggs(task3.fit.mcmc)

###trace plots
ggs_traceplot(HLR)

###autocorrelation
ggs_autocorrelation(HLR)

###Gelman and Rubin's plot
gelman.plot(task3.fit.mcmc)

##cross -correlation of parameters
ggs_crosscorrelation(HLR, absolute_scale = FALSE)


###########################################################################################################
#######################OUTPUT FOR Predicted RESPSONES#######################################################
###########################################################################################################
set.seed(11)
params =c('newresp') 
writeLines(task3,con='task3.bug')

HBLR <- R2jags::jags(model='task3.bug', data=task3.list, inits=NULL,param= params, 
                     n.chain=3,n.iter=70000, n.thin=5, n.burnin=15000)

stats::update(HBLR , 1e3)

####predicted results
HBLR_predict <- as.data.frame(MCMCsummary(HBLR , round = 3))
HBLR_predict <- HBLR_predict[-c(1),]

###########################################################################################################
#######################RANDOM INTERCEPTS#######################################################
###########################################################################################################
set.seed(11)
params =c('b0') 
writeLines(task3,con='task3.bug')

HBLR <- R2jags::jags(model='task3.bug', data=task3.list, inits=NULL,param= params, 
                     n.chain=3,n.iter=70000, n.thin=5, n.burnin=15000)

stats::update(HBLR , 1e3)

###random intercepts
HBLR_b0 <- as.data.frame(MCMCsummary(HBLR , round = 3))
HBLR_b0 <- HBLR_b0[-c(504),]
##histogram of bo

ggplot(HBLR_b0, aes(x=task3_b0$mean))+geom_histogram() + xlab("b0 means") 
ggplot(task3_b0, aes(x=task3_b0$mean))+geom_density()+ xlab("b0 means") 

##quantile plot
ggplot(HBLR_b0, aes(sample = mean))+ stat_qq() + stat_qq_line()

###########################################################################################################
#######################Discrepancy Measures#######################################################
###########################################################################################################
set.seed(11)
params =c('test') 
writeLines(task3,con='task3.bug')

HBLR <- R2jags::jags(model='task3.bug', data=task3.list, inits=NULL,param= params, 
                     n.chain=3,n.iter=70000, n.thin=5, n.burnin=15000)

stats::update(HBLR , 1e3)
###Test
HBLR_test <- as.data.frame(MCMCsummary(HBLR , round = 3))

###########################################################################################################
#######################CPO#######################################################
###########################################################################################################
set.seed(11)
params =c('icpo') 
writeLines(task3,con='task3.bug')
HBLR <- R2jags::jags(model='task2.bug', data=task2.list, inits=NULL,param= params, 
                      n.chain=3,n.iter=70000, n.thin=5, n.burnin=15000)
stats::update(HBLR , 1e3)
####output for cpo
HBLR_cpo<- as.data.frame(MCMCsummary(HBLR, round = 3))
HBLR_cpo <- HBLR_cpo[-c(1), ]

##index
#index<- as.data.frame(index <- seq(1, 1006, by=1))
#names(index)[names(index) == 'index <- seq(1, 1006, by = 1)'] <- 'index'
HBLR_cpo <- cbind.data.frame(HBLR_cpo,index)
##CPO
test_cpo$in_cpo <- 1/(test_cpo$mean)
###Influencial observations using cpo
ggplot(HBLR_cpo, aes(y=mean,x=index))+geom_point()+
  ggtitle("Index Plot of 1/CPO") +
  xlab("index") + ylab("1/CPO") 








###relplicated kurtosis
task3_m4.b0.rep <- as.data.frame(MCMCsummary(task3 , round = 3))
#observed kurtosis=0.673
task3_m4.b0.rep2 <- task3_m4.b0.rep[-c(1,505,506),]

##histogram of kurtosis
ggplot(task3_m4.b0.rep2, aes(x=mean))+geom_histogram()


























