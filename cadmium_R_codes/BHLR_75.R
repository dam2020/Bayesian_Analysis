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

########################################################################################
###########################75 percentile threshold#####################################
##########0.585 is the 75% of the median###################################################
###Compute the 75%-ile of the median exposure values of cadmium. What is the
#probability of exceeding this threshold and what are the factors that impact
#that probability?

cadmium$exceed75 <- ifelse(cadmium$exmd > 0.585, 1, 0)  ### exceeds_75% = 1, below_75% = 0
#cadmium$exceed75 <- as.factor(cadmium$exceed75)
#summary(cadmium$exceed75)

#######centering all numeric variables
task4 <- cadmium
task4 <- within(task4, {
age <- as.numeric(scale(age, scale = FALSE))
weight <- as.numeric(scale(weight, scale = FALSE))
})


#######count of unique househould id

(NSUB1<-length(which(plyr::count(cadmium,vars="houseid") == 1))+
    length(which(plyr::count(cadmium,vars="houseid") == 2))+
    length(which(plyr::count(cadmium,vars="houseid") == 3))+
    length(which(plyr::count(cadmium,vars="houseid") == 4))+
    length(which(plyr::count(cadmium,vars="houseid") == 5))+
    length(which(plyr::count(cadmium,vars="houseid") == 6))+
    length(which(plyr::count(cadmium,vars="houseid") == 7)))


###############
#####################
task4.list <- with(task4, list(y =exceed75,age=age, weight=weight, sex=sex, ur=ur,SUB=houseid,Nsub=NSUB1,Nobs=nrow(task4)))


########
########
task4 <- "model
	{
	for( iobs in 1 : Nobs ) 
        	{	
        	##likelihood
        	
            y[iobs] ~ dbern(p[iobs])
            logit(p[iobs]) <- beta0 +beta1*age[iobs] + beta2*weight[iobs]
                                      +beta3*sex[iobs] + beta4*ur[iobs] + b0[SUB[iobs]]
        
        	}

        	 #Random Effects
          for( isub in 1 : Nsub) {	
            b0[isub] ~ dnorm(0,taub_b02) 
            
            ####distribution of future b0's
                b0_rep[isub] ~ dnorm(0,taub_b02)
        
          }
          
      #priors for Random Effects:
              sigma_b0 ~ dunif(0,100)  ####do sensitivity analysis here
              sigma_b02  <- pow(sigma_b0,2)
              taub_b02 <- pow(sigma_b02,-1)
  
       #Priors Fixed effect
                beta0 ~ dnorm(0,1.0E-4)   
                beta1 ~ dnorm(0,1.0E-4)
                beta2 ~ dnorm(0,1.0E-4)
                beta3 ~ dnorm(0,1.0E-4)
                beta4 ~ dnorm(0,1.0E-4)
        
        
        # intra-class correlation
      #r <- sigma_b02/(sigma_b02+((3.14*3.14)/3))
      
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
######################
set.seed(11)
#params =c('beta0','beta1','beta2','beta3','beta4','sigma_b02','b0','m3b0','m3b0.rep','m4b0') 
params =c('beta0','beta1','beta2','beta3','beta4','sigma_b02') 
writeLines(task4,con='task4.bug')

BHLR_75 <- R2jags::jags(model='task4.bug', data=task4.list, inits=NULL,param= params, 
                      n.chain=3,n.iter=70000, n.thin=5, n.burnin=15000)

stats::update(BHLR_75, 1e3)

BHLR_75_ouput_jags <- BHLR_75

BHLR_75_ouput_results <- as.data.frame(MCMCsummary(BHLR_75_ouput_jags, round = 3))

print(BHLR_75_ouput_jags,digits=3,intervals=c(0.025,0.75,0.95,0.99))


#MCMC error: Method of batch mean
mcmc_BHLR_75 <- batchSE(as.mcmc(BHLR_75_ouput_jags), batchSize=100)
mcmc_BHLR_75<- as.vector(mcmc_BHLR_75)
formatC(mcmc_BHLR_75,format="fg")
###


### Check convergence using classical diagnostics
task4.fit.mcmc <- as.mcmc(BHLR_75_ouput_jags)
#par(mfrow=c(2,4))
BHLR_75 <- ggs(task4.fit.mcmc)

###trace plots
ggs_traceplot(BHLR_75)

###autocorrelation
ggs_autocorrelation(BHLR_75)

###Gelman and Rubin's plot
gelman.plot(task4.fit.mcmc)

##cross -correlation of parameters
ggs_crosscorrelation(BHLR_75, absolute_scale = FALSE)



###########################################################################################################
#######################RANDOM INTERCEPTS#######################################################
###########################################################################################################
set.seed(11)
params =c('b0') 
writeLines(task4,con='task4.bug')

BHLR_75 <- R2jags::jags(model='task4.bug', data=task4.list, inits=NULL,param= params, 
                     n.chain=3,n.iter=70000, n.thin=5, n.burnin=15000)

stats::update(HBLR , 1e3)

###random intercepts
BHLR_75_b0 <- as.data.frame(MCMCsummary(BHLR_75 , round = 3))
BHLR_75_b0 <- BHLR_75_b0[-c(504),]
##histogram of bo

ggplot(BHLR_75_b0, aes(x=BHLR_75_b0$mean))+geom_histogram() + xlab("b0 means") 
ggplot(BHLR_75_b0, aes(x=BHLR_75_b0$mean))+geom_density()+ xlab("b0 means") 

##quantile plot
ggplot(BHLR_75_b0, aes(sample = mean))+ stat_qq() + stat_qq_line()














