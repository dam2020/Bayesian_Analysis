# Bayesian_Analysis
Application of Bayesian methodologies in modelling the exposure to cadmium.

# Objectives
The aim of this study, is to estimate the percentage of Chinese habitants who exhibit exposure to cadmium (Cd) greater than an internationally agreed threshold. Further, it is of interest to know what factors determine this percentage. More specifically, to estimate the impact of province, age, gender and habitat (rural/urban) of the study participant on this percentage.
In order to answer the above objectives, the analysis is divided into sub questions/objectives as follows: 
1. Across all subjects and irrespective of the covariates, what is the 95%-ile and 99%-ile of exposure?
2. What are the factors that influence the distribution of median values taking into account all the available covariates and ignoring clustering?
3. What are the factors that influence the distribution of median values taking into account all the available covariates and taking clustering into household into account?
4. What is the  75%-ile of the median exposure values of cadmium. What is the probability of exceeding this threshold and what are the factors that impact that probability taking into account clustering into household.?

# Methodology:
In order to answer the above objectives, the following models were fitted:
Bayesian mulitple linear regression
Bayesian hierarchical model with random intercepts
Bayesian logistic random intercept model (BGLMM) 

***For every model fitted***, 
* convergence is checked using classical diagnostics.
* Posterior predictive checks are used to evaluate the proposed model.
* Interpretation are provided for parameter(s) of interest.

# Further Details:
* A proper description of the dataset and explanation of the results can be found in the pdf called cadmium_Bayes. The codes used during the analysis are found in the R script called cadmium_Bayes.R .
* R2Jags was mostly used during the analysis.
