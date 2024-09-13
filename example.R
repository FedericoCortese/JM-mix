source("Utils_JMmix.R")

# Simulate mixed data with 10% continuous missing values
TT=1000 # Number of time steps
P=30 # Number of features
Pcat=15 # Number of categorical features
Ktrue=3 # Number of true states
mu=2 # Mean value for the continuous variable 
phi=.8 # Conditional probability for the categorical outcome k in state k
rho=0 # Correlation
pers=.9 # Persistence
pNAs=.1 # Percentage of missing values
typeNA=1 # 0 for random missing, 1 for continuous missing

Y=sim_data_mixed(seed=1,
         TT=TT,
         P=P,
         Pcat=Pcat,
         Ktrue=Ktrue,
         mu=mu,
         phi=phi,
         rho=rho,
         pers=pers,
         pNAs=pNAs,
         typeNA=typeNA)

# Estimation
lambda=.16

est=jump_mixed(Y$SimData.NA,
               n_states=Ktrue,
               jump_penalty = lambda,
               verbose=F)


# Classification accuracy
adj.rand.index(Y$mchain,est$best_s)

# Imputation error
mean(gower.dist(Y$SimData.complete,est$Y))
