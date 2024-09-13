source("Utils_JMmix.R")

# Simulate mixed data with 10% continuous missing values
TT=1000
P=30 
Pcat=15
Ktrue=3
mu=2
phi=.8
rho=0
pers=.9
pNAs=.1 
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

# Estimate
lambda=.16

est=jump_mixed(Y$SimData.NA,
               n_states=Ktrue,
               jump_penalty = lambda,
               verbose=F)


# Classification accuracy
adj.rand.index(Y$mchain,est$best_s)

# Imputation error
mean(gower.dist(Y$SimData.complete,est$Y))
