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
typeNA=1 # 0: MCAR, 1: MAR, 2: MNAR, all other values will turn into no missing imputation

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

str(Y)

# Estimation
lambda=.16

est=jump_mixed(Y$SimData.NA,
               n_states=Ktrue,
               jump_penalty = lambda,
               verbose=F)


# Classification accuracy
adj.rand.index(Y$mchain,est$best_s)


# Data with temporal gaps
Y_gap_sim=sim_data_mixed(seed=1,
                 TT=10^4,
                 P=P,
                 Pcat=Pcat,
                 Ktrue=Ktrue,
                 mu=mu,
                 phi=phi,
                 rho=rho,
                 pers=pers,
                 pNAs=0)
Y_gap=Y_gap_sim$SimData.NA
set.seed(123)
time_gap=sort(sample(1:dim(Y_gap)[1],size=round(.9*dim(Y_gap)[1])))
Y_gap=Y_gap[-time_gap,]
time_vec=1:dim(Y_gap_sim$SimData.NA)[1]
time_vec=time_vec[-time_gap]

est_gap=jump_mixed(Y_gap,
                   n_states=Ktrue,
                   jump_penalty = lambda,
                   verbose=F,
                   time_vec=time_vec)

# Classification accuracy
mc_gap=Y_gap_sim$mchain[-time_gap]
adj.rand.index(mc_gap,est_gap$best_s)

# Results ignoring temporal gaps
est_gap_2=jump_mixed(Y_gap,
                   n_states=Ktrue,
                   jump_penalty = lambda,
                   verbose=F)

adj.rand.index(mc_gap,est_gap_2$best_s)

