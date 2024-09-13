library(dplyr)
library(cluster)
library(StatMatch)
library(pdfCluster)

initialize_states <- function(Y, K) {
  
  # Initialize states for jump_mixed
  
  # Arguments:
  # Y: data.frame with mixed data types. Categorical variables must be factors.
  # K: number of states
  
  
  n <- nrow(Y)
  
  centr_indx=sample(1:n, 1)
  centroids <- Y[centr_indx, , drop = FALSE]  # Seleziona il primo centroide a caso
  
  closest_dist <- as.matrix(daisy(Y, metric = "gower"))
  closest_dist <- closest_dist[centr_indx,]
  
  for (i in 2:K) {
    prob <- closest_dist / sum(closest_dist)
    next_centr_indx <- sample(1:n, 1, prob = prob)
    next_centroid <- Y[next_centr_indx, , drop = FALSE]
    centroids <- rbind(centroids, next_centroid)
  }
  
  dist_matrix <- gower.dist(Y, centroids)
  init_stats <- apply(dist_matrix, 1, which.min)
  
  return(init_stats)
}

jump_mixed <- function(Y, n_states, jump_penalty=1e-5, 
                       initial_states=NULL,
                       max_iter=10, n_init=10, tol=NULL, verbose=FALSE
                     
) {
  # Fit jump model for muxed type data 
  
  # Arguments:
  # Y: data.frame with mixed data types. Categorical variables must be factors.
  # n_states: number of states
  # jump_penalty: penalty for the number of jumps
  # initial_states: initial state sequence
  # max_iter: maximum number of iterations
  # n_init: number of initializations
  # tol: tolerance for convergence
  # verbose: print progress
  
  # Value:
  # best_s: estimated state sequence
  # Y: imputed data
  # Y.orig: original data
  # condMM: state-conditional means and modes
  
  n_states=as.integer(n_states)
  
  n_obs <- nrow(Y)
  n_features <- ncol(Y)
  Gamma <- jump_penalty * (1 - diag(n_states))
  best_loss <- NULL
  best_s <- NULL
  
  # Which vars are categorical and which are numeric
  
  cat.indx=which(sapply(Y, is.factor))
  cont.indx=which(sapply(Y, is.numeric))
  Ycont=Y[,-cat.indx]
  Ycat=Y[,cat.indx]
  
  n_levs=apply(Ycat, 2, function(x)length(unique(x)))
  # n_levs=apply(Ycat, 2, function(x)levels(x))
  
  
  n_cat=length(cat.indx)
  n_cont=n_features-n_cat
  
  # Initialize mu 
  mu <- colMeans(Ycont,na.rm = T)
  
  # Initialize modes
  mo <- apply(Ycat,2,Mode)
  
  # Track missings with 0 1 matrix
  Mcont=ifelse(is.na(Ycont),T,F)
  Mcat=ifelse(is.na(Ycat),T,F)
  
  Ytil=Y
  # Impute missing values with mean of observed states
  for(i in 1:n_cont){
    Ycont[,i]=ifelse(Mcont[,i],mu[i],Ycont[,i])
  }
  for(i in 1:n_cat){
    Ycat[,i]=ifelse(Mcat[,i],mo[i],Ycat[,i])
    Ycat[,i]=factor(Ycat[,i],levels=1:n_levs[i])
  }
  
  
  Y[,-cat.indx]=Ycont
  Y[,cat.indx]=Ycat
  
  # State initialization through kmeans++
  if (!is.null(initial_states)) {
    s <- initial_states
  } else {
    s=initialize_states(Y,n_states)
  }
  
  for (init in 1:n_init) {
    mu <- matrix(0, nrow=n_states, ncol=n_features-length(cat.indx))
    mo <- matrix(0, nrow=n_states, ncol=length(cat.indx))
    loss_old <- 1e10
    for (it in 1:max_iter) {
      
      for (i in unique(s)) {
        mu[i,] <- colMeans(Ycont[s==i,])
        mo[i,]=apply(Ycat[s==i,],2,Mode)
        
      }
      
      mu=data.frame(mu)
      mo=data.frame(mo,stringsAsFactors=TRUE)
      for(i in 1:n_cat){
        mo[,i]=factor(mo[,i],levels=1:n_levs[i])
      }
      
      # Fit state sequence
      s_old <- s
      
      # Re-fill-in missings
      for(i in 1:ncol(Ycont)){
        Ycont[,i]=ifelse(Mcont[,i],mu[s,i],Ycont[,i])
      }
      for(i in 1:ncol(Ycat)){
        Ycat[,i]=ifelse(Mcat[,i],mo[s,i],Ycat[,i])
        Ycat[,i]=factor(Ycat[,i],levels=1:n_levs[i])
      }
      
      Y[,-cat.indx]=Ycont
      Y[,cat.indx]=Ycat
      
      mumo=data.frame(matrix(0,nrow=n_states,ncol=n_features))
      mumo[,cat.indx]=mo
      mumo[,cont.indx]=mu
      
      
      
      # var.weights in gower.dist allows for weighted distance
      
      loss_by_state=gower.dist(Y,mumo)
      
      V <- loss_by_state
      for (t in (n_obs-1):1) {
        V[t-1,] <- loss_by_state[t-1,] + apply(V[t,] + Gamma, 2, min)
      }
      
      s[1] <- which.min(V[1,])
      for (t in 2:n_obs) {
        s[t] <- which.min(V[t,] + Gamma[s[t-1],])
      }
      
      if (length(unique(s)) == 1) {
        break
      }
      loss <- min(V[1,])
      if (verbose) {
        cat(sprintf('Iteration %d: %.6e\n', it, loss))
      }
      if (!is.null(tol)) {
        epsilon <- loss_old - loss
        if (epsilon < tol) {
          break
        }
      } else if (all(s == s_old)) {
        break
      }
      loss_old <- loss
    }
    if (is.null(best_s) || (loss_old < best_loss)) {
      best_loss <- loss_old
      best_s <- s
    }
    #s <- init_states(Y, n_states)+1
    s=initialize_states(Y,n_states)
  }
  
  return(list(best_s=best_s,
              Y=Y,
              Y.orig=Ytil,
              condMM=mumo))
}

sim_data_mixed=function(seed=123,
                        TT,
                        P,
                        Ktrue=3,
                        mu=1,
                        phi=.8,
                        rho=0,
                        Pcat=NULL,
                        pers=.95,
                        pNAs=0,
                        typeNA=2){
  
  # Function to simulate mixed data with fixed parameters for the data generating process
  
  # Arguments:
  # seed: seed for the random number generator
  # TT: number of observations
  # P: number of features
  # Ktrue: number of states (only 3 states are allowed TO BE UPDATED)
  # mu: mean value for the continuous variables
  # phi: conditional probability for the categorical outcome k in state k
  # rho: correlation for the variables
  # Pcat: number of categorical variables
  # pers: self-transition probability
  # pNAs: percentage of missing values
  # typeNA is the type of missing values (0 for random, 1 for continuous, all other values will turn into no missing imputation)
  
  # value:
  # SimData.NA: matrix of simulated data with missing values
  # SimData: matrix of simulated data wihtout missing values
  # mchain: latent Markov chain
  
  mu=c(-mu,0,mu)
  
  if(is.null(Pcat)){
    Pcat=floor(P/2)
  }
  
  # Markov chain simulation
  x <- numeric(TT)
  Q <- matrix(rep((1-pers)/(Ktrue-1),Ktrue*Ktrue), 
              ncol = Ktrue,
              byrow = TRUE)
  diag(Q)=rep(pers,Ktrue)
  init <- rep(1/Ktrue,Ktrue)
  set.seed(seed)
  x[1] <- sample(1:Ktrue, 1, prob = init)
  for(i in 2:TT){
    x[i] <- sample(1:Ktrue, 1, prob = Q[x[i - 1], ])
  }
  
  # Continuous variables simulation
  Sigma <- matrix(rho,ncol=P,nrow=P)
  diag(Sigma)=1
  
  Sim = matrix(0, TT, P * Ktrue)
  SimData = matrix(0, TT, P)
  
  set.seed(seed)
  for(k in 1:Ktrue){
    u = MASS::mvrnorm(TT,rep(mu[k],P),Sigma)
    Sim[, (P * k - P + 1):(k * P)] = u
  }
  
  for (i in 1:TT) {
    k = x[i]
    SimData[i, ] = Sim[i, (P * k - P + 1):(P * k)]
  }
  
  if(Pcat!=0){
    SimData[,1:Pcat]=apply(SimData[,1:Pcat],2,get_cat,mc=x,mu=mu,phi=phi)
    SimData=as.data.frame(SimData)
    SimData[,1:Pcat]=SimData[,1:Pcat]%>%mutate_all(as.factor)
  }
  
  if(typeNA==0|typeNA==1){
    SimData.NA=apply(SimData,2,punct,pNAs=pNAs,type=typeNA)
    SimData.NA=as.data.frame(SimData.NA)
    SimData.NA[,1:Pcat]=SimData.NA[,1:Pcat]%>%mutate_all(as.factor)
    SimData.NA[,-(1:Pcat)]=SimData.NA[,-(1:Pcat)]%>%mutate_all(as.numeric)
  }
  else{
    SimData.NA=SimData
  }
  
  return(list(
    SimData.NA=SimData.NA,
    SimData.complete=SimData,
    mchain=x,
    TT=TT,
    P=P,
    Ktrue=Ktrue,
    pers=pers, 
    seed=seed))
  
}

get_cat=function(y,mc,mu,phi){
  # Function to simulate categorical data
  
  # Arguments:
  # x: continuous variable 
  # mc: Markov chain states
  # mu: numeric mean value
  # phi: conditional probability for the categorical outcome k in state k
  
  mu=c(-mu,0,mu)
  phi1=(1-phi)/2
  
  TT=length(y)
  for(i in 1:TT){
    k=mc[i]
    switch(k,
           "1"={
             threshold=c(qnorm(phi1,mu[1]),qnorm(phi+phi1,mu[1]))
             if(y[i]>threshold[1]&y[i]<threshold[2]){
               y[i]=1
             }
             else if(y[i]<threshold[1]){
               y[i]=2
             }
             else{
               y[i]=3
             }
           },
           "2"={
             threshold=c(qnorm(phi1,mu[2]),qnorm(phi+phi1,mu[2]))
             if(y[i]>threshold[1]&y[i]<threshold[2]){
               y[i]=2
             }
             else if(y[i]<threshold[1]){
               y[i]=3
             }
             else{
               y[i]=1
             }
           },
           "3"={
             threshold=c(qnorm(phi1,mu[3]),qnorm(phi+phi1,mu[3]))
             if(y[i]>threshold[1]&y[i]<threshold[2]){
               y[i]=3
             }
             else if(y[i]<threshold[1]){
               y[i]=1
             }
             else{
               y[i]=2
             }
           }
    )
  }
  return(y)
  
}

punct=function(x,pNAs,typeNA){
  
  # x is a vector (column of the dataset)
  # pNAs is the percentage of missing values
  # typeNA is the type of missing values (0 for random, 1 for continuous, all other values will turn into no missing imputation)
  
  TT=length(x)
  pTT=round(TT*pNAs)
  if(typeNA==0){
    NAindx=sample(1:TT,pTT,replace = F)
    x[NAindx]=NA
  }
  else if(typeNA==1){
    NAindx=sample(1:(TT-pTT),1,replace = F)
    NAindx=seq(NAindx,NAindx+pTT)
    x[NAindx]=NA
  }
  
  return(x)
  
}

Mode <- function(x,na.rm=T) {
  if(na.rm){
    x <- x[!is.na(x)]
  }
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
  
}
