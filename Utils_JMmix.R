# List of required packages
required_packages <- c("dplyr", "cluster", "StatMatch", "pdfCluster","missMethods")

# Check if packages are installed; install if missing
missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
if(length(missing_packages)) {
  install.packages(missing_packages)
}
# Load the packages
lapply(required_packages, library, character.only = TRUE)
rm(required_packages, missing_packages)


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
                       max_iter=10, n_init=10, tol=NULL, verbose=FALSE,
                       time_vec=NULL
                     
) {
  # Fit jump model for mixed type data 
  
  # Arguments:
  # Y: data.frame with mixed data types. Categorical variables must be factors.
  # n_states: number of states
  # jump_penalty: penalty for the number of jumps
  # initial_states: initial state sequence
  # max_iter: maximum number of iterations
  # n_init: number of initializations
  # tol: tolerance for convergence
  # verbose: print progress
  # time_vec is a vector of time points, needed if times are not equally sampled
  
  # Value:
  # best_s: estimated state sequence
  # Y: imputed data
  # Y.orig: original data
  # condMM: state-conditional medians and modes
  
  timeflag=FALSE
  if(!is.null(time_vec)){
    timeflag=TRUE
      if(length(time_vec)!=nrow(Y)){
        stop("time_vec must have the same length of the number of observations")
      }
      else{
        time=sort(unique(time_vec))
        dtime=diff(time)
        dtime=dtime/as.numeric(min(dtime))
        dtime=as.numeric(dtime)
      }
  }
  
  n_states=as.integer(n_states)
  
  n_obs <- nrow(Y)
  n_features <- ncol(Y)
  Gamma <- jump_penalty * (1 - diag(n_states))
  best_loss <- NULL
  best_s <- NULL
  
  # Which vars are categorical and which are numeric
  cat_flag=any(sapply(Y, is.factor))
  
  if(cat_flag){
    cat.indx=which(sapply(Y, is.factor))
    cont.indx=which(sapply(Y, is.numeric))
    Ycont=Y[,-cat.indx]
    Ycat=Y[,cat.indx]
    
    n_levs=apply(Ycat, 2, function(x)length(unique(x[!is.na(x)])))
    # n_levs=apply(Ycat, 2, function(x)levels(x))
    
    
    n_cat=length(cat.indx)
    n_cont=n_features-n_cat
    # Initialize modes
    mo <- apply(Ycat,2,Mode)
    
    Mcat=ifelse(is.na(Ycat),T,F)
    
  }
  else{
    Ycont=Y
    n_cont=dim(Y)[2]
    n_cat=0
  }
  
  
  # Initialize mu 
  mu <- apply(Ycont, 2, median, na.rm = TRUE)
  Mcont=ifelse(is.na(Ycont),T,F)
  Ytil=Y
  
  
  
  # Impute missing values with medians of observed states
  for(i in 1:n_cont){
    Ycont[,i]=ifelse(Mcont[,i],mu[i],Ycont[,i])
  }
  
  if(cat_flag){
    for(i in 1:n_cat){
      Ycat[,i]=ifelse(Mcat[,i],mo[i],Ycat[,i])
      Ycat[,i]=factor(Ycat[,i],levels=1:n_levs[i])
    }
    Y[,-cat.indx]=Ycont
    Y[,cat.indx]=Ycat
  }
  

  # State initialization through kmeans++
  if (!is.null(initial_states)) {
    s <- initial_states
  } else {
    s=initialize_states(Y,n_states)
  }
  
  for (init in 1:n_init) {
    mu <- matrix(0, nrow=n_states, ncol=n_features-n_cat)
    
    if(cat_flag){
      mo <- matrix(0, nrow=n_states, ncol=length(cat.indx))
    }

    loss_old <- 1e10
    for (it in 1:max_iter) {
      
      # for (i in unique(s)) {
      #   
      #   mu[i,] <- apply(Ycont[s==i,], 2, median, na.rm = TRUE)
      #   if(cat_flag){
      #   mo[i,]=apply(Ycat[s==i,],2,Mode)
      #   }
      #   
      # }
      for (i in unique(s)) {
        
        if(sum(s==i)==1){
          mu[i,] <- Ycont[s==i,]
          if(cat_flag){
            mo[i,]=Ycat[s==i,]
          }
        }
        else{
          mu[i,] <- apply(Ycont[s==i,], 2, 
                          median, na.rm = TRUE)
          if(cat_flag){
            mo[i,]=apply(Ycat[s==i,],2,Mode)
          } 
        }
      }
      
      mu=data.frame(mu)
      if(cat_flag){
        mo=data.frame(mo,stringsAsFactors=TRUE)
        for(i in 1:n_cat){
          mo[,i]=factor(mo[,i],levels=1:n_levs[i])
        }
      }
      
      # Fit state sequence
      s_old <- s
      
      # Re-fill-in missings
      for(i in 1:ncol(Ycont)){
        Ycont[,i]=ifelse(Mcont[,i],mu[s,i],Ycont[,i])
      }
      if(cat_flag){
        for(i in 1:ncol(Ycat)){
          Ycat[,i]=ifelse(Mcat[,i],mo[s,i],Ycat[,i])
          Ycat[,i]=factor(Ycat[,i],levels=1:n_levs[i])
        }
        
        Y[,-cat.indx]=Ycont
        Y[,cat.indx]=Ycat
        mumo=data.frame(matrix(0,nrow=n_states,ncol=n_features))
        mumo[,cat.indx]=mo
        mumo[,cont.indx]=mu
      }
      else{
        Y=Ycont
        mumo=mu
      }
      
      
      
      
      
      # var.weights in gower.dist allows for weighted distance
      
      loss_by_state=gower.dist(Y,mumo)
      
      V <- loss_by_state
      for (t in (n_obs-1):1) {
        if(timeflag){
          V[t-1,] <- loss_by_state[t-1,] + apply(V[t,]/dtime[t] + Gamma, 2, min)
        }
        else{
        V[t-1,] <- loss_by_state[t-1,] + apply(V[t,] + Gamma, 2, min)
        }
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
  
  #Y[,cat.indx]=apply(Y[,cat.indx],2,droplevels)
  
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
                        typeNA=3){
  
  # Function to simulate mixed data with fixed parameters for the data generating process
  
  # Arguments:
  # seed: seed for the random number generator
  # TT: number of observations
  # P: number of features
  # Ktrue: number of states
  # mu: mean value for the continuous variables
  # phi: conditional probability for the categorical outcome k in state k
  # rho: correlation for the variables
  # Pcat: number of categorical variables
  # pers: self-transition probability
  # pNAs: percentage of missing values
  # typeNA is the type of missing values (0: MCAR, 1: MAR, 2: MNAR, all other values will turn into no missing imputation)
  
  # value:
  # SimData: matrix of simulated data
  
  MU=mu
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
    #SimDataCat[i, ] = SimCat[i, (Pcat * k - Pcat + 1):(Pcat * k)]
  }
  
  if(Pcat!=0){
    SimData[,1:Pcat]=apply(SimData[,1:Pcat],2,get_cat,mc=x,mu=MU,phi=phi)
    SimData=as.data.frame(SimData)
    SimData[,1:Pcat]=SimData[,1:Pcat]%>%mutate_all(as.factor)
  }
  
  if(typeNA==0){
    SimData.NA=delete_MCAR(SimData, p = pNAs)
  }
  
  else if(typeNA==1){
    vect <- 1:ncol(SimData)
    col_mis <- vect[vect %% 2 != 0]
    cols_ctrl <- vect[vect %% 2 == 0]
    SimData.NA=delete_MAR_censoring(SimData, 
                                    p = 2 * pNAs, 
                                    cols_mis = col_mis, 
                                    cols_ctrl = cols_ctrl)
  }
  
  else if(typeNA==2){
    vect <- 1:ncol(SimData)  
    SimData.NA=delete_MNAR_censoring(SimData, p = 0.2, cols_mis = 1:ncol(SimData))
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


Mode <- function(x,na.rm=T) {
  if(na.rm){
    x <- x[!is.na(x)]
  }
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
  
}
