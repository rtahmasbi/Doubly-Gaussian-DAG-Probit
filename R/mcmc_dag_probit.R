# NOTE: I set
#      gamma = 0 # Rasool
#      Sigma1_true, A = A1_true
#            out1_t = cbind(0, sapply(FUN = causal_y, X = 2:q, Sigma_D = Sigma1_true, A = A1_true, X_mat = X1, g = gamma))
      

library(pcalg)
library(gRbase)
library(truncnorm)
library(fastmatrix)
library(DensParcorr)

source("move_dag_probit.R")
source("marg_like_dag.R")
source("posterior_sigma_dag_probit.R")
source("causal_effect_probit.R")
source("update_DAG.R")


# X is p*p not q*q
make_init = function(X, thr=.3) {
	dec = ldl(solve(cov(X)))$lower
	A = (abs(dec)>=thr)+0
	diag(A) = 0
	return (A)
}


mcmc_dag_probit = function(y1, y2, X1, X2, TT, burn, a = NULL, g1 = NULL, g2 = NULL, w = NULL, causal = TRUE, gamma_0 = 0, debug = FALSE, A1_true=NULL, Sigma1_true=NULL ){
  
  ###########
  ## INPUT ##
  ###########
  
  # y :    (n,1) vector of observations from the binary response variable
  # X :    (n,p) matrix of observations from the p covariates, p=q-1
  # TT :   number of MCMC iterations
  # burn : burn-in period
  
  # a,g : hyperparamters of the DAG-Wishart prior
  # w :   prior probability of edge inclusion
  
  # causal : logical value (if TRUE, BMA causal effects are also computed; if FALSE, they are not)
  
  ############
  ## OUTPUT ##
  ############
  
  # A_chain :     matrix with TT columns, each representing the (vectorized) (q,q) adjacency matrix of one visited DAG
  # Sigma_chain : matrix with TT columns, each representing the (vectorized) (q,q) covariance matrix of one visited DAG
  # Causal_hat :  (n,q) matrix with BMA estimates of subject-specific causal effects (one row for each subject, one column for each intervened node)
  
  
  n1 = nrow(X1)
  n2 = nrow(X2)
  p = ncol(X1)
  q = p + 1
  print(c(n1,n2,q))
  
  # Set DAG-Wishart hyperparameters (if not specified)
  
  if(is.null(a)){a = q}
  
  if(is.null(g1)){g1 = 1/n1}
  if(is.null(g2)){g2 = 1/n2}
  
  if(is.null(w)){w = 0.5}
  
  # Store space for parameters
  
  A1_chain = matrix(0, nrow = q*q, ncol = TT)
  A2_chain = matrix(0, nrow = q*q, ncol = TT)
  # matrix collecting the adjacecncy matrices of the DAGs
  # (each column is the by-column-vectorized adjacency matrix of the accepted DAG)
  
  Sigma1_chain = matrix(0, nrow = q*q, ncol = TT)
  Sigma2_chain = matrix(0, nrow = q*q, ncol = TT)
  # matrix collecting the posterior draws from Sigma
  # (each column is the by-column-vectorized covariance matrix Sigma)

  Partial1_chain = matrix(0, nrow = q*q, ncol = TT)
  Partial2_chain = matrix(0, nrow = q*q, ncol = TT)
  
  L1_chain = matrix(0, nrow = q*q, ncol = TT)
  L2_chain = matrix(0, nrow = q*q, ncol = TT)
  D_chain = matrix(0, nrow = q*q, ncol = TT)
  
  z1_chain = matrix(0, nrow = n1, ncol = TT)
  z2_chain = matrix(0, nrow = n2, ncol = TT)
  # matrix collecting the posterior draws from the latent z
  
  Gamma_chain = rep(0, TT)
  # matrix collecting posterior samples for the cutoff gamma
    
  inf_lim1 = rep(-Inf, n1)
  sup_lim1 = rep(Inf, n1)
  inf_lim2 = rep(-Inf, n2)
  sup_lim2 = rep(Inf, n2)
  
  inf_lim1[y1 == 1] = gamma_0
  sup_lim1[y1 == 0] = gamma_0
  inf_lim2[y2 == 1] = gamma_0
  sup_lim2[y2 == 0] = gamma_0
  
  
  if(causal == TRUE){
    
    Causal1_hat = matrix(0, n1, q)
    Causal2_hat = matrix(0, n2, q)
    # matrix collecting BMA causal effect estimates
    
  }else{
    
    Causal1_hat = NULL
    Causal2_hat = NULL
    
  }
  
  ## Inits values
  
  A_0 = matrix(0, q, q)
  colnames(A_0) = rownames(A_0) = 1:q
  A1 = A2 = A_0
  A1[2:nrow(A1), 2:ncol(A1)] = make_init(X1)
  A2[2:nrow(A2), 2:ncol(A2)] = make_init(X2)
  A1_chain[,1] = A1
  A2_chain[,1] = A2
  
  z1_0 = rtruncnorm(n1, a = inf_lim1, b = sup_lim1, mean = 0, sd = 1);
  z2_0 = rtruncnorm(n2, a = inf_lim2, b = sup_lim2, mean = 0, sd = 1);
  z1 = z1_0
  z2 = z2_0
  z1_chain[,1] = z1
  z2_chain[,1] = z2
  
  
  Chol_Sigma_0 = posterior_sigma(Y1=cbind(z1,X1), Y2=cbind(z2,X2), A1, A2, g1, g2, a)
  Sigma1_chain[,1] = Chol_Sigma_0$Sigma1_post
  Sigma2_chain[,1] = Chol_Sigma_0$Sigma2_post
  Partial1_chain[,1] = prec2part(Chol_Sigma_0$Omega1_post)
  Partial2_chain[,1] = prec2part(Chol_Sigma_0$Omega2_post)
  
  #L_0    = -Chol_Sigma_0$L_post[-1,1]
  #D_1_0  = Chol_Sigma_0$D_post[1,1]
  
  
  gamma = gamma_0
  
  
  ## MCMC iterations
  
  cat("MCMC sampling")
  pb = utils::txtProgressBar(min = 2, max = TT, style = 3)
  
  for(t in 1:TT) {
    
    ## update of DAG D
    A1 = update_DAG(A1, z1, X1, a, g1, w)
    A2 = update_DAG(A2, z2, X2, a, g2, w)
    #A1 = A1_orig
    #A2 = A2_orig

    ## Sample Sigma    
    Chol_Sigma = posterior_sigma(Y1 = cbind(z1,X1), Y2 = cbind(z2,X2), A1, A2, g1, g2, a)
    Sigma1_post = Chol_Sigma$Sigma1_post
    Sigma2_post = Chol_Sigma$Sigma2_post
    L1_post     = Chol_Sigma$L1_post
    L2_post     = Chol_Sigma$L2_post
    D_post      = Chol_Sigma$D_post
    Omega1_post = Chol_Sigma$Omega1_post
    Omega2_post = Chol_Sigma$Omega2_post
    
    ## Sample and update the latent
    Beta1 = -L1_post[-1,1]
    Beta2 = -L2_post[-1,1]
    z1 = rtruncnorm(n1, a = inf_lim1, b = sup_lim1, mean = X1%*%Beta1, sd = sqrt(D_post[1,1]))
    z2 = rtruncnorm(n2, a = inf_lim2, b = sup_lim2, mean = X2%*%Beta2, sd = sqrt(D_post[1,1]))
    #z1 = y1_true
    #z2 = y2_true
    
    # Update the cut-off
    
    g_prop = rnorm(1, gamma, 0.5)
    
    inf_lim_g1 = inf_lim1
    inf_lim_g1[y1 == 1] = g_prop
    sup_lim_g1 = sup_lim1
    sup_lim_g1[y1 == 0] = g_prop
    #
    inf_lim_g2 = inf_lim2
    inf_lim_g2[y2 == 1] = g_prop
    sup_lim_g2 = sup_lim2
    sup_lim_g2[y2 == 0] = g_prop
    
    # (sup_lim_g1, inf_lim_g1)=(theta,-infty) or (infty,theta)
    r_g = sum(log(pnorm(sup_lim_g1, X1%*%Beta1, sd = sqrt(D_post[1,1])) - pnorm(inf_lim_g1, X1%*%Beta1, sd = sqrt(D_post[1,1])))) -
          sum(log(pnorm(sup_lim1,   X1%*%Beta1, sd = sqrt(D_post[1,1])) - pnorm(inf_lim1,   X1%*%Beta1, sd = sqrt(D_post[1,1])))) +
	      sum(log(pnorm(sup_lim_g2, X2%*%Beta2, sd = sqrt(D_post[1,1])) - pnorm(inf_lim_g2, X2%*%Beta2, sd = sqrt(D_post[1,1])))) -
          sum(log(pnorm(sup_lim2,   X2%*%Beta2, sd = sqrt(D_post[1,1])) - pnorm(inf_lim2,   X2%*%Beta2, sd = sqrt(D_post[1,1])))) +
          dnorm(g_prop, gamma, 0.5, log = TRUE) - dnorm(gamma, g_prop, 0.5, log = TRUE) # rasool change it from g to g_prop
    
    #if (is.nan(r_g)) { # if it is Nan, dont accept it
    #	r_g = -Inf
    #}
    # acceptance ratio
    
    ratio_g = min(0, r_g)
    if (debug) {
	    print('-----------------------------------------------------------------')
	    print(c(gamma, g_prop))
	    print(Beta1)
	    print(Beta2)
	    print(r_g)
	    print(ratio_g)
	    print('-----')
	    print(sum(log(pnorm(sup_lim_g1, X1%*%Beta1, sd = sqrt(D_post[1,1])) - pnorm(inf_lim_g1, X1%*%Beta1, sd = sqrt(D_post[1,1])))))
	    print(sum(log(pnorm(sup_lim1,   X1%*%Beta1, sd = sqrt(D_post[1,1])) - pnorm(inf_lim1,   X1%*%Beta1, sd = sqrt(D_post[1,1])))))
	    print(sum(log(pnorm(sup_lim_g2, X2%*%Beta2, sd = sqrt(D_post[1,1])) - pnorm(inf_lim_g2, X2%*%Beta2, sd = sqrt(D_post[1,1])))))
	    print(sum(log(pnorm(sup_lim2,   X2%*%Beta2, sd = sqrt(D_post[1,1])) - pnorm(inf_lim2,   X2%*%Beta2, sd = sqrt(D_post[1,1])))))
	    print('-----')
	    print(c(min(X1%*%Beta1),max(X1%*%Beta1)))
	    print(c(min(X2%*%Beta2),max(X2%*%Beta2)))
	    #print(c(min(sup_lim_g1-inf_lim_g1), max(sup_lim_g1-inf_lim_g1)))
	    #print(c(min(sup_lim_g2-inf_lim_g2), max(sup_lim_g2-inf_lim_g2)))
	    print('-----')
	    bad_idx1 = which(pnorm(sup_lim_g1, X1%*%Beta1, sd = sqrt(D_post[1,1])) - pnorm(inf_lim_g1, X1%*%Beta1, sd = sqrt(D_post[1,1]))==0)
	    print(bad_idx1)
	    print(y1[bad_idx1])
	    print((X1%*%Beta1)[bad_idx1])
	    print(inf_lim_g1[bad_idx1])
	    print(sup_lim_g1[bad_idx1])
    }
    
    # accept move
    if(log(runif(1)+1e-20) < ratio_g){
      gamma   = g_prop
      inf_lim1 = inf_lim_g1
      sup_lim1 = sup_lim_g1
      inf_lim2 = inf_lim_g2
      sup_lim2 = sup_lim_g2
    }
    
    # store chain values
    
    A1_chain[,t]     = A1
    A2_chain[,t]     = A2
    Sigma1_chain[,t] = Sigma1_post
    Sigma2_chain[,t] = Sigma2_post
    Partial1_chain[,t] = prec2part(Omega1_post)
    Partial2_chain[,t] = prec2part(Omega2_post)
    L1_chain[,t]     = L1_post
    L2_chain[,t]     = L2_post
    D_chain[,t]      = D_post
    Gamma_chain[t]   = gamma
    
    ## compute causal effects (if causal = TRUE)
    
    if(causal == TRUE){
      if(t > burn){
        out1_t = cbind(0, sapply(FUN = causal_y, X = 2:q, Sigma_D = Sigma1_post, A = A1, X_mat = X1, g = gamma))
        out2_t = cbind(0, sapply(FUN = causal_y, X = 2:q, Sigma_D = Sigma2_post, A = A2, X_mat = X2, g = gamma))
        Causal1_hat = Causal1_hat + out1_t/(TT - burn)
        Causal2_hat = Causal2_hat + out2_t/(TT - burn)
      }
    }

    utils::setTxtProgressBar(pb, t)

    close(pb)
    
    
  }
 
  
  return(list(A1_chain     = A1_chain,
			  A2_chain     = A2_chain, 
              Sigma1_chain = Sigma1_chain,
              Sigma2_chain = Sigma2_chain,
              Partial1_chain = Partial1_chain,
              Partial2_chain = Partial2_chain,
              L1_chain     = L1_chain,
              L2_chain     = L2_chain,
              D_chain      = D_chain,
              Gamma_chain  = Gamma_chain,
              Causal1_hat  = Causal1_hat,
              Causal2_hat  = Causal2_hat
  ))
}
