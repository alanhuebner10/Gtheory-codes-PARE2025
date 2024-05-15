
####Compute covariance matrix for variance components (e.g. pg 269 of Brennan) and 
   # matrix of disattenuated correlations (e.g. pg 296)
##input 1. vector of covariances by column, i.e. so entries in matrix are (1,1), (1,2), and so on
##      2. number of "strata" (i.e. multivariate responses)
Calc_CovCorrMat <- function(cov_vec, n_v){
  SIGMA_mat <- matrix(0, n_v, n_v) #set up covariance matrix
  CORR_mat <- matrix(0, n_v, n_v) #set up correlation matrix
  ####Cov matrix
  count <- 1
  for (i in 1:n_v){
    for (j in 1:n_v){
        SIGMA_mat[j, i] <- cov_vec[count]
        count <- count + 1
    }
  }
  ####Corr matrix 
  for (i in 1:n_v){
    for (j in 1:n_v){
        CORR_mat[j, i] <- SIGMA_mat[j, i]/sqrt(SIGMA_mat[j, j]*SIGMA_mat[i, i])
        count <- count + 1
    }
  }
  list_return <- list(SIGMA_mat, CORR_mat)
  names(list_return) <- c("Covariance", "Correlation")
  return(list_return)
}

####Calculate composite reliability (pg 305-306 of Brennan)
#First, calculate matrix of weights
CalcWeightMat <- function(ws){
  nv <- length(ws)
  ##Create matrix of weights
  weight_mat <-  matrix(0, nv, nv)
  for (i in 1:nv){
    for (j in 1:nv){
        weight_mat[i, j] <- ws[j]*ws[i]
    }
  }
  return(weight_mat)
}

####Calc composite
##Input vector of universe score variances, absolute error variances, and weights
CalcCompRel <- function(uscore_vec,abserr_vec, ws){
  ##Compute matrix of weights
  weight_mat <- CalcWeightMat(ws)
  #Compute cov matrix of universe scores
  uscore_mat <- Calc_CovCorrMat(uscore_vec, length(ws))[[1]]
  #Compute cov matrix of absolute error variances
  abserr_mat <- Calc_CovCorrMat(abserr_vec, length(ws))[[1]]
  ##Compute composite u score and err var
  uscore_C <- sum(weight_mat*uscore_mat)
  abserr_C <- sum(weight_mat*abserr_mat)
  #compute composite reliability
  phi_C <- uscore_C/(uscore_C + abserr_C)
  return(list(weight_mat, uscore_mat, abserr_mat,uscore_C, abserr_C, phi_C))
}


calcMGtheory_Gstudy <- function(p, i, v, y){
  dat <- data.frame(p, i, v, y)
  colnames(dat) <- c("p", "i", "v", "y")
  #obtain sample sizes
  n_p <- length(unique(dat$p))
  n_i <- length(unique(dat$i))
  n_v <- length(unique(dat$v))
  
  ##Page 291, TP stats
  #TP_p
  Xbar_pv <- dat %>% group_by(v,p) %>%
    summarise(Xbar_pv = mean(y)) 
  Xbar_pv_list <- list()
  for (i in 1:n_v){
    Xbar_pv_list <- append(Xbar_pv_list, list(Xbar_pv[Xbar_pv$v == unique(dat$v)[i], 3]))
    
  }
  TP_p_vec <- NULL
  for (i in 1:length(Xbar_pv_list)){
    for (j in 1:length(Xbar_pv_list)){
      TP_p_vec <- c(TP_p_vec, n_i*sum(Xbar_pv_list[[i]]*Xbar_pv_list[[j]]))
    }
  }
  
  #TP_i
  Xbar_iv <- dat %>% group_by(i, v) %>% 
    summarise(Xbar_iv = mean(y)) 
  Xbar_iv_list <- list()
  for (i in 1:n_v){
    Xbar_iv_list <- append(Xbar_iv_list, list(Xbar_iv[Xbar_iv$v == unique(dat$v)[i], 3]))
    
  }
  TP_i_vec <- NULL
  for (i in 1:length(Xbar_iv_list)){
    for (j in 1:length(Xbar_iv_list)){
      TP_i_vec <- c(TP_i_vec, n_p*sum(Xbar_iv_list[[i]]*Xbar_iv_list[[j]]))
    }
  }
  #TP_pi
  Xbar_piv_list <- list()
  for (i in 1:n_v){
    Xbar_piv_list <- append(Xbar_piv_list, list(dat[dat$v == unique(dat$v)[i], 4]))
    
  }
  TP_piv_vec <- NULL
  for (i in 1:n_v){
    for (j in 1:n_v){
      TP_piv_vec <- c(TP_piv_vec, sum(Xbar_piv_list[[i]]*Xbar_piv_list[[j]]))
    }
  }
  #TP_mu  
  Xbar_v <- dat %>% group_by(v) %>% 
    summarise(Xbar = mean(y)) 
  Xbar_v_list <- list()
  for (i in 1:n_v){
    Xbar_v_list <- append(Xbar_v_list, list(Xbar_v[Xbar_v$v == unique(dat$v)[i], 2]))
  }
  TP_mu_vec <- NULL
  for (i in 1:length(Xbar_v_list)){
    for (j in 1:length(Xbar_v_list)){
      TP_mu_vec <- c(TP_mu_vec, n_i*n_p*sum(Xbar_v_list[[i]]*Xbar_v_list[[j]]))
    }
  }
  ###MP statistics  
  MP_p <- (TP_p_vec - TP_mu_vec)/(n_p - 1)
  MP_i <- (TP_i_vec - TP_mu_vec)/(n_i - 1)
  MP_pi <- (TP_piv_vec - TP_p_vec - TP_i_vec + TP_mu_vec)/((n_p - 1)*(n_i - 1))
  
  ###Sigmas
  sigma_p <- (MP_p - MP_pi)/n_i 
  sigma_i <- (MP_i - MP_pi)/n_p 
  sigma_pi <- MP_pi
  
  ##Create matrices
  SIG_p <- Calc_CovCorrMat(sigma_p, n_v)
  SIG_i <- Calc_CovCorrMat(sigma_i, n_v)
  SIG_pi <- Calc_CovCorrMat(sigma_pi, n_v)
  
  list_return <- list(SIG_p, SIG_i, SIG_pi)
  names(list_return) <- c("Person", "Rater", "Residual")
  return(list_return)
}  
