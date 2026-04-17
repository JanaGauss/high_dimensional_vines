source("functions/helpers.R")

verify_assumption_curvature <- function(d, K, par, dvine, name, error_const, 
                                        alpha_t = function(t) 1,
                                        seed, n_cores, path = "results/assumption_curvature/",
                                        save = TRUE, return = FALSE){
  p <- d*(d-1)/2
  
  par_vec <- c()
  for(t in 1:(d-1)){
    par_vec <- c(par_vec, rep(par(t), (d-t)))
  }
  
  if(dvine){
    RVM <- D2RVine(1:d, family = rep(1, p), par = par_vec)
  } else {
    RVM <- C2RVine(1:d, family = rep(1, p), par = par_vec)
  }
  
  par_true <- RVM$par
  
  Matrix <- RVM$Matrix
  Matrix <- ToLowerTri(Matrix)
  M <- Matrix
  M <- reorderRVineMatrix(M)
  
  MaxMat <- createMaxMat(M)
  CondDistr <- neededCondDistr(M)
  
  pair_cor <- RVinePcor2cor(RVM)
  
  # generate rho_1, ..., rho_K to approximate the suppremum
  set.seed(seed)
  t_vec <- rep(1, (d - 1))
  for(k in 2:(d - 1)){
    t_vec <- c(t_vec, rep(k, (d - k)))
  }
  w_vec <- sapply(t_vec, alpha_t)
  par_list <- lapply(1:K, function(x) RVM$par + par_vec_to_mat(error_const*w_vec*ifelse(rbinom(p, 1, 0.5) == 1, 1, -1)))
  par_list <- append(list(par_true), par_list)
  
  # compute E[phi(U; theta)]
  res <- parallel::mclapply(1:(K+1), function(i) compute_phi_expectation(M, MaxMat, CondDistr, pair_cor, par_list[[i]]),
                            mc.cores = n_cores)
  
  E_phi_true <- max(abs(res[[1]])) # E(phi) at true parameter: should be 0 or close to 0
  
  res <- res[2:(K+1)] # remove phi at true parameter
  
  # matrix with one column for each theta, p x K dimensional, normalize by error_const
  res_matrix <- do.call(cbind, res)/error_const
  res_matrix <- res_matrix/matrix(rep(w_vec, K), ncol = K, byrow = FALSE) # normalize by alpha
  
  # matrix with sign(delta)
  par_list_mat <- ((lapply(par_list, par_mat_to_vec) %>% do.call(cbind, .))[, 2:(K+1)] - par_mat_to_vec(par_true)) %>% sign
  
  res_product <- res_matrix*par_list_mat
  max_E_phi <- max(res_product) # should be negative an bounded away from 0
  
  
  df <- data.frame(d = d, p = d*(d-1)/2, K = K, par = deparse1(par), 
                   dvine = dvine, error_const = error_const, seed = seed, 
                   name = name, alpha_t = deparse1(alpha_t),
                   E_phi_true = E_phi_true, max_E_phi = max_E_phi)
  
  result <- list(RVM = RVM, par_list = par_list, df = df)
  
  if(save){
    saveRDS(result, paste0(path, name, ".rds"))
  }
  
  
  if(return){
    result
  }
}

# compute expectation of phi at a given parameter
# for a gaussian vine
compute_phi_expectation <- function(M, MaxMat, CondDistr, pair_cor, par_eval){
  
  d <- nrow(M)
  p <- d*(d-1)/2
  
  expectation_list <- get_expectation_all(par_eval, M, MaxMat, CondDistr, pair_cor)
  
  phi <- c()
  for(i in d:2){
    for(j in (i-1):1){
      rho <- par_eval[i, j]
      
      x1_squared <- expectation_list[[d-i+1]][[j]][[1]][i, j] # 1 = direct*direct. x1 is always direct
      
      # is x2 direct or indirect?
      m <- MaxMat[i, j]
      if(m == M[i, j]){ #direct
        x2_squared <- expectation_list[[d-i+1]][[d-m+1]][[1]][i, d-m+1]
        x1x2 <- expectation_list[[d-i+1]][[j]][[1]][i, d-m+1] 
      }else{ #indirect
        x2_squared <- expectation_list[[d-i+1]][[d-m+1]][[4]][i, d-m+1] # 4 = indirect*indirect
        x1x2 <- expectation_list[[d-i+1]][[j]][[2]][i, d-m+1] # 2 = direct(=x1)*indirect(=x2)
      }
      
      res <- expectation_phi(rho, x1_squared, x2_squared, x1x2)
      phi <- c(phi, res)
    }
  }
  
  phi
  
}

# expectation of phi for a Gaussian vine, given E(X1^2), E(X2^2), E(X1*X2)
expectation_phi <- function(rho, x1_squared, x2_squared, x1x2){
  rho/(1-rho^2) - rho/(1-rho^2)^2*(x1_squared + x2_squared) + (1+rho^2)/(1-rho^2)^2*x1x2
}

# compute all expectations of products of two conditional variables in a Gaussian vine (with X_{a|D} = Phi^{-1}(U_{{a|D}})), 
# evaluated at a given parameter par_eval
# pair_cor (matrix of pairwise correlations) contains all necessary information from the true parameter
# the expectations are computed recursively using the provided parameter par_eval and the h-function of a Gaussian copula
# returns a list of (d-1) lists. The k-th list contains (d-k+1) lists, 
# containing the expectations of the products of all conditional X_{a|D} belonging to a parameter of the (k-1)-th tree (0-tree = original X_j) with all conditional X_{b|D'}
# for each rho_{a,b; D} of the (k-1)-th tree, there are four lists: 
# direct*direct, direct*indirect (expectations of X_{a|D,b} (direct) with all direct/indirect X_{c|D'})
# indirect*direct, indirect*indirect (expectations of X_{b|D,a} (indirect) with all direct/indirect X_{c|D'})
# the first list, corresponding to expectations with the original X_j contains only two lists
get_expectation_all <- function(par_eval, M, MaxMat, CondDistr, pair_cor){
  
  d <- ncol(M)
  
  expectation_list <- list()
  expectation_list0 <- list()
  for(k in 1:d){
    expectation_list0[[k]] <- get_expectation(d-k+1, par_eval, M, MaxMat, CondDistr, pair_cor) 
  }
  
  # 1. list: expectations with original X_j
  expectation_list[[1]] <- expectation_list0 
  
  for(l in 2:(d-1)){ 
    
    list_l <- list()
    for(k in 1:(d-l+1)){ 
      
      # direct-direct matrix and direct-indirect matrix
      direct1 <- indirect1 <- matrix(NA, ncol = d, nrow = d)
      if(CondDistr$direct[d-l+1, k]){
        
        direct1[d,] <- sapply(1:d, function(m) expectation_list[[1]][[m]]$direct[d-l+1, k]) 
        
        
        for(j in (d-1):2){ # rows, starting at the bottom
          for(i in 1:j){ # columns
            
            x_a <- direct1[(j + 1), i]
            m <- MaxMat[(j + 1), i]
            x_b <- if(m == M[(j + 1), i]){
              direct1[(j + 1), (d - m + 1)]
            } else{
              indirect1[(j + 1), (d - m + 1)]
            }
            
            if(CondDistr$direct[j, i]){ 
              direct1[j, i] <- h(x_a, x_b, par_eval[(j + 1), i])
            }
            if(CondDistr$indirect[j, i]){
              indirect1[j, i] <- h(x_b, x_a, par_eval[(j + 1), i])
            }
          }
        }
        
      }
      
      # indirect-direct matrix and indirect-indirect matrix
      direct2 <- indirect2 <- matrix(NA, ncol = d, nrow = d)
      if(CondDistr$indirect[d-l+1, k]){
        
        direct2[d,] <- sapply(1:d, function(m) expectation_list[[1]][[m]]$indirect[d-l+1, k]) # das unten mit indirect
        
        
        for(j in (d-1):2){ 
          for(i in 1:j){ 
            
            x_a <- direct2[(j + 1), i]
            m <- MaxMat[(j + 1), i]
            x_b <- if(m == M[(j + 1), i]){
              direct2[(j + 1), (d - m + 1)]
            } else{
              indirect2[(j + 1), (d - m + 1)]
            }
            
            if(CondDistr$direct[j, i]){ 
              direct2[j, i] <- h(x_a, x_b, par_eval[(j + 1), i])
            }
            if(CondDistr$indirect[j, i]){
              indirect2[j, i] <- h(x_b, x_a, par_eval[(j + 1), i])
            }
          }
        }
      }
      
      list_l[[k]] <- list(direct1, indirect1, direct2, indirect2) 
      
    }
    expectation_list[[l]] <- list_l
  }
  
  return(expectation_list)
}

# 1. list: expectations with original X_j
# computed recursively from the pairwise correlations
get_expectation <- function(k, par_eval, M, MaxMat, CondDistr, pair_cor){
  d <- ncol(M)
  
  direct <- indirect <- matrix(NA, ncol = d, nrow = d)
  
  direct[d,] <- pair_cor[k,d:1]
  
  for(j in (d-1):2){ 
    for(i in 1:j){ 
      
      x_a <- direct[(j + 1), i]
      m <- MaxMat[(j + 1), i]
      x_b <- if(m == M[(j + 1), i]){
        direct[(j + 1), (d - m + 1)]
      } else{
        indirect[(j + 1), (d - m + 1)]
      }
      
      if(CondDistr$direct[j, i]){
        direct[j, i] <- h(x_a, x_b, par_eval[(j + 1), i])
      }
      if(CondDistr$indirect[j, i]){
        indirect[j, i] <- h(x_b, x_a, par_eval[(j + 1), i])
      }
    }
  }
  
  return(list(direct = direct, indirect = indirect))
}



##### OLD (empirical estimation of expectation) ####
compute_phi_wrapper <- function(u, i, Matrix, par_list, n_cores,
                                temp_folder){
  
  phi_mat <- parallel::mclapply(par_list, function(theta) compute_phi_analytic(u, Matrix, theta, h_G, log_c_deriv_G),
                                mc.cores = n_cores) %>% do.call(cbind, .)
  saveRDS(phi_mat, paste0(temp_folder, "/phi_mat_", i, ".rds"))
  
}

verify_assumption_curvature_old <- function(d, n, K, par, dvine, name, error_const, B = 200, bootstr = FALSE,
                                        alpha_t = function(t) 1,
                                        seed, n_cores, temp_folder = "temp_assumption", 
                                        path = "results/assumption_curvature/",save = TRUE, return = FALSE){
  
  
  if (!file.exists(temp_folder)){
    dir.create(temp_folder)
  } else{ # delete all files in temp_folder
    unlink(paste0(temp_folder, "/*"))
  }
  
  p <- d*(d-1)/2
  
  par_vec <- c()
  for(t in 1:(d-1)){
      par_vec <- c(par_vec, rep(par(t), (d-t)))
  }
  
  if(dvine){
    RVM <- D2RVine(1:d, family = rep(1, p), par = par_vec)
  } else {
    RVM <- C2RVine(1:d, family = rep(1, p), par = par_vec)
  }
  
  par_true <- RVM$par
  
  # simulate data
  set.seed(seed)
  sim <- RVineSim(n, RVM)
  
  # generate theta_1, ..., theta_K
  t_vec <- rep(1, (d - 1))
  for(k in 2:(d - 1)){
    t_vec <- c(t_vec, rep(k, (d - k)))
  }
  w_vec <- sapply(t_vec, alpha_t)
  par_list <- lapply(1:K, function(x) RVM$par + par_vec_to_mat(error_const*w_vec*ifelse(rbinom(p, 1, 0.5) == 1, 1, -1)))
  par_list <- append(list(par_true), par_list)
  
  Matrix <- RVM$Matrix
  
  lapply(1:n, function(i) compute_phi_wrapper(sim[i,], i, Matrix, par_list, n_cores,
                                              temp_folder))
  
  par_list_mat <- ((lapply(par_list, par_mat_to_vec) %>% do.call(cbind, .))[, 2:(K+1)] - par_mat_to_vec(par_true)) %>% sign
  par_list_mat <- par_list_mat/(error_const* matrix(rep(w_vec, K), ncol = K, byrow = FALSE))
    
  
  df <- data.frame(d = d, p = d*(d-1)/2, K = K, n = n, par = deparse1(par), 
                   dvine = dvine, error_const = error_const, seed = seed, 
                   name = name, alpha_t = deparse1(alpha_t))
  
  if(!bootstr){
    E_phi_mat <- matrix(0, nrow = p, ncol = K + 1)
    for(i in 1:n){
      E_phi_mat <- E_phi_mat + readRDS(paste0(temp_folder, "/phi_mat_", i, ".rds"))/n
      }
    
    diff_mat <- (E_phi_mat - E_phi_mat[,1])[, 2:(K+1)] 
    res <- diff_mat*par_list_mat
    
    result <- list(res = res, RVM = RVM, par_list = par_list, sim = sim, df = df)
  } else{
    
    res_array <- array(dim = c(p, K, n))
    
    for(i in 1:n){
      res_i <- readRDS(paste0(temp_folder, "/phi_mat_", i, ".rds"))
      res_array[, , i] <- ((res_i - res_i[, 1])[, 2:(K+1)] )*par_list_mat
    }
    
    res <- mean_matrix <- apply(res_array, c(1,2), mean) 
    
    bootstrap_array <- array(dim = c(p, K, B)) # B = number of bootstrap samples
    
    for(b in 1:B){
      ind_b <- sample(1:n, n, replace = TRUE)
      bootstrap_array[, , b] <- apply(res_array[, , ind_b], c(1,2), mean)
    }
    
    vec_bootstrap <- apply(bootstrap_array, 3, max)
    
    result <- list(res = res, vec_bootstrap = vec_bootstrap, RVM = RVM, par_list = par_list, sim = sim, df = df)
  }
  
  
  
  if(save){
    saveRDS(result, paste0(path, name, ".rds"))
  }
  
  
  if(return){
    result
  }
  
}

compute_phi_analytic <- function(u, Matrix, par_eval, h, log_c_deriv){
  
  d <- nrow(Matrix)
  p <- d*(d-1)/2
  
  ## reorder matrix to natural order
  Matrix <- ToLowerTri(Matrix)
  M <- Matrix
  M <- reorderRVineMatrix(M)
  
  ## create matrices required for selection of h-functions
  MaxMat <- createMaxMat(M)
  CondDistr <- neededCondDistr(M) 
  
  list_u <- get_cond_u(u, h, par_eval, CondDistr, M, MaxMat)
  
  phi <- c()
  for(i in d:2){
    for(j in (i-1):1){
      theta <- par_eval[i, j]
      u1 <- list_u[[1]][i, j] 
      
      m <- MaxMat[i, j]
      if(m == M[i, j]){ #direct
        u2 <- list_u[[1]][i, d - m + 1]
      }else{ #indirect
        u2 <- list_u[[2]][i, d - m + 1]
      }
      phi <- c(phi, log_c_deriv(u1, u2, theta))
    }
  }
  
  phi
}


