source("functions/helpers.R")

compute_phi_wrapper <- function(u, i, Matrix, par_list, n_cores,
                                temp_folder){
  
  phi_mat <- parallel::mclapply(par_list, function(theta) compute_phi_analytic(u, Matrix, theta, h_G, log_c_deriv_G),
                                mc.cores = n_cores) %>% do.call(cbind, .)
  saveRDS(phi_mat, paste0(temp_folder, "/phi_mat_", i, ".rds"))
  
}

verify_assumption_curvature <- function(d, n, K, par, dvine, name, error_const, B = 200, bootstr = FALSE,
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


