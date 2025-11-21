source("functions/helpers.R")

estim_M <- function(sim, RVM, par_list, 
                    h, score_deriv, score_deriv_u, h1, h2, h_d,
                    alpha_t,
                    n_cores, temp_folder = "temp_M", name){
  
  d <- ncol(sim)
  p <- d*(d-1)/2
  n <- nrow(sim)
  
  ## reorder matrix to natural order
  Matrix <- ToLowerTri(RVM$Matrix)
  M <- Matrix
  M <- reorderRVineMatrix(M)
  
  ## create matrices required for selection of h-functions
  MaxMat <- createMaxMat(M)
  CondDistr <- neededCondDistr(M) 
  
  paths_list <- get_paths(M, MaxMat)
  
  
  parallel::mclapply(1:n, function(i) calc_nabla_phi_K(sim[i, ], i, name, par_list, 
                                                       paths_list, 
                                                       h, score_deriv, score_deriv_u, h1, h2, h_d, 
                                                       CondDistr, Matrix, M, MaxMat, temp_folder),
                       mc.cores = n_cores, mc.preschedule = FALSE)
  
  mean_matrix <- matrix(0, nrow = p, ncol = length(par_list))
  D_vec <- c()
  
  alpha_mat <- matrix(0, nrow = p, ncol = p)
  t_vec <- rep(1, (d - 1))
  for(k in 2:(d - 1)){
    t_vec <- c(t_vec, rep(k, (d - k)))
  }
  w_vec <- sapply(t_vec, alpha_t)
  for(k in 1:p){
    for(j in 1:p){
      alpha_mat[k, j] <- w_vec[j]/w_vec[k]
    }
  }
  alpha_tensor <- array(rep(alpha_mat, times = length(par_list)), dim = c(p, p, length(par_list)))
  
  for(i in 1:n){
    res_i <- readRDS(paste0(temp_folder, "/", name, "_", i, ".rds")) %>% abs() 
    # Hadamard product (elementwise multiplication) of matrices in res_i and alpha_mat
    res_i <- res_i * alpha_tensor
    
    sums_D <- apply(res_i, 3, rowSums) # dimension: p x (K+1)
    sums_squared <- sums_D^2 
    
    D_vec <- c(D_vec, max(sums_D))
    mean_matrix <- mean_matrix + sums_squared/n
  }
  
  return(list(mean_matrix = mean_matrix, D_vec = D_vec))
  
}

calc_nabla_phi_K <- function(u, i, name, par_list, paths_list, h, score_deriv, score_deriv_u, h1, h2, h_d, 
                             CondDistr, Matrix, M, MaxMat, temp_folder){
  
  d <- nrow(M)
  p <- d*(d-1)/2
  nabla_array <- array(dim = c(p, p, length(par_list)))
  
  for(k in 1:length(par_list)){
    theta <- par_list[[k]]
    cond_u <- get_cond_u(u, h, theta, CondDistr, M, MaxMat)
    
    nabla_array[, , k] <- nabla_phi(Matrix, MaxMat, M, cond_u, paths_list, score_deriv, score_deriv_u, h1, h2, h_d, theta)
  }
  
  saveRDS(nabla_array, paste0(temp_folder, "/", name, "_", i, ".rds"))
}

wrapper_estim_M <- function(d, n, K, seed, error_const, dvine, bicop_fam = 1, par = function(t){0},
                            name = "test", 
                            h, score_deriv, score_deriv_u, h1, h2, h_d,
                            alpha_t =  function(t) 1,
                            path = "results/estimate_M/", 
                            temp_folder = "temp_M",
                            n_cores, save = FALSE, return = TRUE){
  
  if (!file.exists(temp_folder)){
    dir.create(temp_folder)
  }
  
  p <- d*(d-1)/2
  
  par_vec <- c()
  for(t in 1:(d-1)){
      par_vec <- c(par_vec, rep(par(t), (d-t)))
  }
  
  if(dvine){
    RVM <- D2RVine(1:d, family = rep(bicop_fam, d*(d-1)/2), par = par_vec)
  } else {
    RVM <- C2RVine(1:d, family = rep(bicop_fam, d*(d-1)/2), par = par_vec)
  }
  
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
  par_list <- append(list(RVM$par), par_list)
  
  res <- estim_M(sim, RVM, par_list, 
                 h, score_deriv, score_deriv_u, h1, h2, h_d,
                 alpha_t,
                 n_cores, temp_folder, name)
  
  df <- data.frame(d = d, p = d*(d-1)/2, K = K, n = n, par = deparse1(par), 
                   dvine = dvine, error_const = error_const, seed = seed, 
                   name = name)
  
  result <- list(res = res, RVM = RVM, par_list = par_list, sim = sim, df = df)
  
  if(save){
    saveRDS(result, paste0(path, name, ".rds"))
  }
  
  
  if(return){
    result
  }
  
}
