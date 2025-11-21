
#### create data frames from simulation results ####
create_data <- function(path = "results/estimate_params/", 
                        substr_files = NULL, par_order = NULL, par_2 = FALSE){
  
  files <- list.files(path)
  
  if(!is.null(substr_files)){
    files <- str_subset(files, substr_files)
  }
  
  res_df <- data.frame()
  
  for(f in files){
    
    result <- readRDS(paste0(path, f))
    
    diff_matrix <- result$theta_hat_matrix - t(replicate(nrow(result$theta_hat_matrix), result$theta_true))
    
    norm_vec <- apply(diff_matrix, 1, function(x) norm(x, "2"))
    sum_vec <- rowSums(diff_matrix)
    max_norm_vec <- apply(diff_matrix, 1, function(x) max(abs(x)))
    
    p <- length(result$theta_true)
    last_param_vec <- result$theta_hat_matrix[, p] - result$theta_true[p]
    
    if(!par_2){
      
      df_temp <- data.frame(d = result$df$d, p = result$df$p, n = result$df$n, par = result$df$par, 
                            family = result$df$bicop_fam, dvine = result$df$dvine,
                            norm_error = norm_vec, sum_error = sum_vec, max_norm_error = max_norm_vec, last_param_error = last_param_vec)
      
      res_df <- rbind(res_df, df_temp)
      
      
    } else{
      
      diff_matrix2 <- result$theta_hat_matrix2 - t(replicate(nrow(result$theta_hat_matrix2), result$theta_true2))
      
      norm_vec2 <- apply(diff_matrix2, 1, function(x) norm(x, "2"))
      sum_vec2 <- rowSums(diff_matrix2)
      last_param_vec2 <- result$theta_hat_matrix2[, p] - result$theta_true2[p]
      max_norm_vec2 <- apply(diff_matrix2, 1, function(x) max(abs(x)))
      
      norm_vec_all <- apply(cbind(diff_matrix, diff_matrix2), 1, function(x) norm(x, "2"))
      
      df_temp <- data.frame(d = result$df$d, p = result$df$p, n = result$df$n, par = result$df$par, 
                            par2 = result$df$par2,
                            family = result$df$bicop_fam, dvine = result$df$dvine,
                            norm_error = norm_vec, sum_error = sum_vec, last_param_error = last_param_vec, max_norm_error = max_norm_vec,
                            norm_error2 = norm_vec2, sum_error2 = sum_vec2, last_param_error2 = last_param_vec2, max_norm_error2 = max_norm_vec2,
                            norm_error_all = norm_vec_all)
      
      res_df <- rbind(res_df, df_temp)
      
    }
    
    
  }
  
  if(!is.null(par_order)){
    res_df$par <- factor(res_df$par, levels = unique(res_df$par)[par_order])
  }
  
  return(res_df)
}

#### inter/extrapolate data for verification of consistency ####
create_dat_interpol <- function(data, 
                        substr_files = NULL, par_order = NULL, par_2 = FALSE){
  
  f_1 <- function(d) d*25
  f_2 <- function(d){
    d^2*0.125
  }
  f_3 <- function(d){
    (d-50)^3*0.003
  }
  if(max(data$d) <= 100){
    f_3 <- function(d){
      d^3 * 0.005
    }
  }
  
  f_list <- list(f_1, f_2, f_3)
  
  dat_interpol <- data %>% select(c(d, p, par, family, dvine)) %>% 
    unique() %>% crossing(f_n = 1:length(f_list))
  
  dat_interpol$n_new <- dat_interpol$error_n_new  <- NA
  
  dat_mean <- data %>% group_by(n, d, p, par, family, dvine) %>% 
    summarise(log_error = mean(log(max_norm_error)))
  
  
  for(i in 1:nrow(dat_interpol)){
    d <- dat_interpol$d[i]
    f_n_new <- f_list[[dat_interpol$f_n[i]]]
    dat_interpol$n_new[i] <- n_new <- f_n_new(d) %>% round()
    if(n_new <= 0){
      dat_interpol$error_n_new[i] <- NA
      next
    }
    
    dat_small <- data %>% filter(d == dat_interpol$d[i], par == dat_interpol$par[i],
                                dvine == dat_interpol$dvine[i], 
                                family == dat_interpol$family[i])
    dat_mean_small <- dat_mean %>% filter(d == dat_interpol$d[i], par == dat_interpol$par[i],
                                          dvine == dat_interpol$dvine[i], 
                                          family == dat_interpol$family[i])
    
    n_values <- unique(dat_small$n) %>% sort()
    
    # interpolation or extrapolation
    if(n_new < min(n_values) | n_new > max(n_values)){ # extrapolation
      if(n_new < min(n_values)){
        dat_lm <- dat_small %>% filter(n %in% n_values[c(1,2)])
      } else{
        dat_lm <- dat_small %>% filter(n %in% n_values[c(length(n_values), length(n_values) - 1)])
      }
      
      lm <- lm(log(max_norm_error) ~ log(n), data = dat_lm)
      
      intercept <- lm$coefficients[1]
      beta <- lm$coefficients[2]
      log_error_n_new <- intercept + beta*log(n_new)
      dat_interpol$error_n_new[i] <- exp(log_error_n_new)
      
    } else { # interpolation
      interpol_fun <- approxfun(x = log(dat_mean_small$n), y = dat_mean_small$log_error)
      dat_interpol$error_n_new[i] <- exp(interpol_fun(log(n_new)))
    }
    
  }
  
  dat_interpol
}

#### compute parameter matrix from vector and vice versa ####
par_vec_to_mat <- function(par_vec){
  d <- as.integer((1+sqrt(1+8*length(par_vec)))/2)
  
  par_mat <- matrix(0, nrow = d, ncol = d)
  l <- d - 1
  for(i in d:2){
    par_mat[i, l:1] <- par_vec[1:l]
    par_vec <- par_vec[l+1:length(par_vec)]
    l <- l - 1
  }
  par_mat
}

par_mat_to_vec <- function(par_mat){
  
  d <- nrow(par_mat)
  par_vec <- c()
  l <- d - 1
  for(i in d:2){
    par_vec <- append(par_vec, par_mat[i, l:1])
    l <- l - 1
  }
  par_vec
}


#### compute necessary conditional u values ####
get_cond_u <- function(u, h, par_eval, CondDistr, M, MaxMat){
  d <- ncol(M)
  
  direct <- indirect <- matrix(NA, ncol = d, nrow = d)
  
  direct[d,] <- u[d:1]
  
  for(j in (d-1):2){ # rows, from bottom to top
    for(i in 1:j){ # columns
      
      u_a <- direct[(j + 1), i] 
      
      m <- MaxMat[(j + 1), i]
      u_b <- if(m == M[(j + 1), i]){
        direct[(j + 1), (d - m + 1)]
      } else{
        indirect[(j + 1), (d - m + 1)]
      }
      
      if(CondDistr$direct[j, i]){
        direct[j, i] <- h(u_a, u_b, par_eval[(j + 1), i])
      }
      if(CondDistr$indirect[j, i]){
        indirect[j, i] <- h(u_b, u_a, par_eval[(j + 1), i])
      }
    }
  }
  return(list(direct = direct, indirect = indirect))
}


#### get "paths" to compute derivative of phi ####
get_paths <- function(M, MaxMat){
  
  d <- ncol(M)
  
  n_paths <- 2^(d-2)
  paths_params <- paths_direct <- matrix(data = NA, ncol = (d-1), nrow = n_paths)

  paths_params[, 1] <- 1 
  
  paths_direct[, 1] <- c(rep(TRUE, n_paths/2), rep(FALSE, n_paths/2))
  
  n_paths_distinct <- 1
  for(j in 2:ncol(paths_params)){
    
    n <- n_paths/(n_paths_distinct*2) 
    
    for(k in 1:n_paths_distinct){
      paths_params[((k-1)*n_paths/n_paths_distinct + 1):((2*k-1)*n_paths/(n_paths_distinct*2)), j] <- paths_params[((k-1)*n_paths/n_paths_distinct + 1), (j-1)] 
      paths_direct[((k-1)*n_paths/n_paths_distinct + 1):((2*k-1)*n_paths/(n_paths_distinct*2)), j] <- TRUE
      
      m <- MaxMat[j, paths_params[((k-1)*n_paths/n_paths_distinct + 1), (j-1)]]
      paths_params[((2*k-1)*n_paths/(n_paths_distinct*2) + 1):(k*n_paths/n_paths_distinct), j] <- (d - m + 1)
      paths_direct[((2*k-1)*n_paths/(n_paths_distinct*2) + 1):(k*n_paths/n_paths_distinct), j] <- (m == M[j, paths_params[((k-1)*n_paths/n_paths_distinct + 1), (j-1)]])
      
    }
    
    n_paths_distinct <- n_paths_distinct*2
  }

  paths_params <- cbind(paths_params, rep(NA, n_paths))
  for(i in 1:n_paths){
    m <- MaxMat[d, paths_params[i, (d-1)]]
    paths_params[i, d] <- (d - m + 1)
  }
  
  paths_h <- matrix(data = NA, nrow = nrow(paths_params), ncol = (ncol(paths_params) - 2)) 
  
  for(j in 2:ncol(paths_h)){
    for(i in 1:nrow(paths_h)){
      if(paths_direct[i, j]){
        paths_h[i, j] <- (paths_params[i, j] == paths_params[i, (j+1)])
      } else{
        paths_h[i, j] <- !(paths_params[i, j] == paths_params[i, (j+1)])
      }
      
    }
  }
  
  return(list(paths_params = paths_params, paths_direct = paths_direct, paths_h = paths_h))
}

#### compute matrix nabla phi ####
nabla_phi <- function(Matrix, MaxMat, M, list_u, paths_list, score_deriv, score_deriv_u, h1, h2, h_d, par_eval){
  
  d <- ncol(Matrix)
  H <- matrix(0, nrow = d*(d-1)/2, ncol = d*(d-1)/2)
  
  # diagonal of nabla_phi
  diag <- c()
  for(x in d:2){
    for(y in (x-1):1){
      diag <- append(diag, get_nabla_phi_diag(x, y, par_eval[x, y], score_deriv, list_u, MaxMat, M))
    }
  }
  diag(H) <- diag
  
  i <- nrow(H)
  j <- 1
  for(x in 2:(d-1)){
    for(y in 1:(x-1)){
      for(k in d:(x+1)){
        t <- d - k + 1 # tree
        derivs <- get_derivs_xyt(x, y, t, list_u, paths_list, score_deriv_u, h1, h2, h_d, par_eval, M, MaxMat)
        H[i, j:(j+length(derivs) - 1)] <- derivs[length(derivs):1]
        j <- j + length(derivs)
      }
      i <- i-1
      j <- 1
    }
  }
  
  H
}

# Diagonal of nabla_phi
get_nabla_phi_diag <- function(x, y, theta, score_deriv, list_u, MaxMat, M){
  d <- ncol(M)
  
  u1 <- list_u[[1]][x, y] 
  
  m <- MaxMat[x,y]
  if(m == M[x,y]){ #direct
    u2 <- list_u[[1]][x, d - m + 1]
  }else{ #indirect
    u2 <- list_u[[2]][x, d - m + 1]
  }
  
  res <- score_deriv(u1, u2, theta)
}

# derivatives of estimating equation of parameter (x,y) w.r.t. all parameters of t-th tree
get_derivs_xyt <- function(x, y, t, list_u, paths_list, score_deriv_u, h1, h2, h_d, par_eval, M, MaxMat){
  
  d <- ncol(M)
  paths_params <- paths_list$paths_params
  paths_direct <- paths_list$paths_direct
  paths_h <- as.matrix(paths_list$paths_h[, -1])
  
  ## 1. restrict paths_params, paths_h, paths_direct to necessary rows/columns
  ind <- which(paths_params[, (x - 1)] == y) 
  
  paths_params <- paths_params[ind, (x-1):(ncol(paths_params) - t + 1)]
  paths_direct <- paths_direct[ind, (x-1):(ncol(paths_direct) - t + 1)]
  paths_h <- if((x-1) <= (ncol(paths_h) - t + 1)){
    as.matrix(paths_h[ind, (x-1):(ncol(paths_h) - t + 1)])
  } else{
    matrix(data = 0, ncol = 0, nrow = 0)
  }
  
  
  
  # remove duplicated rows
  e <- 2^(t-1) # keep every k-th row
  ind2 <- e*(1:(2^(d-x)/e))
  paths_params <- paths_params[ind2,]
  paths_direct <- paths_direct[ind2,]
  if(ncol(paths_h) > 0){
    paths_h <- as.matrix(paths_h[ind2,])
  }
  
  paths_direct[, 1] <- c(rep(TRUE, nrow(paths_direct)/2), rep(FALSE, nrow(paths_direct)/2))
  # first half of paths always belongs to first argument in density c, second half belongs to second argument
  # -> first column of paths_direct must indicate this
  
  res <- rep(0, d-t) # d-t parameters in t-th tree
  
  ## 2. compute derivative for each row of paths_params and add to respective entry of res
  for(i in 1:nrow(paths_params)){ 
    
    k <- paths_params[i, (ncol(paths_params) - 1)] # derivative w.r.t. k-th parameter of t-th tree
    
    p1 <- par_eval[x, y] # we compute derivative of estimating equ. of parameter p1
    p2 <- par_eval[nrow(par_eval) - t + 1, k] # derivative w.r.t. p2
    
    # derivative is product of several terms
    # first term: deriv. of log density of p1 w.r.t respective u
    # last term: deriv. of h function w.r.t p2
    
    # u/v value in log density of p1 
    u <- list_u[[1]][x,y] # log density of p1 always contains direct value [x,y]
    m <- MaxMat[x,y]
    y2 <- (d - m + 1) # [x, y2] = second value in log density of p1
    v <- ifelse(m == M[x,y], list_u[[1]][x, y2], list_u[[2]][x, y2]) # value is direct/indirect depending on M and MaxMat
    u1_1 <- ifelse(paths_direct[i, 1], u, v) # paths_direct indicates order of u and v
    u2_1 <- ifelse(paths_direct[i, 1], v, u)
    
    # u/v value in  h function of p2
    u <- list_u[[1]][d - t + 1, k] # h function of p2 always contains direct value [d - t + 1, k]
    m <- MaxMat[d - t + 1, k]
    y2 <- (d - m + 1) # [d - t + 1, y2] = second value in h function of p2
    v <- ifelse(m == M[d - t + 1, k], list_u[[1]][d - t + 1, y2], list_u[[2]][d - t + 1, y2])  # value is direct/indirect depending on M and MaxMat
    u1_2 <- ifelse(paths_direct[i, ncol(paths_direct)], u, v) # paths_direct indicates order of u and v
    u2_2 <- ifelse(paths_direct[i, ncol(paths_direct)], v, u)
    
    prod <- score_deriv_u(u1_1, u2_1, p1)*h_d(u1_2, u2_2, p2)
    # product of last and first term
    
    # product of derivatives of h functions
    if(ncol(paths_h) > 0){
      for(j in 1:ncol(paths_h)){ 
        row_j <- j + x
        col_j <- paths_params[i, (j+1)]
        par <- par_eval[row_j, col_j] # parameter in h function
        u <- list_u[[1]][row_j, col_j] # first value in h function
        m <- MaxMat[row_j, col_j]
        y2 <- (d - m + 1) # [row_j, y2] = second value in h function 
        v <- ifelse(m == M[row_j, col_j], list_u[[1]][row_j, y2], list_u[[2]][row_j, y2]) 
        
        u1 <- ifelse(paths_direct[i, j + 1], u, v)
        u2 <- ifelse(paths_direct[i, j + 1], v, u)
        
        u1
        u2
        
        if(paths_h[i, j]){ # deriv. of h function w.r.t. first or second parameter?
          prod <- prod*h1(u1, u2, par)
        }else{
          prod <- prod*h2(u1, u2, par)
        }
      }
    }
    
    res[k] <- res[k] + prod
    
  }
  
  return(res)
  
}

############# internal functions from VineCopula Package ##############
# https://rdrr.io/cran/VineCopula/src/R/RVineMatrix.R

## function that converts upper triagonal matrix to lower triagonal
ToLowerTri <- function(x) {
  ## only change matrix if not already lower triagonal
  if(all(x[lower.tri(x)] == 0)) {
    x[nrow(x):1, ncol(x):1]
  } else {
    x
  }
}

reorderRVineMatrix <- function(Matrix, oldOrder = NULL) {
  
  if (length(oldOrder) == 0) {
    oldOrder <- diag(Matrix)
  }
  O <- apply(t(1:nrow(Matrix)), 2, "==", Matrix)
  
  for (i in 1:nrow(Matrix)) {
    Matrix[O[, oldOrder[i]]] <- nrow(Matrix) - i + 1
  }
  return(Matrix)
}

createMaxMat <- function(Matrix) {
  
  if (dim(Matrix)[1] != dim(Matrix)[2])
    stop("Structure matrix has to be quadratic.")
  
  MaxMat <- reorderRVineMatrix(Matrix)
  
  n <- nrow(MaxMat)
  
  for (j in 1:(n - 1)) {
    for (i in (n - 1):j) {
      MaxMat[i, j] <- max(MaxMat[i:(i + 1), j])
    }
  }
  
  tMaxMat <- MaxMat
  tMaxMat[is.na(tMaxMat)] <- 0
  
  oldSort <- diag(Matrix)
  oldSort <- oldSort[n:1]
  
  for (i in 1:n) {
    MaxMat[tMaxMat == i] <- oldSort[i]
  }
  
  return(MaxMat)
}

neededCondDistr <- function(Vine) {
  
  if (dim(Vine)[1] != dim(Vine)[2])
    stop("Structure matrix has to be quadratic.")
  
  Vine <- reorderRVineMatrix(Vine)
  
  MaxMat <- createMaxMat(Vine)
  
  d <- nrow(Vine)
  
  M <- list()
  M$direct <- matrix(FALSE, d, d)
  M$indirect <- matrix(FALSE, d, d)
  
  M$direct[2:d, 1] <- TRUE
  
  for (i in 2:(d - 1)) {
    v <- d - i + 1
    
    bw <- as.matrix(MaxMat[i:d, 1:(i - 1)]) == v
    
    direct <- Vine[i:d, 1:(i - 1)] == v
    
    M$indirect[i:d, i] <- apply(as.matrix(bw & (!direct)), 1, any)
    
    M$direct[i:d, i] <- TRUE
    
    M$direct[i, i] <- any(as.matrix(bw)[1, ] & as.matrix(direct)[1, ])
  }
  
  return(M)
}


#### functions for Gaussian Vines ####
log_c_deriv_G <- function(u1, u2, theta){
  theta/(1-theta^2) - theta/(1-theta^2)^2*(qnorm(u1)^2 + qnorm(u2)^2) + (1+theta^2)/(1-theta^2)^2*qnorm(u1)*qnorm(u2)
}

h_G <- function(u1, u2, theta){
  pnorm((qnorm(u1) - theta*qnorm(u2))/sqrt(1 - theta^2))
}
score_deriv_G <- function(u1, u2, theta){
  (1 + theta^2)/(1-theta^2)^2 - (1 + 3*theta^2)/(1-theta^2)^3*(qnorm(u1)^2 + qnorm(u2)^2) + 2*(3*theta + theta^3)/(1-theta^2)^3*qnorm(u1)*qnorm(u2)
}
score_deriv_u_G <- function(u1, u2, theta){
  ((1+theta^2)*qnorm(u2) - 2*theta*qnorm(u1))/(1-theta^2)^2/dnorm(qnorm(u1))
}
h1_G <- function(u1, u2, theta){
  dnorm((qnorm(u1) - theta*qnorm(u2))/sqrt(1 - theta^2)) /sqrt(1-theta^2)/dnorm(qnorm(u1))
}
h2_G <- function(u1, u2, theta){
  - dnorm((qnorm(u1) - theta*qnorm(u2))/sqrt(1 - theta^2)) *theta /sqrt(1-theta^2)/dnorm(qnorm(u2))
}
h_d_G <- function(u1, u2, theta){
  dnorm((qnorm(u1) - theta*qnorm(u2))/sqrt(1 - theta^2)) *(theta*qnorm(u1) - qnorm(u2))/(1-theta^2)^(3/2) 
}

