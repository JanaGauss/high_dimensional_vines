source("functions/helpers.R")

vine_cop_par_matrix <- function(vc){
  
  params <- get_all_parameters(vc)
  
  d <- length(params) + 1
  par_mat <- matrix(0, nrow =d, ncol = d)
  for(i in 1:(d-1)){
    par_mat[d - i + 1, 1:(d - i)] <- unlist(params[[i]])
  }
  
  par_mat
}

vine_cop_par_matrix2 <- function(vc){
  
  params <- get_all_parameters(vc)
  
  d <- length(params) + 1
  par_mat <- par_mat2 <- matrix(0, nrow =d, ncol = d)
  for(i in 1:(d-1)){
    par_mat[d - i + 1, 1:(d - i)] <- unlist(lapply(params[[i]], function(x) x[1,]))
    par_mat2[d - i + 1, 1:(d - i)] <- unlist(lapply(params[[i]], function(x) x[2,]))
  }
  
  list(par_mat, par_mat2)
}

fam_vec <- c("gaussian", "t", "clayton", "gumbel", "frank", "joe", "xtd_gumbel")
names(fam_vec) <- as.character(c(1:6, 1004))

RVM2VinecopDist <- function(RVM){
  
  d <- nrow(RVM$Matrix)
  
  struct <- rvine_matrix(RVM$Matrix[d:1, ])
  pair_cops <- list()
  
  for(i in 1:(d - 1)){
    pair_cops[[i]] <- list()
    for(j in 1:(d - i)){
      fam <- fam_vec[as.character(RVM$family[(d - i + 1), j])]
      par <- RVM$par[(d - i + 1), j]
      if(fam == "t"){
        par2 <- RVM$par2[(d - i + 1), j]
        pair_cops[[i]][[j]] <- bicop_dist(family = fam, parameters = c(par, par2))
      } else{
        pair_cops[[i]][[j]] <- bicop_dist(family = fam, parameters = par)
      }
    }
  }
  vinecop_dist(pair_cops, struct)
}


sample_and_estim <- function(seed, n, vc, n_cores, family = "gaussian", 
                             margins = FALSE, dist_margins = "norm", itau = FALSE){
  set.seed(seed)
  print(seed)
  if(!margins){
    sim <- rvinecop(n, vc,
                    cores = n_cores)
  } else{
    vc_margins <- vine_dist(list(list(distr = dist_margins)), 
                            pair_copulas = vc$pair_copulas, structure = vc$structure)
    sim <- rvine(n, vc_margins, cores = n_cores) %>% pseudo_obs()
  }
  
  estim <- vinecop(sim, family_set = family,
                   structure = vc$structure, 
                   cores = n_cores, allow_rotations = FALSE, show_trace = FALSE, 
                   par_method = ifelse(itau, "itau", "mle"))
  estim %>% vine_cop_par_matrix() %>% par_mat_to_vec()
}


sample_and_estim2 <- function(seed, n, vc, n_cores, family = "t", 
                              margins = FALSE, dist_margins = "norm", itau = FALSE){
  set.seed(seed)
  print(seed)
  if(!margins){
    sim <- rvinecop(n, vc,
                    cores = n_cores)
  } else{
    vc_margins <- vine_dist(list(list(distr = dist_margins)), 
                            pair_copulas = vc$pair_copulas, structure = vc$structure)
    sim <- rvine(n, vc_margins, cores = n_cores) %>% pseudo_obs()
  }
  estim <- vinecop(sim, family_set = family,
                   structure = vc$structure, 
                   cores = n_cores, allow_rotations = FALSE, show_trace = FALSE, 
                   par_method = ifelse(itau, "itau", "mle"))
  list(par_mat_to_vec(vine_cop_par_matrix2(estim)[[1]]), par_mat_to_vec(vine_cop_par_matrix2(estim)[[2]]))
}

estim_joint_ML <- function(seed, n, RVM, res_seq_est){
  set.seed(seed)
  sim <- RVineSim(n, RVM) 
  estim <- RVineMLE(sim, RVM, start = res_seq_est)
  par_mat_to_vec(estim$RVM$par)
  
}

estim_joint_ML2 <- function(seed, n, RVM, res_seq_est, res_seq_est2){
  set.seed(seed)
  sim <- RVineSim(n, RVM) 
  estim <- RVineMLE(sim, RVM, start = res_seq_est, start2 = res_seq_est2, max.df = 6)
  list(par_mat_to_vec(estim$RVM$par), par_mat_to_vec(estim$RVM$par2))
}


estimate_param <- function(d = 10, n = 100, N = 2, bicop_fam = 1, par = function(t){0}, par2 = 4,
                           ktau = FALSE, two_arg = FALSE, dvine = TRUE,
                           seed = 1234, name = "test", save = FALSE, return = TRUE, n_cores = 1,
                           path = "results/", folder = NA,
                           RVM = NULL, jointML = FALSE, margins = FALSE, dist_margins = "norm", itau = FALSE){
  
  if(is.null(RVM)){
    par_vec <- c()
    if(!two_arg){
      for(t in 1:(d-1)){
        par_vec <- c(par_vec, rep(par(t), (d-t)))
        }
      } else{
        for(t in 1:(d-1)){
          par_vec <- c(par_vec, rep(par(t, d), (d-t)))
          }
      }
    
    if(bicop_fam %in% c(3, 5, 6)){
      par_vec <- sapply(par_vec, function(x) max(c(x, 1e-5)))
    }
    
    if(ktau){
      par_vec <- BiCopTau2Par(bicop_fam, par_vec)
    }
    
    
    if(bicop_fam != 2){
      if(dvine){
      RVM <- D2RVine(1:d, family = rep(bicop_fam, d*(d-1)/2), par = par_vec)
      } else {
        RVM <- C2RVine(1:d, family = rep(bicop_fam, d*(d-1)/2), par = par_vec)
      }
    }else{
      if(dvine){
        RVM <- D2RVine(1:d, family = rep(bicop_fam, d*(d-1)/2), par = par_vec, par2 = rep(par2, length(par_vec)))
      } else {
        RVM <- C2RVine(1:d, family = rep(bicop_fam, d*(d-1)/2), par = par_vec, par2 = rep(par2, length(par_vec)))
      }
      
      theta_true2 <- RVM$par2
    }
  } 
  
  vc <- RVM2VinecopDist(RVM)
  
  family <- fam_vec[as.character(bicop_fam)]
  
  theta_true <- RVM$par
  
  # simulate data
  set.seed(seed)
  seed_list <- as.integer(sample(1:1e7, N, replace = FALSE))
  if(bicop_fam != 2){
    
    res_list <- lapply(seed_list, function(s) sample_and_estim(s, n, vc, n_cores, family,
                                                               margins, dist_margins, itau))
    theta_hat_matrix <- do.call(rbind, res_list)
    
    df <- data.frame(d = d, p = d*(d-1)/2, n = n, N = N, par = deparse1(par), 
                   bicop_fam = bicop_fam,
                   dvine = dvine, seed = seed,
                   name = name)
    result <- list(theta_hat_matrix = theta_hat_matrix,
                 theta_true = par_mat_to_vec(RVM$par), RVM = RVM,
                 df = df) 
    
    
    if(jointML){
      res_list <- parallel::mclapply(1:N, function(i) estim_joint_ML(seed_list[i], n, RVM, 
                                                                     par_vec_to_mat(result$theta_hat_matrix[i, ])),
                                     mc.cores = n_cores)
      theta_hat_joint_matrix <- do.call(rbind, res_list)
    
      result <- append(result, list(theta_hat_joint_matrix = theta_hat_joint_matrix), after = 1)
      }
    
    
  } else{
    res_list <- lapply(seed_list, function(s) sample_and_estim2(s, n, vc, n_cores, family = "t",
                                                                margins, dist_margins, itau))
    theta_hat_matrix <- sapply(1:N, function(i) res_list[[i]][1]) %>% do.call(rbind, .)
    theta_hat_matrix2 <- sapply(1:N, function(i) res_list[[i]][2]) %>% do.call(rbind, .)
    
    df <- data.frame(d = d, p = d*(d-1)/2, n = n, N = N, par = deparse1(par), par2 = par2,
                     bicop_fam = bicop_fam,
                     dvine = dvine, seed = seed,
                     name = name)
    result <- list(theta_hat_matrix = theta_hat_matrix, theta_hat_matrix2 = theta_hat_matrix2,
                   theta_true = par_mat_to_vec(RVM$par), theta_true2 = par_mat_to_vec(RVM$par2),
                   RVM = RVM,
                   df = df)
    
    if(jointML){
      res_list <- parallel::mclapply(1:N, function(i) estim_joint_ML2(seed_list[i], n, RVM, 
                                                                     par_vec_to_mat(result$theta_hat_matrix[i, ]),
                                                                     par_vec_to_mat(result$theta_hat_matrix2[i,])),
                                     mc.cores = n_cores)
      theta_hat_joint_matrix <- sapply(1:N, function(i) res_list[[i]][1]) %>% do.call(rbind, .)
      theta_hat_joint_matrix2 <- sapply(1:N, function(i) res_list[[i]][2]) %>% do.call(rbind, .)
      result <- append(result, list(theta_hat_joint_matrix = theta_hat_joint_matrix), after = 1)
      result <- append(result, list(theta_hat_joint_matrix = theta_hat_joint_matrix2), after = 3)
    }
  }

  
  if(save){
    if(is.na(folder)){
      if(margins){
      saveRDS(result, paste0(path, "estimate_params_margins", ifelse(itau, "_itau", ""), "/", name, ".rds"))
    }else{
      saveRDS(result, paste0(path, "estimate_params", ifelse(itau, "_itau", ""), "/", name, ".rds"))
    }
    } else{
      saveRDS(result, paste0(path, folder, name, ".rds"))
    }
    
    
  }

  if(return){
    result
  }
}

