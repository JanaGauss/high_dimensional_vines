
fam_vec <- c("gaussian", "t", "clayton", "gumbel", "frank", "joe", "xtd_gumbel")
names(fam_vec) <- as.character(c(1:6, 1004))

sample_and_estim_trunc <- function(seed, n, vc, n_cores, family = "gaussian", tr_lvl,
                                   margins = FALSE, dist_margins = "norm"){
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
                   par_method = "mle",
                   trunc_lvl = tr_lvl)
  unlist(get_all_parameters(estim))
}


estimate_param_trunc <- function(d = 10, n = 100, N = 2, bicop_fam = 1, par = function(t){1/(t+1)}, tr_lvl = 2,
                                 dvine = FALSE, 
                                 margins = FALSE, dist_margins = "norm",
                                 seed = 1234, name = "test", save = FALSE, return = TRUE, n_cores = 1,
                                 path = "results/"){
  
  
  if(dvine){
    struct <- dvine_structure(d:1, trunc_lvl = tr_lvl)
  } else{
    struct <- cvine_structure(d:1, trunc_lvl = tr_lvl)
  }
  
  pair_cops <- list()
  family <- fam_vec[as.character(bicop_fam)]
  if(tr_lvl == Inf){
    tr_lvl <- d-1
  }
  for(i in 1:tr_lvl){
    pair_cops[[i]] <- list()
    param <- par(as.numeric(i))
    for(j in 1:(d - i)){
      pair_cops[[i]][[j]] <- bicop_dist(family = family, parameters = param)
    }
  }
  vc <- vinecop_dist(pair_cops, struct)
  theta_true <- unlist(get_all_parameters(vc))

  
  set.seed(seed)
  seed_list <- as.integer(sample(1:1e7, N, replace = FALSE))
  
  res_list <- lapply(seed_list, function(s) sample_and_estim_trunc(s, n, vc, n_cores, family, tr_lvl,
                                                                   margins, dist_margins))
  
  
  theta_hat_matrix <- do.call(rbind, res_list)
  
  df <- data.frame(d = d, p = d*(d-1)/2, n = n, N = N, par = deparse1(par), 
                   bicop_fam = bicop_fam,
                   dvine = dvine, seed = seed,
                   name = name)
  result <- list(theta_hat_matrix = theta_hat_matrix,
                 theta_true = theta_true, 
                 df = df)


  
  if(save){
    if(margins){
      saveRDS(result, paste0(path, "estimate_params_trunc_margins/", name, ".rds"))
    } else{
      saveRDS(result, paste0(path, "estimate_params_trunc/", name, ".rds"))
    }
    
  }
  
  if(return){
    result
  }
}

