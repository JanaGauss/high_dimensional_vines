
source("packages.R")
source("functions/functions_assumption_curvature.R")

d_vector <- c(5, 8, 10, 15, 20, 30, 40, 50)
par <- c("function (t)  0", "function (t)  0.5/sqrt(t+1)", 
         "function (t)  1/(t+1)", "function (t)  0.5^t")
dvine <- c(TRUE, FALSE)
error_const <- 0.005

params_df <- expand.grid(d_vector, par, dvine, error_const)
names(params_df) <- c("d", "par", "dvine", "error_const")
params_df$n <- round(2000*log(params_df$d), 0)  
params_df$K <- 50

set.seed(12345)
params_df$seed <- sample(1:100000, nrow(params_df))


names_df <- data.frame(par = par,
                       name = c("constant_0", "harmonicRoot05", "harmonic",
                                "exponential"))
params_df <- left_join(params_df, names_df)
params_df$name <- paste(params_df$name, ifelse(params_df$dvine, "dvine", "cvine"), params_df$d, sep = "_")

params_df <- params_df[order(params_df$d), ]

n_cores <- parallel::detectCores()


lapply(split(params_df, seq(nrow(params_df))), 
       function(params) verify_assumption_curvature(d = as.numeric(params$d), n = as.numeric(params$n),
                                                    K = as.numeric(params$K), par = eval(parse(text = params$par)),
                                                    name = params$name, dvine = params$dvine, error_const = as.numeric(params$error_const),
                                                    seed = params$seed, n_cores = n_cores, save = TRUE, return = FALSE,
                                                    bootstr = FALSE))



#### modifications for 0.5/sqrt(t+1), D-vine ####
params_df_d <- params_df %>% filter(str_detect(name, "harmonicRoot05_")) %>% filter(dvine)

alpha_t <- function (t)  {     
  res <- 0     
  for (i in 1:t) {         
    res <- res + i^(-1.1)     
  }     
  res
}

alpha_list <- list(alpha_t, function(t) t, function(t) t^2, function(t) t^3) 
alpha_ind <- 1:length(alpha_list)

params_df_d <- merge(params_df_d, data.frame(alpha_ind = alpha_ind))

params_df_d$name <- paste0(params_df_d$name, "_", params_df_d$alpha_ind + 1)

params_df_d <- params_df_d %>% mutate(error_const = case_when(
  alpha_ind == 1 ~ 0.005,
  alpha_ind == 2 ~ 1e-7,
  alpha_ind == 3 ~ 1e-9,
  alpha_ind == 4 ~ 1e-10
))



lapply(split(params_df_d, seq(nrow(params_df_d))), 
       function(params) verify_assumption_curvature(d = as.numeric(params$d), n = as.numeric(params$n),
                                                    K = as.numeric(params$K), par = eval(parse(text = params$par)),
                                                    name = params$name, dvine = params$dvine, error_const = as.numeric(params$error_const),
                                                    seed = params$seed, n_cores = n_cores, save = TRUE, return = FALSE,
                                                    alpha_t = alpha_list[[params$alpha_ind]],
                                                    bootstr = FALSE))
