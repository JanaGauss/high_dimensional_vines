
source("packages.R")
source("functions/M_n_estimation.R")

d_vector <- c(4, 6, 8, 10, 12, 15)
par_functions <- c("function (t)  0", 
                   "function (t)  1/(t+1)",
                   "function (t) 0.5^t",
                   "function (t) 0.5/sqrt(t+1)")
dvine <- c(TRUE, FALSE)
error_const <- 0.005

params_df <- expand.grid(d_vector, par_functions, dvine, error_const)
names(params_df) <- c("d", "par", "dvine", "error_const")

params_df$n <- round(200*log(params_df$d), 0) 
params_df$K <- 30

params_df <- params_df %>% filter(par != par_functions[4]) %>% rbind(filter(params_df, par == par_functions[4]))
# obtain order of settings used when generating seeds

set.seed(1234)
params_df$seed <- sample(1:100000, nrow(params_df))


names_df <- data.frame(par = par_functions,
                       name = c("constant_0","harmonic", "exponential", "harmonic_root05"))
params_df <- left_join(params_df, names_df)
params_df$name <- paste(params_df$name, ifelse(params_df$dvine, "dvine", "cvine"), params_df$d, sep = "_")

params_df <- params_df[order(params_df$d), ]

params_df_d <- params_df %>% filter(str_detect(name, "harmonic_root05_d")) 

params_df <- params_df %>% filter(!str_detect(name, "harmonic_root05_d"))

n_cores <- parallel::detectCores()

lapply(split(params_df, seq(nrow(params_df))), 
       function(params) wrapper_estim_M(d = as.numeric(params$d), n = as.numeric(params$n),
                                        K = as.numeric(params$K), par = eval(parse(text = params$par)),
                                        name = params$name, dvine = params$dvine, error_const = as.numeric(params$error_const),
                                        seed = params$seed, n_cores = n_cores, save = TRUE, return = FALSE,
                                        h = h_G, score_deriv = score_deriv_G, score_deriv_u = score_deriv_u_G,
                                        h1 = h1_G, h2 = h2_G, h_d = h_d_G))

#### modify alpha_t for D-vine, 0.5/sqrt(t+1) ####
alpha_t <- function(t) t
params_df_d$error_const <- 1e-7

lapply(split(params_df_d, seq(nrow(params_df_d))), 
       function(params) wrapper_estim_M(d = as.numeric(params$d), n = as.numeric(params$n),
                                        K = as.numeric(params$K), par = eval(parse(text = params$par)),
                                        name = params$name, dvine = params$dvine, error_const = as.numeric(params$error_const),
                                        seed = params$seed, n_cores = n_cores, save = TRUE, return = FALSE,
                                        h = h_G, score_deriv = score_deriv_G, score_deriv_u = score_deriv_u_G,
                                        h1 = h1_G, h2 = h2_G, h_d = h_d_G,
                                        alpha_t = alpha_t))

