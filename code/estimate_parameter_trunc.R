
path <- "results/" 

source("packages.R")
source("functions/estimate_params_trunc.R")


d <- c(100, 200, 300, 400, 500, 750, 1000, 1250, 1500, 1750, 2000)
family <- 1
par <- "function (t)  1/(t+1)"
N <- 100
n <- c(100, 1000, 5000)
dvine <- FALSE
params_df <- expand.grid(d, par, N, n, dvine, family)
names(params_df) <- c("d", "par", "N", "n", "dvine", "family")

names_df <- data.frame(par = par,
                       name = c("harmonic"))
params_df <- left_join(params_df, names_df)
params_df$name <- paste(params_df$name, params_df$family, 
                        ifelse(params_df$dvine, "dvine", "cvine"), params_df$d, params_df$n, sep = "_")
set.seed(12345)
params_df$seed <- sample.int(10000, size = nrow(params_df))  

params_df <- params_df[order(params_df$d, params_df$n), ]

n_cores <- parallel::detectCores()


lapply(split(params_df, seq(nrow(params_df))), 
       function(params) estimate_param_trunc(d = as.numeric(params$d), n = as.numeric(params$n),
                                       N = as.numeric(params$N), bicop_fam = as.numeric(params$family),
                                       par = eval(parse(text = params$par)),
                                       name = params$name, seed = as.numeric(params$seed), 
                                       dvine = params$dvine, 
                                       path = path, n_cores = n_cores,
                                       save = TRUE, return = FALSE))
