
path <- "results/" 

source("packages.R")
source("functions/estimate_params.R")

d <- c(10, 25, 50, 100, 200)
family <- c(1, 2, 1004) 
par <- c("function (t)  1/(t+1)", "function (t)  0.5^t", "function (t)  0.5/sqrt(t+1)")
N <- 100
n <- c(100, 1000, 5000)
dvine <- c(TRUE, FALSE)
params_df <- expand.grid(d, par, N, n, dvine, family)
names(params_df) <- c("d", "par", "N", "n", "dvine", "family")
params_df$ktau <- FALSE

names_df <- data.frame(par = par,
                       name = c("harmonic",
                                "exponential", "harmonicRoot05")) 
params_df <- left_join(params_df, names_df)
params_df$name <- paste(params_df$name, params_df$family, 
                        ifelse(params_df$dvine, "dvine", "cvine"), params_df$d, params_df$n, sep = "_")
params_df$two_arg <- FALSE

set.seed(123456)
params_df$seed <- sample.int(10000, size = nrow(params_df))  

params_df <- params_df[order(params_df$d, params_df$n), ]

params_df <- params_df %>% filter(!(family == 2 & (d == 200 | (d == 100 & n == 5000)))) # computational reasons
params_df$N[params_df$family == 2] <- 50

n_cores <- parallel::detectCores()

lapply(split(params_df, seq(nrow(params_df))), 
       function(params) estimate_param(d = as.numeric(params$d), n = as.numeric(params$n),
                                       N = as.numeric(params$N), bicop_fam = as.numeric(params$family),
                                       par = eval(parse(text = params$par)),
                                       name = params$name, seed = as.numeric(params$seed), 
                                       two_arg = params$two_arg, dvine = params$dvine, 
                                       path = path, n_cores = n_cores,
                                       ktau = params$ktau,
                                       margins = TRUE,
                                       save = TRUE, return = FALSE))
