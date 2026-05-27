
path <- "results/" 

source("packages.R")
source("functions/estimate_params.R")


d <- c(10, 25, 50, 75, 100, 150)
family <- "mix" 
par <- c("function (t)  1/(t+1)", "function (t)  0.5^t", "function (t)  0.5/sqrt(t+1)") # par1, par2 = par1
N <- 100
n <- c(100, 500, 1000, 2500, 5000)
dvine <- c(TRUE, FALSE)
params_df <- expand.grid(d, par, N, n, dvine, family)
names(params_df) <- c("d", "par", "N", "n", "dvine", "family")

names_df <- data.frame(par = par,
                       name = c("harmonic",
                                "exponential", "harmonicRoot05"))  
params_df <- left_join(params_df, names_df)
params_df$name <- paste(params_df$name, params_df$family, 
                        ifelse(params_df$dvine, "dvine", "cvine"), params_df$d, params_df$n, sep = "_")

set.seed(123456)
params_df$seed <- sample.int(10000, size = nrow(params_df))  

params_df <- params_df[order(params_df$d, params_df$n), ]

n_cores <- 11

lapply(split(params_df, seq(nrow(params_df))), 
       function(params) estimate_param_mixgumbel(d = as.numeric(params$d), n = as.numeric(params$n),
                                                N = as.numeric(params$N),
                                                par1 = eval(parse(text = params$par)),
                                                name = params$name, seed = as.numeric(params$seed), 
                                                dvine = params$dvine, 
                                                path = path, folder = "estimate_params/",
                                                n_cores = n_cores,
                                                save = TRUE, return = FALSE))
