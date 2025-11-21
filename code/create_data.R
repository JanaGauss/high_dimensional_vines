# create data from raw results

source("packages.R")
source("functions/helpers.R")

#### assumption curvature ####
path <- "results/assumption_curvature"
files <- list.files(path)

res <- data.frame()

for(f in files){
  
  res_f <- readRDS(paste0(path, "/", f))
  df <- res_f$df
  
  df$value <- max(res_f$res)
  
  res <- rbind(res, df)
  
}

write_csv(res, "results/assumption_curvature.csv")

#### assumption M_n ####
path <- "results/estimate_M"
files <- list.files(path) 

res <- data.frame()

for(f in files){
  
  res_f <- readRDS(paste0(path, "/", f))
  df <- res_f$df
  
  df$M_squared <- max(res_f$res$mean_matrix) # M^2
  
  res <- rbind(res, df)
  
}

# estimate D
res$obs_D <- round(1/res$p^2*15*res$n, 0)
ind_dvine <- res$dvine & res$par == "function (t)  1/sqrt(t + 1)"
res$obs_D[ind_dvine] <- round(1/res$p[ind_dvine]*res$n[ind_dvine], 0)
res$D <- NA
for(i in 1:nrow(res)){
  
  res_f <- readRDS(paste0(path, "/", res$name[i], ".rds"))
  obs <- res$obs_D[i]
  
  res$D[i] <- sort(res_f$res$D_vec, decreasing = TRUE)[obs]
  
}

write_csv(res, "results/assumption_M.csv")

#### parameter estimation ####
par_order <- 4:1
gaussian_df <- create_data(substr_files = "_1_c|_1_d",
                           par_order = par_order) %>% filter(n!= 10)

gumbel_df <- create_data(substr_files = "_1004_",
                         par_order = par_order) %>% filter(n!= 10)

student_df <- create_data(substr_files = "_2_",
                          par_order = par_order, par_2 = TRUE) %>% filter(n!= 10)

write_csv(gaussian_df, "results/gaussian_df.csv")
write_csv(gumbel_df, "results/gumbel_df.csv")
write_csv(student_df, "results/student_df.csv")

#### parameter estimation with margins ####
gaussian_df_margins <- create_data(path = "results/estimate_params_margins/",
                                   substr_files = "_1_c|_1_d",
                                   par_order = par_order)

gumbel_df_margins <- create_data(path = "results/estimate_params_margins/",
                                 substr_files = "_1004_",
                                 par_order = par_order)

student_df_margins <- create_data(path = "results/estimate_params_margins/",
                                  substr_files = "_2_",
                                  par_order = par_order, par_2 = TRUE)

write_csv(gaussian_df_margins, "results/gaussian_df_margins.csv")
write_csv(gumbel_df_margins, "results/gumbel_df_margins.csv")
write_csv(student_df_margins, "results/student_df_margins.csv")

#### truncated C-Vine ####
trunc_df <- create_data(path = "results/estimate_params_trunc/")
trunc_df$p <- trunc_df$d*2 - 3
write_csv(trunc_df, "results/trunc_df.csv")

trunc_df_margins <- create_data(path = "results/estimate_params_trunc_margins/")
trunc_df_margins$p <- trunc_df_margins$d*2 - 3
write_csv(trunc_df_margins, "results/trunc_df_margins.csv")
