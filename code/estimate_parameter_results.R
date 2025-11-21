
source("packages.R")
source("functions/helpers.R")
source("functions/plot_functions.R")


#### parameter estimation without margins #### 
gaussian_df <- read_csv("results/gaussian_df.csv") 
gumbel_df <- read_csv("results/gumbel_df.csv") 
student_df <- read_csv("results/student_df.csv") 

# maximum norm of error, Gauss and Gumbel
p_gauss_norm <- plot_max_norm_error(gaussian_df, title = "Gaussian")
p_gumbel_norm <- plot_max_norm_error(gumbel_df, title = TeX("Gumbel"))

# maximum norm of error, student 
plot_max_norm_error(student_df, title = TeX("Student's t ($\\rho$)"))
plot_max_norm_error(student_df %>% mutate(max_norm_error = max_norm_error2), 
                    title = TeX("Student's t (degrees of freedom $\\nu$, $\\nu^* = 4$)"))

# sum of error, Gauss and Gumbel
p_gauss <- plot_sum_error(gaussian_df, title = "Gaussian")
p_gumbel <- plot_sum_error(gumbel_df, 
                           title = "Gumbel")

ggsave("figures/par_estim_sum_gauss_gumbel.pdf", grid.arrange(p_gauss, p_gumbel, ncol = 1),
       height = 30, width = 25, units = "cm")



p_student1 <- plot_sum_error(student_df, title = TeX("Student's t ($\\rho$)"))
p_student2 <- plot_sum_error(student_df %>% mutate(sum_error = sum_error2), 
                             title = TeX("Student's t (degrees of freedom $\\nu$, $\\nu^* = 4$)"),
                             student_nu = TRUE) 

ggsave("figures/par_estim_sum_student.pdf", p_student1,
       height = 15, width = 25, units = "cm")


#### interpolation to verify consistency ####

dat_interpol_gaussian <- create_dat_interpol(gaussian_df)
dat_interpol_gumbel <- create_dat_interpol(gumbel_df)
dat_interpol_student <- create_dat_interpol(student_df)
dat_interpol_student_nu <- create_dat_interpol(student_df %>% 
                                                 mutate(max_norm_error = max_norm_error2))


plot_cons_gauss <- plot_consistency(dat_interpol_gaussian, title = "Gaussian")
plot_cons_gumbel <- plot_consistency(dat_interpol_gumbel, title = "Gumbel")
plot_cons_student1 <- plot_consistency(dat_interpol_student, title = TeX("Student's t ($\\rho$)"))
plot_cons_student2 <- plot_consistency(dat_interpol_student_nu %>% filter(d != 10), 
                 title = TeX("Student's t (degrees of freedom $\\nu$, $\\nu^* = 4$)"),
                 student_nu = TRUE)

ggsave("figures/par_estim_cons_gauss_gumbel.pdf", grid.arrange(plot_cons_gauss, plot_cons_gumbel, ncol = 1),
       height = 23, width = 30, units = "cm")

ggsave("figures/par_estim_cons_student.pdf", grid.arrange(plot_cons_student1, plot_cons_student2, ncol = 1),
       height = 23, width = 30, units = "cm")


### D-Vine with slower rate of convergence ####
dat_dvine_root <- dat_interpol_gaussian %>% 
  filter(dvine & par == "function (t)  0.5/sqrt(t + 1)" & n_new > 0) %>%
  mutate(error_1 = error_n_new*sqrt(n_new/log(d)),
         error_2 = error_n_new*sqrt(n_new/log(d)/d),
         error_3 = error_n_new*sqrt(n_new/log(d)/d^2))


dat_dvine_root_long <- dat_dvine_root %>% pivot_longer(cols = paste0("error_", 1:3))

dat_dvine_root_long$label <- factor(dat_dvine_root_long$name)
levels(dat_dvine_root_long$label) <- TeX(c("$r_n = \\sqrt{\\log ( d) / n}$", 
                                           "$r_n = \\sqrt{d \\cdot \\log ( d) / n}$", 
                                           "$r_n = \\sqrt{d^2 \\cdot \\log ( d) / n}$"))

plot_dvine <- ggplot(dat_dvine_root_long, aes(x = d, y = value, col = factor(f_n))) +
  geom_line() + geom_point() + 
  facet_wrap(~label, scales = "free", labeller = label_parsed) +
  labs(x = "d", col = "d vs n", y = TeX("$ r_n^{-1} \\ \\| \\hat{\\theta} - \\theta^* \\|_{\\infty}$")) + 
  theme_bw() +
  expand_limits(y = 0) + 
  scale_color_brewer(palette = "Set1",
                     labels = unname(TeX(c("$n \\sim \\, d$", "$n \\sim \\, d^2$", "$n \\sim \\, d^3$"))))
plot_dvine

ggsave("figures/par_estim_dvine.pdf", plot_dvine,
       height = 5, width = 25, units = "cm")

#### estimation of margins ####

gaussian_df$margins <- FALSE
gumbel_df$margins <- FALSE
student_df$margins <- FALSE

gaussian_df_margins <- read_csv("results/gaussian_df_margins.csv")
gaussian_df_margins$margins <- TRUE
gaussian_df <- rbind(gaussian_df, gaussian_df_margins)


gumbel_df_margins <- read_csv("results/gumbel_df_margins.csv")
gumbel_df_margins$margins <- TRUE
gumbel_df <- rbind(gumbel_df, gumbel_df_margins)

student_df_margins <- read_csv("results/student_df_margins.csv")
student_df_margins$margins <- TRUE
student_df <- rbind(student_df, student_df_margins)



p_gauss_margins <- plot_sum_error_margins(gaussian_df %>% filter(d > 10), title = "Gaussian")
ggsave("figures/par_estim_margins_gauss.pdf", p_gauss_margins,
       height = 15, width = 25, units = "cm")

# not shown in paper: Gumbel, student
p_gumbel_margins <- plot_sum_error_margins(gumbel_df, title = "Gumbel")
p_student_margins1 <- plot_sum_error_margins(student_df %>% filter(d <= 50), title = TeX("Student's t ($\\rho$)"))
p_student_margins2 <- plot_sum_error_margins(student_df %>% filter(d <= 50) %>% mutate(sum_error = sum_error2), 
                                             title = TeX("Student's t (degrees of freedom $\\nu$, $\\nu^* = 4$)"),
                                             student_nu = TRUE)

plot_max_norm_error_margins(gaussian_df, title = "Gaussian")
plot_max_norm_error_margins(gumbel_df, title = TeX("Gumbel"))
plot_max_norm_error_margins(student_df, title = TeX("Student's t ($\\rho$)"))
plot_max_norm_error_margins(student_df %>% mutate(sum_error = sum_error2), title = TeX("Student's t (degrees of freedom $\\nu$, $\\nu^* = 4$)"))



#### truncated C Vine ####

trunc_df <- read_csv("results/trunc_df.csv")

trunc_df <- trunc_df %>% filter(!(d %in% c(100, 300, 400)))

p_max_norm_trunc <- ggplot(data = trunc_df) + geom_boxplot(aes(x = as.factor(d), y = max_norm_error*(sqrt(n/log(d))), fill = as.factor(n)), outlier.size = 0.5) + 
  theme_bw() + scale_fill_brewer(palette = "Set1") +
  labs(x = "d", fill = "n",  y = TeX("$\\sqrt{n/ \\ln(d)} \\ \\| \\hat{\\theta} - \\theta^* \\|_{\\infty}$")) +
  guides(fill = "none")
p_max_norm_trunc

p_sum_trunc <- ggplot(data = trunc_df) + geom_boxplot(aes(x = as.factor(d), y = sum_error*(sqrt(n))/d, fill = as.factor(n)), outlier.size = 0.5) +
  geom_hline(yintercept = 0, col = "grey30") + theme_bw() +
  scale_fill_brewer(palette = "Set1") + labs(x = "d", fill = "n", y = TeX("$\\sqrt{n}/d\\ \\sum \\hat{\\theta}_k - \\theta_k^*$")) 
p_sum_trunc

ggsave("figures/par_estim_trunc.pdf", 
       p_sum_trunc, 
       height = 7, width = 12, units = "cm")


# with estimation of margins, not shown in paper
trunc_df$margins <- FALSE
trunc_df_marg <- read_csv("results/trunc_df_margins.csv")  %>% filter(!(d %in% c(100, 300, 400)))
trunc_df_marg$margins <- TRUE

trunc_df <- rbind(trunc_df, trunc_df_marg)

ggplot(data = trunc_df) + geom_boxplot(aes(x = as.factor(d), y = norm_error*(sqrt(n))/sqrt(p), col = as.factor(n), fill = margins), outlier.size = 0.5, lwd = 0.8) + 
  theme_bw() + scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "grey70"))
ggplot(data = trunc_df) + geom_boxplot(aes(x = as.factor(d), y = sum_error*(sqrt(n))/p, col = as.factor(n), fill = margins), outlier.size = 0.5, lwd = 0.8) +
  geom_hline(yintercept = 0, col = "grey30") + theme_bw() + 
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "grey70"))


