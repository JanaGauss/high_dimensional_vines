

plot_sum_error <- function(res_df, title, outlier.size = 0.5, 
                           d_vec = c(10, 25, 50, 100, 200),
                           n_vec = c(100, 1000, 5000),
                           student_nu = FALSE,
                           par_levels = c("function (t)  0.5/sqrt(t + 1)",
                                          "function (t)  1/(t + 1)",
                                          "function (t)  0.5^t"),
                           par_levels_tex = TeX(c("$\\theta_t = 0.5/\\sqrt{t + 1}$", 
                                                  "$\\theta_t =1/(t + 1)$" , 
                                                  "$\\theta_t = 0.5^t$" ))){
  
  
  res_df <- res_df %>% filter(n %in% n_vec & d %in% d_vec & par %in% par_levels)
  
  res_df$par <- factor(res_df$par, levels = par_levels)
  
  levels(res_df$par) <- par_levels_tex
  
  
  plot1 <- ggplot(data = res_df %>% filter(dvine)) + geom_boxplot(aes(x = as.factor(d), y = sum_error*sqrt(n)/d, fill = as.factor(n)), outlier.size = outlier.size) + 
    facet_wrap(~par, ncol = 3, scales = "free_y", labeller = label_parsed) + theme_bw() + 
    labs(title = "D-Vine", fill = "n", x = "d", y = TeX("$\\sqrt{n}/d \\ \\sum \\hat{\\theta}_k - \\theta_k^*$")) +
    geom_hline(yintercept = 0, color = "grey30") +
    scale_fill_brewer(palette = "Set1")
  plot2 <- ggplot(data = res_df %>% filter(!dvine)) + geom_boxplot(aes(x = as.factor(d), y = sum_error*sqrt(n)/d, fill = as.factor(n)), outlier.size = outlier.size) + 
    facet_wrap(~par, ncol = 3, scales = "free_y", labeller = label_parsed) + theme_bw() + 
    labs(title = "C-Vine", fill = "n", x = "d", y = TeX("$\\sqrt{n}/d \\ \\sum \\hat{\\theta}_k - \\theta_k^*$")) +
    geom_hline(yintercept = 0, color = "grey30") +
    scale_fill_brewer(palette = "Set1")
  if(student_nu){
    plot1 <- plot1 + labs(y = TeX("$\\sqrt{n}/d \\ \\sum \\hat{\\nu}_k - \\nu_k^*$"))
    plot2 <- plot2 + labs(y = TeX("$\\sqrt{n}/d \\ \\sum \\hat{\\nu}_k - \\nu_k^*$"))
  }
  grid.arrange(plot1, plot2, top = textGrob(title, gp = gpar(fontsize=20,font=1)))
}

plot_norm_error <- function(res_df, title, outlier.size = 0.5, 
                            d_vec = c(10, 20, 50, 100, 200),
                            n_vec = c(100, 1000, 5000),
                            par_levels = c("function (t)  0.5/sqrt(t + 1)",
                                           "function (t)  1/(t + 1)",
                                           "function (t)  0.5^t"),
                            par_levels_tex = TeX(c("$\\theta_t = 0.5/\\sqrt{t + 1}$", 
                                                   "$\\theta_t =1/(t + 1)$" , 
                                                   "$\\theta_t = 0.5^t$" ))){
  
  
  res_df <- res_df %>% filter(n %in% n_vec & d %in% d_vec & par %in% par_levels)
  
  res_df$par <- factor(res_df$par, levels = par_levels)
  
  levels(res_df$par) <- par_levels_tex
  
  plot1 <- ggplot(data = res_df %>% filter(dvine)) + geom_boxplot(aes(x = as.factor(d), y = norm_error*sqrt(n)/d, fill = as.factor(n)), outlier.size = outlier.size) + 
    facet_wrap(~par, ncol = 3, scales = "free_y", labeller = label_parsed) + theme_bw() + 
    labs(title = "D-Vine", fill = "n", x = "d", y = TeX("$\\frac{\\sqrt{n}}{d} \\| \\hat{\\theta} - \\theta^* \\|$")) +
    geom_hline(yintercept = 0, color = "grey30") +
    scale_fill_brewer(palette = "Set1")
  plot2 <- ggplot(data = res_df %>% filter(!dvine)) + geom_boxplot(aes(x = as.factor(d), y = norm_error*sqrt(n)/d, fill = as.factor(n)), outlier.size = outlier.size) + 
    facet_wrap(~par, ncol = 3, scales = "free_y", labeller = label_parsed) + theme_bw() + 
    labs(title = "C-Vine", fill = "n", x = "d", y = TeX("$\\frac{\\sqrt{n}}{d} \\| \\hat{\\theta} - \\theta^* \\|$")) +
    geom_hline(yintercept = 0, color = "grey30") +
    scale_fill_brewer(palette = "Set1")
  grid.arrange(plot1, plot2, top = textGrob(title, gp = gpar(fontsize=20,font=1)))
}

plot_max_norm_error <- function(res_df, title, outlier.size = 0.5, 
                                d_vec = c(10, 20, 50, 100, 200),
                                n_vec = c(100, 1000, 5000),
                                par_levels = c("function (t)  0.5/sqrt(t + 1)",
                                               "function (t)  1/(t + 1)",
                                               "function (t)  0.5^t"),
                                par_levels_tex = TeX(c("$\\theta_t = 0.5/\\sqrt{t + 1}$", 
                                                       "$\\theta_t =1/(t + 1)$" , 
                                                       "$\\theta_t = 0.5^t$" ))){
  
  
  res_df <- res_df %>% filter(n %in% n_vec & d %in% d_vec & par %in% par_levels)
  
  res_df$par <- factor(res_df$par, levels = par_levels)
  
  levels(res_df$par) <- par_levels_tex
  
  plot1 <- ggplot(data = res_df %>% filter(dvine)) + geom_boxplot(aes(x = as.factor(d), y = max_norm_error*sqrt(n/log(d)), fill = as.factor(n)), outlier.size = outlier.size) + 
    facet_wrap(~par, ncol = 3, scales = "free_y", labeller = label_parsed) + theme_bw() + 
    labs(title = "D-Vine", fill = "n", x = "d", y = TeX("$\\sqrt{n/\\ln(d)} \\| \\hat{\\theta} - \\theta^* \\|_{\\infty}$")) +
    geom_hline(yintercept = 0, color = "grey30") +
    scale_fill_brewer(palette = "Set1")
  plot2 <- ggplot(data = res_df %>% filter(!dvine)) + geom_boxplot(aes(x = as.factor(d), y = max_norm_error*sqrt(n/log(d)), fill = as.factor(n)), outlier.size = outlier.size) + 
    facet_wrap(~par, ncol = 3, scales = "free_y", labeller = label_parsed) + theme_bw() + 
    labs(title = "C-Vine", fill = "n", x = "d", y = TeX("$\\sqrt{n/\\ln(d)} \\| \\hat{\\theta} - \\theta^* \\|_{\\infty}$")) +
    geom_hline(yintercept = 0, color = "grey30") +
    scale_fill_brewer(palette = "Set1")
  grid.arrange(plot1, plot2, top = textGrob(title, gp = gpar(fontsize=20,font=1)))
}

plot_last_param_error <- function(res_df, title, outlier.size = 0.5, 
                                  d_vec = c(10, 20, 50, 100, 200),
                                  n_vec = c(100, 1000, 5000),
                                  par_levels = c("function (t)  0.5/sqrt(t + 1)",
                                                 "function (t)  1/(t + 1)",
                                                 "function (t)  0.5^t"),
                                  par_levels_tex = TeX(c("$\\theta_t = 0.5/\\sqrt{t + 1}$", 
                                                         "$\\theta_t =1/(t + 1)$" , 
                                                         "$\\theta_t = 0.5^t$" ))){
  
  
  res_df <- res_df %>% filter(n %in% n_vec & d %in% d_vec & par %in% par_levels)
  
  res_df$par <- factor(res_df$par, levels = par_levels)
  
  levels(res_df$par) <- par_levels_tex
  
  
  plot1 <- ggplot(data = res_df %>% filter(dvine)) + geom_boxplot(aes(x = as.factor(d), y = last_param_error*sqrt(n), fill = as.factor(n)), outlier.size = outlier.size) + 
    facet_wrap(~par, ncol = 3, scales = "free_y", labeller = label_parsed) + theme_bw() + 
    labs(title = "D-Vine", fill = "n", x = "d", y = TeX("$\\sqrt{n} ( \\hat{\\theta}_p - \\theta_p^*)$")) +
    geom_hline(yintercept = 0, color = "grey30") +
    scale_fill_brewer(palette = "Set1")
  plot2 <- ggplot(data = res_df %>% filter(!dvine)) + geom_boxplot(aes(x = as.factor(d), y = last_param_error*sqrt(n), fill = as.factor(n)), outlier.size = outlier.size) + 
    facet_wrap(~par, ncol = 3, scales = "free_y", labeller = label_parsed) + theme_bw() + 
    labs(title = "C-Vine", fill = "n", x = "d", y = TeX("$\\sqrt{n} ( \\hat{\\theta}_p - \\theta_p^*)$")) +
    geom_hline(yintercept = 0, color = "grey30") +
    scale_fill_brewer(palette = "Set1")
  
  grid.arrange(plot1, plot2, top = textGrob(title, gp = gpar(fontsize=20,font=1)))
}

plot_last_param_error_norm <- function(res_df, title, outlier.size = 0.5, 
                                       d_vec = c(10, 20, 50, 100, 200),
                                       n_vec = c(100, 1000, 5000),
                                       par_levels = c("function (t)  0.5/sqrt(t + 1)",
                                                      "function (t)  1/(t + 1)",
                                                      "function (t)  0.5^t"),
                                       par_levels_tex = TeX(c("$\\theta_t = 0.5/\\sqrt{t + 1}$", 
                                                              "$\\theta_t =1/(t + 1)$" , 
                                                              "$\\theta_t = 0.5^t$" ))){
  
  
  res_df <- res_df %>% filter(n %in% n_vec & d %in% d_vec & par %in% par_levels)
  
  res_df$par <- factor(res_df$par, levels = par_levels)
  
  levels(res_df$par) <- par_levels_tex
  
  
  plot1 <- ggplot(data = res_df %>% filter(dvine)) + geom_boxplot(aes(x = as.factor(d), y = abs(last_param_error)*sqrt(n), fill = as.factor(n)), outlier.size = outlier.size) + 
    facet_wrap(~par, ncol = 3, scales = "free_y", labeller = label_parsed) + theme_bw() + 
    labs(title = "D-Vine", fill = "n", x = "d", y = TeX("$\\sqrt{n} | \\hat{\\theta}_p - \\theta_p^*|$")) +
    geom_hline(yintercept = 0, color = "grey30") +
    scale_fill_brewer(palette = "Set1")
  plot2 <- ggplot(data = res_df %>% filter(!dvine)) + geom_boxplot(aes(x = as.factor(d), y = abs(last_param_error)*sqrt(n), fill = as.factor(n)), outlier.size = outlier.size) + 
    facet_wrap(~par, ncol = 3, scales = "free_y", labeller = label_parsed) + theme_bw() + 
    labs(title = "C-Vine", fill = "n", x = "d", y = TeX("$\\sqrt{n} | \\hat{\\theta}_p - \\theta_p^*|$")) +
    geom_hline(yintercept = 0, color = "grey30") +
    scale_fill_brewer(palette = "Set1")
  
  grid.arrange(plot1, plot2, top = textGrob(title, gp = gpar(fontsize=20,font=1)))
}

plot_sum_error_margins <- function(res_df, title, outlier.size = 0.5, 
                                   d_vec = c(10, 25, 50, 100, 200),
                                   n_vec = c(100, 1000, 5000),
                                   student_nu = FALSE,
                                   par_levels = c("function (t)  0.5/sqrt(t + 1)",
                                                  "function (t)  1/(t + 1)",
                                                  "function (t)  0.5^t"),
                                   par_levels_tex = TeX(c("$\\theta_t = 0.5/\\sqrt{t + 1}$", 
                                                          "$\\theta_t =1/(t + 1)$" , 
                                                          "$\\theta_t = 0.5^t$" ))){
  
  
  res_df <- res_df %>% filter(n %in% n_vec & d %in% d_vec & par %in% par_levels)
  
  res_df$par <- factor(res_df$par, levels = par_levels)
  
  levels(res_df$par) <- par_levels_tex
  
  plot1 <- ggplot(data = res_df %>% filter(dvine)) + geom_boxplot(aes(x = as.factor(d), y = sum_error*sqrt(n)/d, color = as.factor(n), fill = margins), 
                                                                  outlier.size = outlier.size, lwd = 0.7) + 
    facet_wrap(~par, scales = "free_y", labeller = label_parsed) + theme_bw() + 
    labs(title = "D-Vine", fill = "estim. of \nmargins", col = "n", x = "d", y = TeX("$\\sqrt{n}/d\\ \\sum \\hat{\\theta}_k - \\theta_k^*$")) +
    geom_hline(yintercept = 0, color = "grey30") +
    scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "grey70"), labels = c("FALSE" = "no", "TRUE" = "yes")) + 
    scale_color_brewer(palette = "Set1")
  plot2  <- ggplot(data = res_df %>% filter(!dvine)) + geom_boxplot(aes(x = as.factor(d), y = sum_error*sqrt(n)/d, color = as.factor(n), fill = margins), 
                                                                    outlier.size = outlier.size, lwd = 0.7) + 
    facet_wrap(~par, scales = "free_y", labeller = label_parsed) + theme_bw() + 
    labs(title = "C-Vine", fill = "estim. of \nmargins", col = "n", x = "d", y = TeX("$\\sqrt{n}/d \\ \\sum \\hat{\\theta}_k - \\theta_k^*$")) +
    geom_hline(yintercept = 0, color = "grey30") +
    scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "grey70"), labels = c("FALSE" = "no", "TRUE" = "yes")) + 
    scale_color_brewer(palette = "Set1")
  if(student_nu){
    plot1 <- plot1 + labs(y = TeX("$\\frac{\\sqrt{n}}{d}\\sum \\hat{\\nu}_k - \\nu_k^*$"))
    plot2 <- plot2 + labs(y = TeX("$\\frac{\\sqrt{n}}{d}\\sum \\hat{\\nu}_k - \\nu_k^*$"))
  }
  grid.arrange(plot1, plot2, top = textGrob(title, gp = gpar(fontsize=20,font=1)))
  
}

plot_norm_error_margins <- function(res_df, title, outlier.size = 0.5, 
                                    d_vec = c(10, 20, 50, 100, 200),
                                    n_vec = c(100, 1000, 5000),
                                    par_levels = c("function (t)  0.5/sqrt(t + 1)",
                                                   "function (t)  1/(t + 1)",
                                                   "function (t)  0.5^t"),
                                    par_levels_tex = TeX(c("$\\theta_t = 0.5/\\sqrt{t + 1}$", 
                                                           "$\\theta_t =1/(t + 1)$" , 
                                                           "$\\theta_t = 0.5^t$" ))){
  
  
  res_df <- res_df %>% filter(n %in% n_vec & d %in% d_vec & par %in% par_levels)
  
  res_df$par <- factor(res_df$par, levels = par_levels)
  
  levels(res_df$par) <- par_levels_tex
  
  plot1 <- ggplot(data = res_df %>% filter(dvine)) + geom_boxplot(aes(x = as.factor(d), y = norm_error*sqrt(n)/d, color = as.factor(n), fill = margins), 
                                                                  outlier.size = outlier.size, lwd = 0.7) + 
    facet_wrap(~par, scales = "free_y", labeller = label_parsed) + theme_bw() + 
    labs(title = "D-Vine", fill = "estimation of margins", x = "d", y = TeX("$\\frac{\\sqrt{n}}{d} \\| \\hat{\\theta} - \\theta^* \\|$")) +
    geom_hline(yintercept = 0, color = "grey30") +
    scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "grey70")) + 
    scale_color_brewer(palette = "Set1")
  plot2  <- ggplot(data = res_df %>% filter(!dvine)) + geom_boxplot(aes(x = as.factor(d), y = norm_error*sqrt(n)/d, color = as.factor(n), fill = margins), 
                                                                    outlier.size = outlier.size, lwd = 0.7) + 
    facet_wrap(~par, scales = "free_y", labeller = label_parsed) + theme_bw() + 
    labs(title = "C-Vine", fill = "estimation of margins", x = "d", y = TeX("$\\frac{\\sqrt{n}}{d} \\| \\hat{\\theta} - \\theta^* \\|$")) +
    geom_hline(yintercept = 0, color = "grey30") +
    scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "grey70")) + 
    scale_color_brewer(palette = "Set1")
  
  grid.arrange(plot1, plot2, top = textGrob(title, gp = gpar(fontsize=20,font=1)))
  
}


plot_max_norm_error_margins <- function(res_df, title, outlier.size = 0.5, 
                                        d_vec = c(10, 20, 50, 100, 200),
                                        n_vec = c(100, 1000, 5000),
                                        par_levels = c("function (t)  0.5/sqrt(t + 1)",
                                                       "function (t)  1/(t + 1)",
                                                       "function (t)  0.5^t"),
                                        par_levels_tex = TeX(c("$\\theta_t = 0.5/\\sqrt{t + 1}$", 
                                                               "$\\theta_t =1/(t + 1)$" , 
                                                               "$\\theta_t = 0.5^t$" ))){
  
  
  res_df <- res_df %>% filter(n %in% n_vec & d %in% d_vec & par %in% par_levels)
  
  res_df$par <- factor(res_df$par, levels = par_levels)
  
  levels(res_df$par) <- par_levels_tex
  
  plot1 <- ggplot(data = res_df %>% filter(dvine)) + geom_boxplot(aes(x = as.factor(d), y = max_norm_error*sqrt(n/log(d)), color = as.factor(n), fill = margins), 
                                                                  outlier.size = outlier.size, lwd = 0.7) + 
    facet_wrap(~par, scales = "free_y", labeller = label_parsed) + theme_bw() + 
    labs(title = "D-Vine", fill = "estimation of margins", x = "d", y = TeX("$\\sqrt{n/\\ln(d)} \\| \\hat{\\theta} - \\theta^* \\|_{\\infty}$")) +
    geom_hline(yintercept = 0, color = "grey30") +
    scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "grey70")) + 
    scale_color_brewer(palette = "Set1")
  plot2  <- ggplot(data = res_df %>% filter(!dvine)) + geom_boxplot(aes(x = as.factor(d), y = max_norm_error*sqrt(n/log(d)), color = as.factor(n), fill = margins), 
                                                                    outlier.size = outlier.size, lwd = 0.7) + 
    facet_wrap(~par, scales = "free_y", labeller = label_parsed) + theme_bw() + 
    labs(title = "C-Vine", fill = "estimation of margins", x = "d", y = TeX("$\\sqrt{n/\\ln(d)} \\| \\hat{\\theta} - \\theta^* \\|_{\\infty}$")) +
    geom_hline(yintercept = 0, color = "grey30") +
    scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "grey70")) + 
    scale_color_brewer(palette = "Set1")
  
  grid.arrange(plot1, plot2, top = textGrob(title, gp = gpar(fontsize=20,font=1)))
  
}


plot_consistency <- function(res_df, title, 
                             student_nu = FALSE,
                             par_levels = c("function (t)  0.5/sqrt(t + 1)",
                                          "function (t)  1/(t + 1)",
                                          "function (t)  0.5^t"),
                             par_levels_tex = TeX(c("$\\theta_t = 0.5/\\sqrt{t + 1}$", 
                                                    "$\\theta_t =1/(t + 1)$" , 
                                                    "$\\theta_t = 0.5^t$" ))){
  
  
  res_df$par <- factor(res_df$par, levels = par_levels)
  
  levels(res_df$par) <- par_levels_tex
  
  
  plot1 <- ggplot(data = res_df %>% filter(dvine), aes(d, y = error_n_new*sqrt(n_new/log(d)), col = factor(f_n))) + geom_line() + geom_point() +
    facet_wrap(~par, ncol = 3, labeller = label_parsed) + theme_bw() + 
    labs(title = "D-Vine", col =  "d vs n", x = "d", y = TeX("$\\sqrt{n/\\ln(d)} \\  \\| \\hat{\\theta} - \\theta^* \\|_{\\infty}$")) +
    scale_color_brewer(palette = "Set1", 
                       labels = unname(TeX(c("$n \\sim \\, d$", "$n \\sim \\, d^2$", "$n \\sim \\, d^3$"))))
  plot2 <- ggplot(data = res_df %>% filter(!dvine), aes(d, y = error_n_new*sqrt(n_new/log(d)), col = factor(f_n))) + geom_line() + geom_point() +
    facet_wrap(~par, ncol = 3, labeller = label_parsed) + theme_bw() + 
    labs(title = "C-Vine", col =  "d vs n", x = "d", y = TeX("$\\sqrt{n/\\ln(d)} \\  \\| \\hat{\\theta} - \\theta^* \\|_{\\infty}$")) +
    scale_color_brewer(palette = "Set1",
                       labels = unname(TeX(c("$n \\sim \\, d$", "$n \\sim \\, d^2$", "$n \\sim \\, d^3$"))))
  if(student_nu){
    plot1 <- plot1 + labs(y = TeX("$\\sqrt{n/\\ln(d)} \\  \\| \\hat{\\nu} - \\nu^* \\|_{\\infty}$"))
    plot2 <- plot2 + labs(y = TeX("$\\sqrt{n/\\ln(d)} \\  \\| \\hat{\\nu} - \\nu^* \\|_{\\infty}$"))
  }
  grid.arrange(plot1, plot2, top = textGrob(title, gp = gpar(fontsize=20,font=1)))
}
