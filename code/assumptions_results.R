
source("packages.R")

#### curvature ####
res <- read.csv("results/assumption_curvature.csv")
res$dvine_label <- factor(ifelse(res$dvine,
                                 "D-Vine",
                                 "C-Vine"), levels = c("D-Vine", "C-Vine"))

res$par_alpha <- paste0(res$par, ", ", res$alpha_t)


res <- res %>% filter(!(dvine & par_alpha == "function (t)  0.5/sqrt(t + 1), function (t)  1") &
                        !(str_sub(name, -1) == "2") & !(alpha_t == "function (t)  t^3") & !(alpha_t == "function (t)  t^2"))

res$par_alpha <- factor(res$par_alpha, levels = unique(res$par_alpha))
label_tex <- c("$\\theta_t = 0$", "$\\theta_t = 0.5^t$", "$\\theta_t =1/(t + 1)$", 
               "$\\theta_t = 0.5/\\sqrt{t + 1 }$", "$\\theta_t = 0.5/\\sqrt{t + 1 }, \\alpha(t) = t$")
names(label_tex) <- unique(res$par_alpha)

plot <- ggplot(data = res, aes(x = d, y = value, col = par_alpha)) + geom_line() + geom_point() +
  facet_wrap(~dvine_label) + geom_hline(yintercept = 0) + theme_bw() + 
  scale_color_brewer(palette = "Set1" , name = "Parameter",  labels = function(l) latex2exp::TeX(label_tex[l])) +
  labs(y = "Estimation A3")
plot 
ggsave("figures/estimation_A3.pdf", plot,
       height = 5, width = 15, units = "cm")


#### M_n and D_n ####
res <- read.csv("results/assumption_M.csv")
res <- res %>% filter(!(par == "function (t)  0"))

res$par <- factor(res$par, levels = c("function (t)  0.5/sqrt(t + 1)", "function (t)  1/(t + 1)", 
                                      "function (t)  0.5^t"))
levels(res$par) <- TeX(c("$\\theta_t = 0.5/\\sqrt{t + 1}$", "$\\theta_t =1/(t + 1)$", 
                         "$\\theta_t = 0.5^t$"))

res_long <- res %>% pivot_longer(cols = c("M_squared", "D"), names_to = "name_2")

res_long$label <- factor(res_long$name_2, levels = c("M_squared", "D"))
levels(res_long$label) <- TeX(c("$M_n^2$", "$D_n$"))

plot1 <- ggplot(data = res_long %>% filter(dvine), aes(x = d^2, y = value)) + 
  facet_wrap(par~label, scales = "free", nrow = 2, labeller = \(x) label_parsed(x, multi_line = FALSE), dir = "v") +
  geom_line() + geom_point() + 
  labs(title = "D-Vine", y = "Estimation", x = TeX("$d^2$")) + 
  theme_bw()
plot1

plot2 <-  ggplot(data = res_long %>% filter(!dvine), aes(x = d^2, y = value)) + geom_line() + geom_point() + 
  facet_wrap(par~label, scales = "free", nrow = 2, labeller = \(x) label_parsed(x, multi_line = FALSE), dir = "v") +
  labs(title = "C-Vine", y = "Estimation", x = TeX("$d^2$")) + 
  theme_bw()
plot2
grid.arrange(plot1, plot2)
ggsave("figures/estimation_M_D.pdf", grid.arrange(plot1, plot2),
       height = 18, width = 18, units = "cm")

