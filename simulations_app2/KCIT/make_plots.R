
# Install necessary packages
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}

pak::pkg_install(c("data.table", "ggplot2", "here"))

library(data.table)
library(ggplot2)
library(here)

load(here("simulations_app2/KCIT/results_df.Rdata"))

dt <- as.data.table(res_all)

# Summarise rejection rates for alpha = 0.01, 0.05
alphas <- c(0.01, 0.05)

sum_dt <- rbindlist(lapply(alphas, function(alpha) {
  tmp <- copy(dt)
  tmp[, reject := (p_value < alpha)]
  out <- tmp[, .(num = sum(reject), k = .N), by = .(case, CI, n, d, alg_name)]
  out[, `:=`(alpha = alpha,
             prop_reject = num/k,
             scenario = ifelse(CI, "Null", "Alt."))]
  out
}))

sum_dt[, alg_name := as.character(alg_name)]
sum_dt[, algorithm := alg_name]
sum_dt[alg_name == "GKCM_KRR", algorithm := "GKCM KRR"]
sum_dt[alg_name == "GKCM_RF",  algorithm := "GKCM RF"]
sum_dt[, algorithm := factor(algorithm,
                             levels = c("GCM","wGCM","PCM","KCIT","RCIT","RCoT","GKCM KRR","GKCM RF"))]
sum_dt[, n := factor(n)]
sum_dt[, alpha_f := factor(alpha, levels = c(0.01,0.05), labels = c("0.01","0.05"))]
sum_dt[, scenario := factor(scenario, levels = c("Null", "Alt."))]


plot_caseI <- ggplot(sum_dt[case == "I"],
                   aes(x = d, y = algorithm, fill = prop_reject)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", prop_reject)), size = 2.6) +
  scale_x_continuous(breaks = 1:5) +
  scale_fill_gradientn(
    colours = c("green3", "yellow", "orange", "red3"),
    values  = scales::rescale(c(0.00, 0.15, 0.5, 1.00)),
    limits  = c(0, 1),
    name    = "Rejection rate"
  ) + 
  facet_grid(scenario ~ n + alpha_f,
             labeller = labeller(n = function(x) paste0("n = ", x),
                                 alpha_f = function(x) paste0("alpha = ", x))) +
  labs(title = paste0("Case I"),
       x = "No. of covariates",
       y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
        strip.placement = "outside",
        legend.position = "right")

ggsave(filename = here("simulations_app2/KCIT/plot_case_I.pdf"), plot = plot_caseI, 
       width = 8, height = 5.5, units = "in", device = "pdf")


plot_caseII <- ggplot(sum_dt[case == "II"],
                   aes(x = d, y = algorithm, fill = prop_reject)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", prop_reject)), size = 2.6) +
  scale_x_continuous(breaks = 1:5) +
  scale_fill_gradientn(
    colours = c("green3", "yellow", "orange", "red3"),
    values  = scales::rescale(c(0.00, 0.15, 0.5, 1.00)),
    limits  = c(0, 1),
    name    = "Rejection rate"
  ) + 
  facet_grid(scenario ~ n + alpha_f,
             labeller = labeller(n = function(x) paste0("n = ", x),
                                 alpha_f = function(x) paste0("alpha = ", x))) +
  labs(title = paste0("Case II"),
       x = "No. of covariates",
       y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
        strip.placement = "outside",
        legend.position = "right")

ggsave(filename = here("simulations_app2/KCIT/plot_case_II.pdf"), plot = plot_caseII, 
       width = 8, height = 5.5, units = "in", device = "pdf")

