
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}

pak::pkg_install(c("batchtools", "here", "ggplot2", "patchwork"))

library(batchtools)
library(here)
library(data.table)
library(ggplot2)
library(patchwork)

alpha <- 0.05

# Load registry with precomputed experiments
reg_dir <- here("registry")
reg <- loadRegistry(reg_dir, writeable = T)
# Save results
dt <- flatten(ijoin(reduceResultsDataTable(), getJobPars()))

# Compute proportion of rejections over the iterations
dt <- dt[, .(num_rjct = sum(result.1 < alpha), k = .N), 
   by = .(problem, n, algorithm)
][, prop_significant := num_rjct/k]

# Wilson 95% confidence interval for binom proportion
dt[, c("ci_low", "ci_high") := {
  ci <- binom::binom.confint(x = num_rjct, n = k, methods = "wilson")
  list(ci$lower, ci$upper)
}]

dt[, problem := factor(problem,   
                         levels = c("null_1", "null_2", "null_3", "null_4",
                                    "alt_1", "alt_2", "alt_3"), 
                         labels = c("Null 1", "Null 2", "Null 3", "Null 4",
                                    "Alt. 1", "Alt. 2", "Alt. 3"))]

dt[, n := factor(n, levels = sort(unique(n)))]

dt[, algorithm := factor(algorithm, 
                         levels = c("GCM","wGCM","PCM",
                                    "KCIT","RCIT","RCoT",
                                    "GKCM_KRR", "GKCM_RF"),
                         labels = c("GCM","wGCM","PCM",
                                    "KCIT","RCIT","RCoT",
                                    "GKCM KRR", "GKCM RF"))]

cols_alg <- c(
  "GCM"   = "#6F4A9F", 
  "wGCM"  = "#51127C", 
  "PCM"   = "#3A0A64",
  "KCIT"  = "#FDB8A0", 
  "RCIT"  = "#FB9471",
  "RCoT"  = "#E77D59",
  "GKCM KRR" = "#CE6E9F",
  "GKCM RF" = "#B63679" 
)

pd <- position_dodge(width = 0.4)

plot_null <- ggplot(
  dt[problem %chin% c("Null 1", "Null 2", "Null 3", "Null 4")],
  aes(x = algorithm, y = prop_significant, group = algorithm, color = algorithm)
) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_errorbar(
    aes(ymin = pmax(ci_low, 0), ymax = pmin(ci_high, 1)),
    position = pd,
    linewidth = 0.6,
    width = 0.6,
    alpha = 1
  ) +
  geom_point(
    position = pd,
    size = 1.5,
    stroke = 0.4,
    alpha = 1
  ) +
  scale_color_manual(values = cols_alg, name = "Algorithm") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.05, 0.05))) +
  labs(
    x = NULL,
    y = NULL,
    title = NULL 
  ) +
  theme_minimal(base_size = 9) +
  facet_grid(
    problem ~ n,
    labeller = labeller(
      n       = function(x) paste0("n = ", x),
      problem = label_value
    ),
    switch = "y" 
  ) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
    strip.placement   = "outside",  
    strip.text.y.left = element_text(face = "bold",
                                     angle = 90, 
                                     size = 11),
    legend.position   = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    panel.spacing.y   = unit(0.5, "lines"), 
    plot.margin       = margin(t = 0, r = 20, b = 0, l = 0), 
    panel.border     = element_rect(color = "grey70", fill = NA, linewidth = 0.4),
    panel.grid.major.y = element_line(color = "grey85", size = 0.3),
    panel.grid.minor.y = element_line(color = "grey85", size = 0.3),
    panel.background = element_blank()
  )

ggsave(filename = here("plot_null.pdf"), plot = plot_null, 
       width = 8, height = 5.5, units = "in", device = "pdf")


plot_alt <- ggplot(
  dt[problem %chin% c("Alt. 1", "Alt. 2", "Alt. 3")],
  aes(x = algorithm, y = prop_significant, group = algorithm, color = algorithm)
) +
  geom_errorbar(
    aes(ymin = pmax(ci_low, 0), ymax = pmin(ci_high, 1)),
    position = pd,
    linewidth = 0.6,
    width = 0.6,
    alpha = 1
  ) +
  geom_point(
    position = pd,
    size = 1.5,
    stroke = 0.4,
    alpha = 1
  ) +
  scale_color_manual(values = cols_alg, name = "Algorithm") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.05, 0.05))) +
  labs(
    x = NULL,
    y = NULL,
    title = NULL 
  ) +
  theme_minimal(base_size = 9) +
  facet_grid(
    problem ~ n,
    labeller = labeller(
      n       = function(x) paste0("n = ", x),
      problem = label_value
    ),
    switch = "y" 
  ) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        strip.placement   = "outside",
        strip.text.y.left = element_text(face = "bold",
                                         angle = 90, 
                                         size = 11),
        legend.position   = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
        panel.spacing.y   = unit(0.5, "lines"),
        plot.margin       = margin(t = 0, r = 20, b = 0, l = 0),
        panel.border     = element_rect(color = "grey70", fill = NA, linewidth = 0.4),
        panel.grid.major.y = element_line(color = "grey85", size = 0.3),
        panel.grid.minor.y = element_line(color = "grey85", size = 0.3),
        panel.background = element_blank()
  )

ggsave(filename = here("plot_alt.pdf"), plot = plot_alt, 
       width = 8, height = 4.3, units = "in", device = "pdf")
