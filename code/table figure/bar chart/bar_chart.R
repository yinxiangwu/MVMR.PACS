# Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

# ----------------------------
# 1) Enter the table as a data frame
# ----------------------------
dat <- tribble(
  ~Estimator,        ~N,     ~MSE,   ~CorrectSparsity, ~Sensitivity, ~FalsePositive,
  "MVMR-IVW",       "1e5",  1.038,  0.656,            0.222,        0.055,
  "MVMR-IVW",       "2e5",  0.711,  0.729,            0.376,        0.035,
  "MVMR-IVW",       "3e5",  0.625,  0.767,            0.458,        0.027,
  
  "SRIVW",          "1e5",  5.874,  0.735,            0.341,        0.003,
  "SRIVW",          "2e5",  8.972,  0.712,            0.283,        0.002,
  "SRIVW",          "3e5",  8.051,  0.703,            0.262,        0.003,
  
  "IVW-LASSO",      "1e5",  1.127,  0.848,            0.701,        0.054,
  "IVW-LASSO",      "2e5",  0.941,  0.899,            0.789,        0.028,
  "IVW-LASSO",      "3e5",  0.834,  0.943,            0.896,        0.025,
  
  "MVMR-dLASSO",    "1e5",  5.717,  0.732,            0.331,        0.000,
  "MVMR-dLASSO",    "2e5",  5.858,  0.770,            0.427,        0.000,
  "MVMR-dLASSO",    "3e5",  4.305,  0.818,            0.545,        0.001,
  
  "MVMR-PACS",      "1e5",  0.279,  0.887,            0.834,        0.077,
  "MVMR-PACS",      "2e5",  0.225,  0.945,            0.887,        0.017,
  "MVMR-PACS",      "3e5",  0.100,  0.973,            0.941,        0.005,
  
  "MVMR-PACS-0.8",  "1e5",  0.254,  0.912,            0.784,        0.002,
  "MVMR-PACS-0.8",  "2e5",  0.250,  0.944,            0.866,        0.004,
  "MVMR-PACS-0.8",  "3e5",  0.069,  0.976,            0.942,        0.001,
  
  "MRBMA-0.2",      "1e5",  0.559,  0.415,            0.984,        0.964,
  "MRBMA-0.2",      "2e5",  0.345,  0.499,            0.999,        0.834,
  "MRBMA-0.2",      "3e5",  0.284,  0.581,            1.000,        0.699,
  
  "MRBMA-0.5",      "1e5",  0.559,  0.752,            0.859,        0.318,
  "MRBMA-0.5",      "2e5",  0.345,  0.853,            0.972,        0.226,
  "MRBMA-0.5",      "3e5",  0.284,  0.892,            0.994,        0.176,
  
  "MRBMA-0.8",      "1e5",  0.559,  0.816,            0.697,        0.104,
  "MRBMA-0.8",      "2e5",  0.345,  0.902,            0.855,        0.066,
  "MRBMA-0.8",      "3e5",  0.284,  0.936,            0.921,        0.053,
  
  "MVMR-cML-SuSiE",          "1e5",  NA,     0.527,            0.397,        0.386,
  "MVMR-cML-SuSiE",          "2e5",  NA,     0.663,            0.493,        0.223,
  "MVMR-cML-SuSiE",          "3e5",  NA,     0.745,            0.565,        0.136
)


# Order methods for plotting
method_order <- c("MVMR-PACS","MVMR-PACS-0.8","MVMR-dLASSO",
                  "IVW-LASSO","MVMR-IVW","SRIVW",
                  "MRBMA-0.2","MRBMA-0.5","MRBMA-0.8",
                  "MVMR-cML-SuSiE")
dat <- dat %>%
  mutate(Estimator = factor(Estimator, levels = method_order),
         N = factor(N, levels = c("1e5","2e5","3e5")))

# ----------------------------
# 2) Plot helper
# ----------------------------
# Theme
# --- assume 'dat' is already defined ---

# Nature-style discrete color palette (muted and colorblind-safe)
method_colors <- c(
  "MVMR-PACS"       = "#E64B35", # warm red
  "MVMR-PACS-0.8"   = "#F39B7F", # orange (related to PACS)
  "MVMR-dLASSO"     = "#7E6148", # brown
  "IVW-LASSO"       = "#00A087", # teal/green
  "MVMR-IVW"        = "#AE3D63", # deep magenta-violet
  "SRIVW"           = "#8491B4", # purple-gray (distinct from beige)
  "MRBMA-0.2"       = "#91D1C2", # mint-blue
  "MRBMA-0.5"       = "#4DBBD5", # sky blue
  "MRBMA-0.8"       = "#3C5488", # navy blue
  "MVMR-cML-SuSiE"  = "#B09C85"  # beige/neutral
)


# Consistent theme for Nature-style plots
bar_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
  )

# ----------------------------
# 1) MSE Panel (baseline = log(0.01))
# ----------------------------
mse_floor <- 0.01

p_mse <- ggplot(dat %>% mutate(MSE_clamped = pmax(MSE, mse_floor)),
                aes(x = N, fill = Estimator)) +
  geom_rect(aes(
    xmin = as.numeric(factor(N)) - 0.35 + as.numeric(factor(Estimator)) * 0.07,
    xmax = as.numeric(factor(N)) - 0.35 + as.numeric(factor(Estimator)) * 0.07 + 0.06,
    ymin = log10(mse_floor),
    ymax = log10(MSE_clamped)
  ),
  colour = "grey20", linewidth = 0.2) +
  scale_y_continuous(
    breaks = log10(c(0.1, 1, 10)),
    labels = c( "0.1", "1", "10"),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = method_colors) +
  labs(x = "Sample size (N)", y = "MSE (log10 scale)") +
  bar_theme

# ----------------------------
# 2) Other Panels
# ----------------------------
make_bar <- function(df, yvar, ylab) {
  ggplot(df, aes(x = N, y = .data[[yvar]], fill = Estimator)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7,
             colour = "grey20", linewidth = 0.2) +
    scale_fill_manual(values = method_colors) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Sample size (N)", y = ylab) +
    bar_theme
}

p_cs   <- make_bar(dat, "CorrectSparsity", "Correct sparsity")
p_sens <- make_bar(dat, "Sensitivity", "Sensitivity")
p_fn   <- make_bar(dat %>% mutate(FalseNegative = 1 - Sensitivity),
                   "FalsePositive", "False Positives")

# ----------------------------
# 3) Combine all panels
# ----------------------------
combined_plot <- 
  (p_mse | p_cs) / (p_sens | p_fn) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combined_plot

jpeg("~/OneDrive - UW/UW biostats/2023-24 RA/MR with highly correlated exposure/manuscript/code/table figure/bar chart/sim_figure1.jpeg",width = 6400, height = 5600, res = 600)
combined_plot
dev.off()