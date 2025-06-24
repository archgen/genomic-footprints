# ──────────────────────────────────────────────────────────────────────────────
# 1. Build your ggplot (assign to an object 'p' for clarity)
# ──────────────────────────────────────────────────────────────────────────────
setwd("/Users/matthewwilliams/Documents/Research/PSU/Admixture-Migration-Perspectives/genomic-footprints/")
library(ggplot2)
library(viridis)
library(tidyr)

M <- seq(0.01, 1, length.out = 1e6)
theta <- 10000
df <- data.frame(
  M,
  f3_stepping      =  (2    / (7 * M)) * theta/2,
  f4_stepping      = (-8    / (7 * M)) * theta/2,
  f3_hierarchical  = (-0.06 / M) * theta/2,
  f4_hierarchical  = (14    / (55 * M)) * theta/2
)

df_long <- pivot_longer(
  df,
  cols      = -M,
  names_to  = "Statistic",
  values_to = "Value"
)

df_long$Statistic <- factor(
  df_long$Statistic,
  levels = c("f3_stepping", "f4_stepping", "f3_hierarchical", "f4_hierarchical"),
  labels = c(
    "Stepping Stone f3(PX; P1, P2) = (2/7M) * θ/2",
    "Stepping Stone f4(P1, PX; P2, P3) = (-8/7M) * θ/2",
    "Hierarchical Stepping Stone f3(PX; P1, P2) = (-0.06/M) * θ/2",
    "Hierarchical Stepping Stone f4(P1, PX; P2, P3) = (14/55M) * θ/2"
  )
)

p <- ggplot(df_long, aes(x = M, y = Value,
                         color    = Statistic,
                         linetype = Statistic,
                         shape    = Statistic)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_color_viridis_d() +
  scale_linetype_manual(values = c("solid","solid","dashed","dashed")) +
  scale_shape_manual(   values = c(16,17,15,3)) +
  scale_x_continuous(breaks = seq(0,1,0.1)) +
  scale_y_continuous(expand = expansion(mult = c(0.1,0.1))) +
  labs(
    title = "F-Statistics under Stepping Stone Models",
    x     = "Migration Rate (M)",
    y     = "F-Statistic Value"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major    = element_line(color = "gray80", linetype = "dotted"),
    panel.grid.minor    = element_blank(),
    legend.position     = c(0.98, 0.3),
    legend.justification= c("right","top"),
    legend.background   = element_rect(fill = "white", color = NA),
    legend.key          = element_rect(fill = "white", color = NA),
    axis.text.y         = element_text(size=15),
    axis.text.x         = element_text(size=15),
    axis.text           = element_text(size=30),
    legend.text         = element_text(size=15)
  )

# Print to screen (optional)
#print(p)




