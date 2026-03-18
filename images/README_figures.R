#%%%%%%%%%%%%%%%%%%%%%%#
#### README Figures ####
#%%%%%%%%%%%%%%%%%%%%%%#

# Load packages
library(ggplot2)

# Set Theta
Theta <- c(seq(-5, -0.10, 0.01), seq(-0.10, 0.10, 0.0001), seq(0.10, 5, 0.01))

# Compute derivatives ----

# Compute derivatives
derivative_df <- data.frame(
  Derivative = rep(
    c(
      "Atan", "EXP", "Log", "Weibull (k = 0.7)", "Weibull (k = 1.7)"
    ), each = length(Theta)
  ),
  Theta = rep(Theta, times = 5),
  Values = c(
    abs(L0ggm:::atan_derivative(x = Theta, lambda = 1)),
    abs(L0ggm:::exp_derivative(x = Theta, lambda = 1)),
    abs(L0ggm:::log_derivative(x = Theta, lambda = 1)),
    abs(L0ggm:::weibull_derivative(x = Theta, lambda = 1, shape = 0.7)),
    abs(L0ggm:::weibull_derivative(x = Theta, lambda = 1, shape = 1.7))
  )
)

# Create factor
derivative_df$Derivative <- factor(derivative_df$Derivative)

# Create plot
derivative_plot <- ggplot(
  data = derivative_df, aes(
    x = Theta, y = Values,
    color = Derivative, group = Derivative,
    linetype = Derivative
  )
) +
  geom_segment(
    aes(x = 0, y = 0, xend = 5), lineend = "square",
    linewidth = 2, show.legend = FALSE, color = "black"
  ) +
  geom_segment(
    aes(x = 0, y = 0, yend = 1), lineend = "square",
    linewidth = 2, show.legend = FALSE, color = "black"
  ) +
  geom_segment(
    aes(x = 0, xend = 5, y = 1),
    lineend = "square", linewidth = 1, color = "black", linetype = "dashed"
  ) +
  geom_line(lineend = "square", linewidth = 1.5) +
  scale_x_continuous(limits = c(0, 5)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(
    name = "Penalty",
    values = c(
      "Atan" = "orange",
      "EXP" = "#53B3CB",
      "Log" = "#AC729F",
      "Weibull (k = 0.7)" = "#C4DF2A",
      "Weibull (k = 1.7)" = "#8BA018"
    ),
    breaks = c("Atan", "EXP", "LASSO", "L0", "Log", "Weibull (k = 0.7)", "Weibull (k = 1.7)")
  ) +
  scale_linetype_manual(
    name = "Penalty",
    values = c(
      "Atan" = "solid",
      "EXP" = "solid",
      "Log" = "solid",
      "Weibull (k = 0.7)" = "solid",
      "Weibull (k = 1.7)" = "solid"
    ),
    breaks = c("Atan", "EXP", "LASSO", "L0", "Log", "Weibull (k = 0.7)", "Weibull (k = 1.7)")
  ) +
  labs(x = "|Theta|") +
  theme_minimal() +
  theme(
    # legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 12)
  )

# Save
ggsave(
  derivative_plot, file = "derivative.png",
  height = 8, width = 10, dpi = 600, bg = "white"
)
