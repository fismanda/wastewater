# =============================================================================
# Figure 2: Linear and log-log regression scatter plots
# Fisman et al. — Wastewater Validation Paper
# =============================================================================
# Four-panel figure (2x2):
#   A: Linear scale — Reported cases ~ WWS
#   B: Linear scale — Adjusted cases ~ WWS
#   C: Log-log scale — Reported cases ~ WWS
#   D: Log-log scale — Adjusted cases ~ WWS
#
# Omicron weeks (Dec 27, 2021 onward) highlighted in orange to show
# how testing collapse drives divergence between reported and adjusted cases.
#
# Requires: df (merged_wastewater_cases.csv, already loaded and log-transformed)
#           lm_rep_linear, lm_adj_linear, lm_rep_log, lm_adj_log (already fitted)
# =============================================================================

library(ggplot2)
library(dplyr)
library(patchwork)

# --- Define Omicron period ---------------------------------------------------
omicron_start <- as.Date("2021-12-27")

df_plot <- df %>%
  mutate(period = ifelse(date >= omicron_start, "Omicron", "Pre-Omicron"))

# --- Colour palette ----------------------------------------------------------
col_reported  <- "#2171B5"   # blue  — reported cases
col_adjusted  <- "#238B45"   # green — adjusted cases
col_omicron   <- "#E6550D"   # orange — Omicron weeks
col_line      <- "#CC0000"   # red regression line
alpha_main    <- 0.65
alpha_omicron <- 0.90
size_main     <- 1.8
size_omicron  <- 2.2

# --- Helper: R² label --------------------------------------------------------
r2_label <- function(mod) {
  r2 <- round(summary(mod)$r.squared, 3)
  paste0("R\u00B2 = ", r2)
}

# --- Helper: build one scatter panel -----------------------------------------
make_panel <- function(data, xvar, yvar, xlab, ylab,
                       point_colour, lm_mod, tag_label,
                       x_breaks = NULL, y_breaks = NULL,
                       x_labels = NULL, y_labels = NULL) {

  # Regression line data
  x_seq  <- seq(min(data[[xvar]], na.rm = TRUE),
                max(data[[xvar]], na.rm = TRUE), length.out = 200)
  pred   <- predict(lm_mod,
                    newdata = setNames(data.frame(x_seq), xvar),
                    interval = "confidence")
  line_df <- data.frame(x = x_seq, fit = pred[,"fit"],
                        lwr = pred[,"lwr"], upr = pred[,"upr"])

  p <- ggplot(data, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    # Confidence band
    geom_ribbon(data = line_df, aes(x = x, ymin = lwr, ymax = upr),
                inherit.aes = FALSE,
                fill = col_line, alpha = 0.12) +
    # Regression line
    geom_line(data = line_df, aes(x = x, y = fit),
              inherit.aes = FALSE,
              colour = col_line, linewidth = 0.9) +
    # Pre-Omicron points
    geom_point(data = data %>% filter(period == "Pre-Omicron"),
               colour = point_colour, alpha = alpha_main,
               size = size_main, show.legend = FALSE) +
  
    # Omicron points on top
    geom_point(data = data %>% filter(period == "Omicron"),
               colour = col_omicron, alpha = alpha_omicron,
               size = size_omicron, shape = 17, show.legend = FALSE) +  # triangles for Omicron
    # R² annotation
    annotate("text",
             x = min(data[[xvar]], na.rm = TRUE) +
               0.05 * diff(range(data[[xvar]], na.rm = TRUE)),
             y = max(data[[yvar]], na.rm = TRUE) * 0.97,
             label = r2_label(lm_mod),
             hjust = 0, vjust = 1, size = 3.5,
             colour = col_line, fontface = "italic") +
    labs(x = xlab, y = ylab, tag = tag_label) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position   = "none",  
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(colour = "#eeeeee"),
      axis.title        = element_text(colour = "#555555", size = 9.5),
      plot.tag          = element_text(face = "bold", size = 13),
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA)
    )

  if (!is.null(x_breaks)) p <- p + scale_x_continuous(breaks = x_breaks,
                                                        labels = x_labels)
  if (!is.null(y_breaks)) p <- p + scale_y_continuous(breaks = y_breaks,
                                                        labels = y_labels)
  p
}

# =============================================================================
# BUILD FOUR PANELS
# =============================================================================

# Panel A: Linear — Reported
pa <- make_panel(
  data         = df_plot,
  xvar         = "normalized_wws",
  yvar         = "reported_cases",
  xlab         = "Normalized WWS signal",
  ylab         = "Reported cases",
  point_colour = col_reported,
  lm_mod       = lm_rep_linear,
  tag_label    = "A",
  y_breaks     = c(0, 25000, 50000, 75000, 100000),
  y_labels     = c("0", "25k", "50k", "75k", "100k")
)

# Panel B: Linear — Adjusted
pb <- make_panel(
  data         = df_plot,
  xvar         = "normalized_wws",
  yvar         = "adjusted_cases",
  xlab         = "Normalized WWS signal",
  ylab         = "Adjusted cases",
  point_colour = col_adjusted,
  lm_mod       = lm_adj_linear,
  tag_label    = "B",
  y_breaks     = c(0, 50000, 100000, 150000, 200000),
  y_labels     = c("0", "50k", "100k", "150k", "200k")
)

# Panel C: Log-log — Reported
pc <- make_panel(
  data         = df_plot,
  xvar         = "log_wws",
  yvar         = "log_reported",
  xlab         = "log(WWS signal)",
  ylab         = "log(Reported cases)",
  point_colour = col_reported,
  lm_mod       = lm_rep_log,
  tag_label    = "C"
)

# Panel D: Log-log — Adjusted
pd <- make_panel(
  data         = df_plot,
  xvar         = "log_wws",
  yvar         = "log_adjusted",
  xlab         = "log(WWS signal)",
  ylab         = "log(Adjusted cases)",
  point_colour = col_adjusted,
  lm_mod       = lm_adj_log,
  tag_label    = "D"
)

# =============================================================================
# LEGEND — added directly to panel A
# =============================================================================

pa <- pa +
  annotate("point", x = 0.08, y = 91000,
           colour = "#555555", size = 2.5, shape = 16) +
  annotate("text", x = 0.13, y = 91000,
           label = "Pre-Omicron", hjust = 0, size = 3.2,
           colour = "#444444") +
  annotate("point", x = 0.08, y = 85000,
           colour = col_omicron, size = 2.5, shape = 17) +
  annotate("text", x = 0.13, y = 85000,
           label = "Omicron (Dec 2021 onwards)", hjust = 0, size = 3.2,
           colour = "#444444")
# =============================================================================
# COMBINE AND SAVE — patchwork only, no cowplot
# =============================================================================

fig2 <- (pa + pb) / (pc + pd)

ggsave("figure2.png", fig2, width = 10, height = 9,
       dpi = 300, bg = "white")
ggsave("figure2.pdf", fig2, width = 10, height = 9,
       bg = "white")

cat("Saved: figure2.png and figure2.pdf\n")

# =============================================================================
# Figure 2 caption text (for manuscript)
# =============================================================================
cat("\nFigure 2 caption:\n")
cat(strwrap(paste0(
  "Figure 2. Association between normalized wastewater-based surveillance (WWS) signal ",
  "and COVID-19 case counts, Ontario, July 2020 to August 2022. ",
  "Panels A and B show linear-scale regression of reported and test-adjusted cases ",
  "on the WWS signal, respectively. Panels C and D show the corresponding log-log ",
  "regressions. Circles indicate pre-Omicron weeks; triangles (orange) indicate ",
  "Omicron-period weeks (December 27, 2021 onwards), during which population-level ",
  "diagnostic testing contracted sharply. Shaded bands indicate 95% confidence ",
  "intervals for the regression line. R\u00B2 values are shown for each model."
), width = 100), sep = "\n")

print(fig2)

