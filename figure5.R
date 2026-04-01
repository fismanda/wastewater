# =============================================================================
# Figure 3: Lagged Spearman Correlations and DLNM Lag-Response
# Fisman et al. — Wastewater Validation Paper
# =============================================================================
# Two-panel figure:
#   A: Spearman rank correlation between WWS at lag k and cases (k = 0 to 8)
#      for both reported and test-adjusted cases
#   B: DLNM lag-response curves (lags 0 to 8) on log(RR) scale,
#      at 75th percentile WWS exposure, reported vs adjusted overlaid
#
# Requires: df, dlnm_reported, dlnm_adjusted, cb_wws, percentiles_wws
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# --- Colours -----------------------------------------------------------------
col_reported <- "#2171B5"   # blue
col_adjusted <- "#238B45"   # green

# --- Exposure for DLNM predictions ------------------------------------------
wws_p75 <- percentiles_wws["75%"]
cat("75th percentile WWS:", round(wws_p75, 4), "\n")

# =============================================================================
# PANEL A: Lagged Spearman correlations, lags 0 to 8
# =============================================================================

max_lag <- 8
n       <- nrow(df)

spearman_df <- purrr::map_dfr(0:max_lag, function(k) {
  wws_t  <- df$normalized_wws[1:(n - k)]
  rep_tk <- df$reported_cases[(1 + k):n]
  adj_tk <- df$adjusted_cases[(1 + k):n]

  rho_rep <- cor(wws_t, rep_tk, method = "spearman")
  rho_adj <- cor(wws_t, adj_tk, method = "spearman")

  n_obs <- length(wws_t)
  se    <- 1 / sqrt(n_obs - 3)
  z_rep <- atanh(rho_rep)
  z_adj <- atanh(rho_adj)

  data.frame(
    lag     = k,
    rho_rep = rho_rep,
    rho_adj = rho_adj,
    rep_low = tanh(z_rep - 1.96 * se),
    rep_hi  = tanh(z_rep + 1.96 * se),
    adj_low = tanh(z_adj - 1.96 * se),
    adj_hi  = tanh(z_adj + 1.96 * se)
  )
})

spearman_long <- spearman_df %>%
  dplyr::select(lag, rho_rep, rho_adj, rep_low, rep_hi, adj_low, adj_hi) %>%
  tidyr::pivot_longer(
    cols      = c(rho_rep, rho_adj),
    names_to  = "model",
    values_to = "rho"
  ) %>%
  mutate(
    ci_low = ifelse(model == "rho_rep", rep_low, adj_low),
    ci_hi  = ifelse(model == "rho_rep", rep_hi,  adj_hi),
    model  = ifelse(model == "rho_rep", "Reported cases", "Adjusted cases")
  )

panel_a <- ggplot(spearman_long,
                  aes(x = lag, y = rho, colour = model, fill = model)) +
  geom_hline(yintercept = 0, colour = "#999999", linewidth = 0.5) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_hi),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c("Reported cases" = col_reported,
                                  "Adjusted cases" = col_adjusted),
                      name = NULL) +
  scale_fill_manual(values   = c("Reported cases" = col_reported,
                                  "Adjusted cases" = col_adjusted),
                    name = NULL) +
  scale_x_continuous(breaks = 0:8) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x   = "Lag (weeks; WWS leads cases)",
       y   = "Spearman rank correlation (\u03c1)",
       tag = "A") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position   = c(0.72, 0.88),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.text       = element_text(size = 9),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(colour = "#eeeeee"),
    axis.title        = element_text(colour = "#555555", size = 9.5),
    plot.tag          = element_text(face = "bold", size = 13),
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.background  = element_rect(fill = "white", colour = NA)
  )

# =============================================================================
# PANEL B: DLNM lag-response on log(RR) scale, lags 0 to 8
# =============================================================================

pred_rep_p75 <- crosspred(cb_wws, dlnm_reported,
                           at = wws_p75, cen = 0, bylag = 0.1)
pred_adj_p75 <- crosspred(cb_wws, dlnm_adjusted,
                           at = wws_p75, cen = 0, bylag = 0.1)

extract_lag_df <- function(pred, label) {
  raw  <- colnames(pred$matRRfit)
  lags <- as.numeric(sub("^lag", "", raw))
  keep <- !is.na(lags) & lags >= 0
  data.frame(
    lag    = lags[keep],
    RR     = as.numeric(pred$matRRfit)[keep],
    RR_low = as.numeric(pred$matRRlow)[keep],
    RR_hi  = as.numeric(pred$matRRhigh)[keep],
    model  = label
  )
}

df_rep <- extract_lag_df(pred_rep_p75, "Reported cases")
df_adj <- extract_lag_df(pred_adj_p75, "Adjusted cases")
df_rr  <- bind_rows(df_rep, df_adj)

panel_b <- ggplot(df_rr, aes(x = lag, y = log(RR),
                               colour = model, fill = model)) +
  geom_hline(yintercept = 0, colour = "#999999", linewidth = 0.5) +
  geom_ribbon(aes(ymin = log(RR_low), ymax = log(RR_hi)),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = c("Reported cases" = col_reported,
                                  "Adjusted cases" = col_adjusted),
                      name = NULL) +
  scale_fill_manual(values   = c("Reported cases" = col_reported,
                                  "Adjusted cases" = col_adjusted),
                    name = NULL) +
  scale_x_continuous(breaks = 0:8) +
  scale_y_continuous(breaks = seq(-1, 3, 1)) +
  labs(x   = "Lag (weeks; WWS leads cases)",
       y   = paste0("log(Rate ratio) vs. zero signal\n(75th percentile WWS = ",
                    round(wws_p75, 3), ")"),
       tag = "B") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#eeeeee"),
    axis.title       = element_text(colour = "#555555", size = 9.5),
    plot.tag         = element_text(face = "bold", size = 13),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

# =============================================================================
# COMBINE AND SAVE
# =============================================================================
# Note: parentheses around (panel_a | panel_b) required so plot_annotation
# applies to both panels, not just panel_b

fig3 <- (panel_a | panel_b) +
  plot_annotation(
    caption = stringr::str_wrap(paste0(
      "Panel A: Spearman rank correlation between the normalized wastewater signal ",
      "at lag k and case counts k weeks later, for k = 0 to 8 weeks. ",
      "Shaded bands indicate 95% confidence intervals derived via Fisher z-transformation. ",
      "Panel B: lag-response curves from negative binomial distributed lag nonlinear models ",
      "(DLNM) on the log rate ratio scale, evaluated at the 75th percentile of the ",
      "normalized wastewater signal (WWS = ", round(wws_p75, 3), "). ",
      "Cross-basis natural splines with df = 3 in both exposure and lag dimensions; ",
      "maximum lag = 8 weeks; secular trend modelled with natural spline (df = 2). ",
      "Shaded bands indicate 95% confidence intervals. ",
      "Blue = crude reported cases; green = test-adjusted cases."
    ), width = 130),
    theme = theme(
      plot.caption = element_text(size = 8, colour = "#666666",
                                  hjust = 0, lineheight = 1.3)
    )
  )

ggsave("figure3.png", fig3, width = 12, height = 5.5,
       dpi = 300, bg = "white")
ggsave("figure3.pdf", fig3, width = 12, height = 5.5,
       bg = "white")

cat("Saved: figure3.png and figure3.pdf\n")
print(fig3)

