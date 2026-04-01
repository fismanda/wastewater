# =============================================================================
# Wastewater Validation Analysis — Fisman et al.
# =============================================================================
# Data: merged_wastewater_cases.csv
# Columns:
#   date            — week start date (case reporting week)
#   series_week     — sequential week number (29–139)
#   adjusted_cases  — test-adjusted case estimates (Bosco et al. method)
#   reported_cases  — crude reported COVID-19 cases
#   tests           — number of tests performed
#   normalized_wws  — SD-normalized wastewater signal (week_norm_mn)
#   wws_week_start  — start of matched wastewater measurement period
#   wws_week_end    — end of matched wastewater measurement period
#
# Sections:
#   1. Packages and data
#   2. Spearman rank correlations
#   3. Linear regression R² (linear-linear and log-log)
#   4. Distributed lag nonlinear models (DLNM)
#      4a. Model fit and pseudo-R²
#      4b. Lag-response plots
#      4c. 3D perspective and contour plots (base R)
#      4d. ggplot2 contour plots (publication quality)
#      4e. Plotly interactive 3D surfaces
#      4f. Cumulative RR table at percentiles
#   5. Granger causality
#   6. Time series overview plot
#
# Install packages if needed:
# install.packages(c("tidyverse","dlnm","tseries","lmtest","vars",
#                    "patchwork","viridis","plotly","htmlwidgets","MASS"))
# =============================================================================


# =============================================================================
# SECTION 1: PACKAGES AND DATA
# =============================================================================

library(tidyverse)
library(dlnm)
library(splines)
library(tseries)
library(lmtest)
library(vars)
library(patchwork)
library(viridis)
library(plotly)
library(htmlwidgets)
library(MASS)

# --- Load data ---------------------------------------------------------------
df <- read.csv(
  "~/Dropbox/Family Room/Wastewater/Final version of wws with all files/merged_wastewater_cases.csv",
  stringsAsFactors = FALSE
)
df$date <- as.Date(df$date)

cat("Rows:", nrow(df), "\n")
cat("Weeks:", min(df$series_week), "to", max(df$series_week), "\n")
cat("Dates:", as.character(min(df$date)), "to", as.character(max(df$date)), "\n")
cat("Any missing:", any(is.na(df)), "\n")

# Small offset to avoid log(0) for WWS values of zero in early weeks
wws_offset <- min(df$normalized_wws[df$normalized_wws > 0]) / 2

df <- df %>%
  mutate(
    log_wws      = log(normalized_wws + wws_offset),
    log_reported = log(reported_cases),
    log_adjusted = log(adjusted_cases)
  )


# =============================================================================
# SECTION 2: SPEARMAN RANK CORRELATIONS
# =============================================================================
# Spearman is rank-based and robust to skew and monotone nonlinearity.
# Log transformation is neither needed nor helpful — ranks are invariant
# to monotone transformations. We use raw values throughout.

cat("\n--- Spearman Correlations ---\n")

rho_reported <- cor.test(df$normalized_wws, df$reported_cases,
                         method = "spearman", exact = FALSE)
rho_adjusted <- cor.test(df$normalized_wws, df$adjusted_cases,
                         method = "spearman", exact = FALSE)

cat("WWS vs reported cases:  rho =", round(rho_reported$estimate, 3),
    "  p =", format(rho_reported$p.value, scientific = TRUE, digits = 3), "\n")
cat("WWS vs adjusted cases:  rho =", round(rho_adjusted$estimate, 3),
    "  p =", format(rho_adjusted$p.value, scientific = TRUE, digits = 3), "\n")
# =============================================================================
# SECTION 2b: STRATIFIED SPEARMAN AND R² — PRE-OMICRON VS OMICRON
# =============================================================================
# We stratify at December 27, 2021 — the onset of the Omicron wave in Ontario,
# coinciding with the sharp contraction in population-level diagnostic testing.
# Pre-Omicron includes weeks with limited WWS geographic coverage (Hamilton and
# Ottawa only) prior to provincial scale-up, which contributes to weaker
# correlations in that period alongside variable testing intensity.

omicron_start <- as.Date("2021-12-27")

df_pre     <- df %>% dplyr::filter(date <  omicron_start)
df_omicron <- df %>% dplyr::filter(date >= omicron_start)

cat("\nStratified analysis:\n")
cat("  Pre-Omicron: n =", nrow(df_pre),
    "(", as.character(min(df_pre$date)), "to", as.character(max(df_pre$date)), ")\n")
cat("  Omicron+:    n =", nrow(df_omicron),
    "(", as.character(min(df_omicron$date)), "to", as.character(max(df_omicron$date)), ")\n")

# Helper: Spearman with Fisher z 95% CI
spearman_ci <- function(x, y) {
  r   <- cor.test(x, y, method = "spearman", exact = FALSE)
  rho <- as.numeric(r$estimate)
  n   <- length(x)
  se  <- 1 / sqrt(n - 3)
  z   <- atanh(rho)
  list(rho = rho, p = r$p.value,
       low = tanh(z - 1.96 * se),
       hi  = tanh(z + 1.96 * se))
}

# Print results for each period
for (period_label in c("Full", "Pre-Omicron", "Omicron+")) {
  d <- switch(period_label,
              "Full"        = df,
              "Pre-Omicron" = df_pre,
              "Omicron+"    = df_omicron)

  s_rep <- spearman_ci(d$normalized_wws, d$reported_cases)
  s_adj <- spearman_ci(d$normalized_wws, d$adjusted_cases)
  r2_rep <- summary(lm(reported_cases ~ normalized_wws, data = d))$r.squared
  r2_adj <- summary(lm(adjusted_cases ~ normalized_wws, data = d))$r.squared

  cat(sprintf("\n%s (n=%d):\n", period_label, nrow(d)))
  cat(sprintf("  Reported: rho = %.3f (95%% CI %.3f-%.3f)  R2 = %.3f\n",
              s_rep$rho, s_rep$low, s_rep$hi, r2_rep))
  cat(sprintf("  Adjusted: rho = %.3f (95%% CI %.3f-%.3f)  R2 = %.3f\n",
              s_adj$rho, s_adj$low, s_adj$hi, r2_adj))
}


# =============================================================================
# SECTION 3: LINEAR REGRESSION R²
# =============================================================================
# We compare two specifications:
#
#   (A) Linear-linear: cases ~ wws (raw scale)
#       R² = proportion of variance in raw cases explained by raw WWS.
#       The Omicron testing collapse drives this down for reported cases.
#
#   (B) Log-log: log(cases) ~ log(wws)
#       Appropriate when the relationship is multiplicative / power-law.
#       R² = proportion of variance in LOG cases explained by LOG WWS.
#       NOTE: log-scale R² is not comparable to linear-scale R².
#       Log transformation compresses the Omicron peak, partially
#       "rescuing" the reported cases fit — interpret with caution.

cat("\n--- Linear Regression R² ---\n")

lm_rep_linear <- lm(reported_cases ~ normalized_wws, data = df)
lm_adj_linear <- lm(adjusted_cases ~ normalized_wws, data = df)
lm_rep_log    <- lm(log_reported   ~ log_wws,        data = df)
lm_adj_log    <- lm(log_adjusted   ~ log_wws,        data = df)

cat("\n(A) Linear scale (cases ~ wws):\n")
cat("  Reported cases  R² =", round(summary(lm_rep_linear)$r.squared, 4), "\n")
cat("  Adjusted cases  R² =", round(summary(lm_adj_linear)$r.squared, 4), "\n")

cat("\n(B) Log-log scale (log(cases) ~ log(wws)):\n")
cat("  Reported cases  R² =", round(summary(lm_rep_log)$r.squared, 4), "\n")
cat("  Adjusted cases  R² =", round(summary(lm_adj_log)$r.squared, 4), "\n")

# --- Scatter plots -----------------------------------------------------------
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

plot(df$normalized_wws, df$reported_cases,
     main = "Linear: Reported ~ WWS", xlab = "WWS signal",
     ylab = "Reported cases", pch = 16, col = "steelblue", cex = 0.7)
abline(lm_rep_linear, col = "red", lwd = 2)
legend("topleft", bty = "n", text.col = "red",
       legend = paste0("R² = ", round(summary(lm_rep_linear)$r.squared, 3)))

plot(df$normalized_wws, df$adjusted_cases,
     main = "Linear: Adjusted ~ WWS", xlab = "WWS signal",
     ylab = "Adjusted cases", pch = 16, col = "darkgreen", cex = 0.7)
abline(lm_adj_linear, col = "red", lwd = 2)
legend("topleft", bty = "n", text.col = "red",
       legend = paste0("R² = ", round(summary(lm_adj_linear)$r.squared, 3)))

plot(df$log_wws, df$log_reported,
     main = "Log-log: Reported ~ WWS", xlab = "log(WWS signal)",
     ylab = "log(Reported cases)", pch = 16, col = "steelblue", cex = 0.7)
abline(lm_rep_log, col = "red", lwd = 2)
legend("topleft", bty = "n", text.col = "red",
       legend = paste0("R² = ", round(summary(lm_rep_log)$r.squared, 3)))

plot(df$log_wws, df$log_adjusted,
     main = "Log-log: Adjusted ~ WWS", xlab = "log(WWS signal)",
     ylab = "log(Adjusted cases)", pch = 16, col = "darkgreen", cex = 0.7)
abline(lm_adj_log, col = "red", lwd = 2)
legend("topleft", bty = "n", text.col = "red",
       legend = paste0("R² = ", round(summary(lm_adj_log)$r.squared, 3)))

par(mfrow = c(1, 1))


# =============================================================================
# SECTION 4: DISTRIBUTED LAG NONLINEAR MODELS (DLNM)
# =============================================================================
# DLNMs allow the effect of WWS to be distributed across multiple lags
# AND allow the exposure-response to be nonlinear.
#
# Model structure:
#   - Outcome: raw case counts (negative binomial handles count scale + overdispersion)
#   - Cross-basis for WWS: natural spline (df=3) in both exposure and lag dimensions
#   - Maximum lag: 8 weeks
#   - Secular trend: natural spline on week index (df=2)
#     df=2 chosen to avoid competing with WWS cross-basis for variance
#
# Pseudo-R² (McFadden) = 1 - (residual deviance / null deviance)

cat("\n--- DLNM Models ---\n")

max_lag  <- 8
week_idx <- df$series_week - min(df$series_week)   # centred week index

cb_wws <- crossbasis(
  df$normalized_wws,
  lag    = max_lag,
  argvar = list(fun = "ns", df = 3),
  arglag = list(fun = "ns", df = 3)
)

# --- 4a. Model fit and pseudo-R² ---------------------------------------------
# We use negative binomial regression rather than quasi-Poisson because:
#   - NB has a proper likelihood, enabling AIC-based model comparison
#   - NB explicitly models overdispersion via the theta parameter
#     (lower theta = more overdispersion)
#   - quasi-Poisson uses an approximation; NB is more principled when
#     overdispersion is severe, as it clearly is here
# The cross-basis and crosspred() calls are identical — only the family changes.


dlnm_reported <- MASS::glm.nb(
  reported_cases ~ cb_wws + splines::ns(week_idx, df = 2),
  data = df
)

dlnm_adjusted <- MASS::glm.nb(
  adjusted_cases ~ cb_wws + splines::ns(week_idx, df = 2),
  data = df
)

# Trend-only null models (no WWS) — to isolate WWS contribution
dlnm_trend_rep <- MASS::glm.nb(
  reported_cases ~ splines::ns(week_idx, df = 2),
  data = df
)

dlnm_trend_adj <- MASS::glm.nb(
  adjusted_cases ~ splines::ns(week_idx, df = 2),
  data = df
)

# McFadden pseudo-R² = 1 - (residual deviance / null deviance)
pseudo_r2 <- function(mod) 1 - (mod$deviance / mod$null.deviance)

cat("\n                         Reported    Adjusted\n")
cat("Theta (overdispersion):  ",
    round(dlnm_reported$theta, 3), "      ",
    round(dlnm_adjusted$theta, 3), "\n")
cat("Trend-only pseudo-R²:   ",
    round(pseudo_r2(dlnm_trend_rep), 4), "     ",
    round(pseudo_r2(dlnm_trend_adj), 4), "\n")
cat("Full model pseudo-R²:   ",
    round(pseudo_r2(dlnm_reported), 4), "     ",
    round(pseudo_r2(dlnm_adjusted), 4), "\n")
cat("WWS increment:          ",
    round(pseudo_r2(dlnm_reported) - pseudo_r2(dlnm_trend_rep), 4), "     ",
    round(pseudo_r2(dlnm_adjusted) - pseudo_r2(dlnm_trend_adj), 4), "\n")
cat("AIC:                    ",
    round(AIC(dlnm_reported), 1), "  ",
    round(AIC(dlnm_adjusted), 1), "\n")
cat("AIC (trend only):       ",
    round(AIC(dlnm_trend_rep), 1), "  ",
    round(AIC(dlnm_trend_adj), 1), "\n")

# Likelihood ratio test — does adding WWS significantly improve fit?
# (Valid because NB has a proper likelihood, unlike quasi-Poisson)
lrt_rep <- anova(dlnm_trend_rep, dlnm_reported, test = "Chisq")
lrt_adj <- anova(dlnm_trend_adj, dlnm_adjusted, test = "Chisq")
cat("\nLRT p-value (WWS vs trend-only):\n")
cat("  Reported cases: p =", format(lrt_rep[2, "Pr(Chi)"], scientific = TRUE, digits = 3), "\n")
cat("  Adjusted cases: p =", format(lrt_adj[2, "Pr(Chi)"], scientific = TRUE, digits = 3), "\n")

# --- 4b. Lag-response plots (at median WWS exposure) -------------------------

pred_lag_rep <- crosspred(cb_wws, dlnm_reported,
                          at = median(df$normalized_wws), cen = 0, cumul = TRUE)
pred_lag_adj <- crosspred(cb_wws, dlnm_adjusted,
                          at = median(df$normalized_wws), cen = 0, cumul = TRUE)

par(mfrow = c(1, 2))
plot(pred_lag_rep, "slices", var = median(df$normalized_wws),
     main = "Lag-response: Reported cases",
     xlab = "Lag (weeks)", ylab = "RR vs. zero signal",
     col = "steelblue", lwd = 2,
     ci.arg = list(col = adjustcolor("steelblue", 0.2)))

plot(pred_lag_adj, "slices", var = median(df$normalized_wws),
     main = "Lag-response: Adjusted cases",
     xlab = "Lag (weeks)", ylab = "RR vs. zero signal",
     col = "darkgreen", lwd = 2,
     ci.arg = list(col = adjustcolor("darkgreen", 0.2)))
par(mfrow = c(1, 1))

# Note: bidirectional DLNM (lags -4 to +8) was explored as a sensitivity
# analysis but is not reported in the primary manuscript. The lag-response
# peak at lag 0 for both series confirmed contemporaneous tracking; negative
# lags did not improve model fit. See manuscript for discussion.

# --- 4c. Full prediction objects for surface plots ---------------------------
# Two prediction objects:
#   pred_rep / pred_adj — fine lag grid (bylag=0.1) for smooth base R plots,
#                         includes percentiles for slice/overall plots
#   pred_gg_rep / pred_gg_adj — integer lags (bylag=1) for ggplot2 tiles

percentiles_wws <- round(quantile(df$normalized_wws,
                                  c(0.05, 0.25, 0.50, 0.75, 0.95)), 4)
wws_seq <- seq(min(df$normalized_wws), max(df$normalized_wws), length.out = 100)

cat("\nWWS percentiles:\n")
print(percentiles_wws)

pred_rep <- crosspred(cb_wws, dlnm_reported,
                      at = c(percentiles_wws, wws_seq),
                      cen = 0, bylag = 0.1, cumul = TRUE)

pred_adj <- crosspred(cb_wws, dlnm_adjusted,
                      at = c(percentiles_wws, wws_seq),
                      cen = 0, bylag = 0.1, cumul = TRUE)

pred_gg_rep <- crosspred(cb_wws, dlnm_reported,
                         at = wws_seq, cen = 0, bylag = 1)

pred_gg_adj <- crosspred(cb_wws, dlnm_adjusted,
                         at = wws_seq, cen = 0, bylag = 1)

# --- 4d. Base R surface plots → saved to PDF ---------------------------------

pdf("dlnm_surfaces.pdf", width = 10, height = 7)

# 3D perspective
plot(pred_rep, cumul = FALSE, xlab = "WWS signal", zlab = "RR",
     main = "3D surface: WWS vs. Reported cases",
     theta = 200, phi = 35, lphi = 30, col = "lightblue")

plot(pred_adj, cumul = FALSE, xlab = "WWS signal", zlab = "RR",
     main = "3D surface: WWS vs. Adjusted cases",
     theta = 200, phi = 35, lphi = 30, col = "lightgreen")

# Contour — use palette= not col= to avoid argument conflict
plot(pred_rep, ptype = "contour", cumul = FALSE,
     key.title  = title("RR"),
     plot.title = title("Contour: WWS vs. Reported cases",
                        xlab = "WWS signal", ylab = "Lag (weeks)"),
     palette    = topo.colors(12))

plot(pred_adj, ptype = "contour", cumul = FALSE,
     key.title  = title("RR"),
     plot.title = title("Contour: WWS vs. Adjusted cases",
                        xlab = "WWS signal", ylab = "Lag (weeks)"),
     palette    = topo.colors(12))

# Slice plots — lag-response at 75th and 95th percentile exposure
for (pct in c("75%", "95%")) {
  plot(pred_rep, "slices", var = percentiles_wws[pct], cumul = FALSE,
       ylab = "RR vs. zero signal",
       main = paste0("Lag-response at ", pct, " WWS (",
                     round(percentiles_wws[pct], 3), "): Reported"),
       col = "steelblue", lwd = 2,
       ci.arg = list(col = adjustcolor("steelblue", 0.2)))

  plot(pred_adj, "slices", var = percentiles_wws[pct], cumul = FALSE,
       ylab = "RR vs. zero signal",
       main = paste0("Lag-response at ", pct, " WWS (",
                     round(percentiles_wws[pct], 3), "): Adjusted"),
       col = "darkgreen", lwd = 2,
       ci.arg = list(col = adjustcolor("darkgreen", 0.2)))
}

# Overall cumulative association
plot(pred_rep, "overall", cumul = TRUE,
     ylab = "Cumulative RR vs. zero signal", xlab = "WWS signal",
     main = "Overall cumulative association: Reported cases",
     col = "steelblue", lwd = 2,
     ci.arg = list(col = adjustcolor("steelblue", 0.2)))
abline(h = 1, lty = 2, col = "grey50")

plot(pred_adj, "overall", cumul = TRUE,
     ylab = "Cumulative RR vs. zero signal", xlab = "WWS signal",
     main = "Overall cumulative association: Adjusted cases",
     col = "darkgreen", lwd = 2,
     ci.arg = list(col = adjustcolor("darkgreen", 0.2)))
abline(h = 1, lty = 2, col = "grey50")

dev.off()
cat("Saved: dlnm_surfaces.pdf\n")

# --- 4e. ggplot2 contour plots -----------------------------------------------

make_surface_df <- function(pred, wws_vals, label) {
  rr_mat <- pred$matRRfit
  colnames(rr_mat) <- seq(0, ncol(rr_mat) - 1)   # force clean 0,1,2...8
  as.data.frame(rr_mat) %>%
    mutate(exposure = wws_vals) %>%
    pivot_longer(-exposure, names_to = "lag", values_to = "RR") %>%
    mutate(lag = as.numeric(lag), log_RR = log(RR)) %>%
    filter(is.finite(log_RR)) %>%
    mutate(model = label)
}

df_gg_rep <- make_surface_df(pred_gg_rep, wws_seq, "Reported cases")
df_gg_adj <- make_surface_df(pred_gg_adj, wws_seq, "Adjusted cases")
df_gg_all <- bind_rows(df_gg_rep, df_gg_adj)

clim <- min(max(abs(df_gg_all$log_RR), na.rm = TRUE), 3)

make_gg_contour <- function(data, title_text, clim) {
  ggplot(data, aes(x = exposure, y = lag)) +
    geom_tile(aes(fill = pmin(pmax(log_RR, -clim), clim))) +
    geom_contour(aes(z = log_RR), colour = "white",
                 alpha = 0.35, linewidth = 0.25, breaks = c(-2, -1, 0, 1, 2)) +
    geom_contour(aes(z = log_RR), colour = "black",
                 linewidth = 0.7, breaks = 0) +
    scale_fill_distiller(
      palette   = "RdYlBu",
      direction = -1,
      limits    = c(-clim, clim),
      name      = "RR",
      labels    = function(x) round(exp(x), 1),
      breaks    = c(-clim, -clim/2, 0, clim/2, clim),
      guide     = guide_colorbar(barwidth = 1, barheight = 7)
    ) +
    scale_x_continuous(expand = c(0, 0),
                       labels = scales::number_format(accuracy = 0.1)) +
    scale_y_continuous(breaks = seq(0, 8, 2), expand = c(0, 0)) +
    labs(title = title_text,
         x = "Normalized wastewater signal", y = "Lag (weeks)") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(size = 12, face = "bold", hjust = 0.5),
      panel.grid       = element_blank(),
      legend.position  = "right",
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

p_gg_rep <- make_gg_contour(df_gg_rep, "Reported cases", clim)
p_gg_adj <- make_gg_contour(df_gg_adj, "Adjusted cases", clim)

p_gg_combined <- p_gg_rep + p_gg_adj +
  plot_annotation(
    title    = "DLNM exposure-lag-response surface: wastewater signal vs. COVID-19 cases",
    subtitle = "Colour = rate ratio vs. zero WWS signal; black contour = RR 1.0 (null)",
    caption  = "Cross-basis: natural spline (df=3 exposure, df=3 lag, max lag=8 weeks); negative binomial GLM with df=2 secular trend",
    theme    = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, colour = "#666666"),
      plot.caption  = element_text(size =  9, colour = "#888888")
    )
  )

ggsave("dlnm_contour_ggplot.png", p_gg_combined, width = 12, height = 5,
       dpi = 300, bg = "white")
ggsave("dlnm_contour_ggplot.pdf", p_gg_combined, width = 12, height = 5,
       bg = "white")
cat("Saved: dlnm_contour_ggplot.png and .pdf\n")

# --- 4f. Plotly interactive 3D surfaces → saved as HTML ----------------------

wws_seq_50 <- seq(min(df$normalized_wws), max(df$normalized_wws), length.out = 50)
lag_seq    <- 0:8

pred_plotly_rep <- crosspred(cb_wws, dlnm_reported,
                             at = wws_seq_50, cen = 0, bylag = 1)
pred_plotly_adj <- crosspred(cb_wws, dlnm_adjusted,
                             at = wws_seq_50, cen = 0, bylag = 1)

make_plotly_surface <- function(pred, title_text) {
  rr_mat <- pmin(t(pred$matRRfit), 20)   # rows=lags, cols=exposures; cap at 20
  plot_ly(
    x = ~wws_seq_50, y = ~lag_seq, z = ~rr_mat,
    type = "surface", colorscale = "RdBu", reversescale = TRUE,
    colorbar = list(title = "RR")
  ) %>%
    layout(
      title = title_text,
      scene = list(
        xaxis = list(title = "WWS signal"),
        yaxis = list(title = "Lag (weeks)"),
        zaxis = list(title = "Rate ratio")
      )
    )
}

fig_rep_3d <- make_plotly_surface(pred_plotly_rep, "DLNM Surface — Reported Cases")
fig_adj_3d <- make_plotly_surface(pred_plotly_adj, "DLNM Surface — Adjusted Cases")

saveWidget(fig_rep_3d, "dlnm_surface_reported.html", selfcontained = TRUE)
saveWidget(fig_adj_3d, "dlnm_surface_adjusted.html", selfcontained = TRUE)
cat("Saved: dlnm_surface_reported.html and dlnm_surface_adjusted.html\n")

# --- 4g. Cumulative RR table at percentiles ----------------------------------

cat("\n--- Cumulative RR at WWS percentiles (vs. zero signal) ---\n")
cat("\nReported cases:\n")
cat("  RR:      "); print(round(pred_rep$allRRfit[as.character(percentiles_wws)], 3))
cat("  95% low: "); print(round(pred_rep$allRRlow[as.character(percentiles_wws)], 3))
cat("  95% high:"); print(round(pred_rep$allRRhigh[as.character(percentiles_wws)], 3))

cat("\nAdjusted cases:\n")
cat("  RR:      "); print(round(pred_adj$allRRfit[as.character(percentiles_wws)], 3))
cat("  95% low: "); print(round(pred_adj$allRRlow[as.character(percentiles_wws)], 3))
cat("  95% high:"); print(round(pred_adj$allRRhigh[as.character(percentiles_wws)], 3))


# =============================================================================
# SECTION 5: GRANGER CAUSALITY
# =============================================================================
# Tests whether past values of WWS improve prediction of cases beyond
# what past cases alone provide ("Granger-cause").
# We work on log-differenced series:
#   - Log scale compresses the Omicron spike for linear VAR
#   - First differences achieve stationarity (verified by ADF test)
#   - H0 for ADF: series has a unit root (non-stationary)

cat("\n--- Granger Causality ---\n")

log_wws      <- log(df$normalized_wws + wws_offset)
log_reported <- log(df$reported_cases)
log_adjusted <- log(df$adjusted_cases)

# ADF tests — levels
cat("\nADF tests — levels (H0: non-stationary):\n")
for (nm in c("log(WWS)", "log(reported)", "log(adjusted)")) {
  s <- list(log_wws, log_reported, log_adjusted)[[match(nm, c("log(WWS)", "log(reported)", "log(adjusted)"))]]
  res <- adf.test(s, alternative = "stationary")
  cat(sprintf("  %-16s p = %.4f  %s\n", nm,
              res$p.value,
              ifelse(res$p.value < 0.05, "→ stationary", "→ non-stationary")))
}

# First differences
d_log_wws      <- diff(log_wws)
d_log_reported <- diff(log_reported)
d_log_adjusted <- diff(log_adjusted)

cat("\nADF tests — first differences:\n")
for (nm in c("Δlog(WWS)", "Δlog(reported)", "Δlog(adjusted)")) {
  s <- list(d_log_wws, d_log_reported, d_log_adjusted)[[match(nm, c("Δlog(WWS)", "Δlog(reported)", "Δlog(adjusted)"))]]
  res <- adf.test(s, alternative = "stationary")
  cat(sprintf("  %-16s p = %.4f  %s\n", nm,
              res$p.value,
              ifelse(res$p.value < 0.05, "→ stationary", "→ non-stationary")))
}

# Granger tests — forward and reverse directions
run_granger <- function(y, x, label) {
  cat("\n---", label, "---\n")
  for (k in 1:8) {
    gt <- grangertest(y ~ x, order = k)
    cat(sprintf("  Lag %d: F = %.3f, df = (%d,%d), p = %.4f %s\n",
                k, gt$F[2], gt$Df[2], gt$Res.Df[2], gt$`Pr(>F)`[2],
                ifelse(gt$`Pr(>F)`[2] < 0.05, "*", "")))
  }
}

run_granger(d_log_reported, d_log_wws, "WWS → Reported cases")
run_granger(d_log_adjusted, d_log_wws, "WWS → Adjusted cases")
run_granger(d_log_wws, d_log_reported, "Reported cases → WWS (reverse)")
run_granger(d_log_wws, d_log_adjusted, "Adjusted cases → WWS (reverse)")


# =============================================================================
# SECTION 6: TIME SERIES OVERVIEW PLOT
# =============================================================================

par(mar = c(4, 5, 2, 5))

plot(df$date, df$normalized_wws,
     type = "l", col = "grey50", lwd = 1.5,
     xlab = "Date", ylab = "Normalized WWS signal",
     main = "Wastewater signal vs. case counts over time",
     ylim = c(0, max(df$normalized_wws) * 1.1))

par(new = TRUE)
plot(df$date, df$reported_cases / 1000,
     type = "l", col = "steelblue", lwd = 1.5, lty = 2,
     axes = FALSE, xlab = "", ylab = "",
     ylim = c(0, max(df$adjusted_cases) / 1000 * 1.1))
lines(df$date, df$adjusted_cases / 1000, col = "darkgreen", lwd = 1.5)

axis(side = 4)
mtext("Cases (thousands)", side = 4, line = 3)
legend("topleft", bty = "n", lwd = 1.5, cex = 0.85,
       legend = c("WWS signal", "Reported cases", "Adjusted cases"),
       col    = c("grey50", "steelblue", "darkgreen"),
       lty    = c(1, 2, 1))

cat("\nAnalysis complete.\n")

