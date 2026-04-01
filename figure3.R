# =============================================================================
# Figure 1: Wastewater signal, case counts, and test-adjustment ratio
# Fisman et al. — Wastewater Validation Paper
# =============================================================================
# Two-panel figure:
#   Panel A: Normalized WWS signal, reported cases, and test-adjusted cases
#   Panel B: Ratio of test-adjusted to reported cases over time
#
# Requires: df (merged_wastewater_cases.csv, already loaded)
#           ratio_data (weekly ratio file, loaded below)
# =============================================================================

library(ggplot2)
library(dplyr)
library(lubridate)
library(patchwork)
library(scales)

# --- Load ratio data ----------------------------------------------------------
# Expects columns: date (Stata format "19jul2020"), tot_adj_case, tot_cases, ratio
ratio_raw <- read.csv(
  "~/Dropbox/Family Room/Wastewater/Final version of wws with all files/ratio_data.csv",
  stringsAsFactors = FALSE
)

# Parse Stata-format dates
ratio_raw$date <- as.Date(ratio_raw$date, format = "%d%b%Y")

# Restrict to WWS study period only
wws_start <- as.Date("2020-07-19")
wws_end   <- as.Date("2022-08-28")

ratio <- ratio_raw %>%
  filter(date >= wws_start & date <= wws_end)

# --- Variant onset dates -----------------------------------------------------
variants <- data.frame(
  date  = as.Date(c("2021-01-15", "2021-07-05", "2021-12-27")),
  label = c("Alpha", "Delta", "Omicron")
)

# --- Shared x-axis limits ----------------------------------------------------
x_limits <- c(wws_start, wws_end)

# =============================================================================
# PANEL A: WWS signal + reported + adjusted cases
# =============================================================================
# Dual y-axis in ggplot2 requires scaling the secondary axis manually.
# We scale cases to the WWS signal range for visual alignment.

wws_max  <- max(df$normalized_wws, na.rm = TRUE)
case_max <- max(df$adjusted_cases, na.rm = TRUE)
scale_factor <- wws_max / case_max

df_long <- df %>%
  dplyr::select(date, normalized_wws, reported_cases, adjusted_cases) %>%
  mutate(
    reported_scaled = reported_cases * scale_factor,
    adjusted_scaled = adjusted_cases * scale_factor
  )

panel_a <- ggplot(df_long, aes(x = date)) +
  # Variant lines
  geom_vline(data = variants, aes(xintercept = date),
             linetype = "dashed", colour = "#888888",
             linewidth = 0.7, alpha = 0.7) +
  # Variant labels
  geom_text(data = variants, aes(x = date + c(14, 0, 0), 
                                 y = wws_max * 1.02, label = label),
            hjust = 0.5, vjust = 0, size = 3.0,
            colour = "#555555", fontface = "italic") +
  # WWS signal
  geom_line(aes(y = normalized_wws, colour = "WWS signal"),
            linewidth = 0.8, alpha = 0.7) +
  # Reported cases (scaled)
  geom_line(aes(y = reported_scaled, colour = "Reported cases"),
            linewidth = 0.9, linetype = "dashed") +
  # Adjusted cases (scaled)
  geom_line(aes(y = adjusted_scaled, colour = "Adjusted cases"),
            linewidth = 0.9) +
  # Primary y-axis (WWS)
  scale_y_continuous(
    name   = "Normalized WWS signal",
    limits = c(0, wws_max * 1.08),
    expand = c(0, 0),
    sec.axis = sec_axis(
      transform = ~ . / scale_factor / 1000,
      name      = "Cases (thousands)"
    )
  ) +
  scale_x_date(
    limits       = x_limits,
    date_breaks  = "3 months",
    date_labels  = "%b\n%Y",
    expand       = c(0.01, 0.01)
  ) +
  scale_colour_manual(
    values = c(
      "WWS signal"      = "#888888",
      "Reported cases"  = "#2171B5",
      "Adjusted cases"  = "#238B45"
    ),
    breaks = c("WWS signal", "Reported cases", "Adjusted cases")
  ) +
  labs(x = NULL, colour = NULL,
       tag = "A") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = c(0.45, 0.6),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key.width  = unit(1.5, "cm"),
    legend.text       = element_text(size = 9),
    plot.margin      = margin(t = 10, r = 5, b = 5, l = 5), 
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(colour = "#eeeeee"),
    axis.title.y.left  = element_text(colour = "#555555", size = 10),
    axis.title.y.right = element_text(colour = "#555555", size = 10),
    plot.tag           = element_text(face = "bold", size = 13),
    plot.background    = element_rect(fill = "white", colour = NA),
    panel.background   = element_rect(fill = "white", colour = NA)
  )

# =============================================================================
# PANEL B: Test-adjustment ratio over time
# =============================================================================

# Identify sub-unity period for annotation
sub_unity <- ratio %>% filter(ratio < 1)

panel_b <- ggplot(ratio, aes(x = date, y = ratio)) +
  # Variant lines
  geom_vline(data = variants, aes(xintercept = date),
             linetype = "dashed", colour = "#888888",
             linewidth = 0.7, alpha = 0.7) +
  # Reference line at ratio = 1
  geom_hline(yintercept = 1, linetype = "solid",
             colour = "#999999", linewidth = 0.6) +
  # Shade sub-unity region
  geom_ribbon(
    data = ratio %>% mutate(ymin = pmin(ratio, 1), ymax = 1),
    aes(ymin = ymin, ymax = ymax),
    fill = "#2171B5", alpha = 0.5
  ) +
  # Ratio line
  geom_line(colour = "#6A3D9A", linewidth = 1.0) +
  # Annotate ratio=1 line
  annotate("text", x = as.Date("2020-10-01"), y = 1.4,
           label = "Ratio = 1", size = 3.0,
           colour = "#999999", hjust = 0) +
  scale_y_continuous(
    name   = "Adjusted : reported case ratio",
    limits = c(0, NA),
    breaks = c(0, 2, 4, 6, 8, 10, 12),
    expand = expansion(mult = c(0, 0.05)),
    labels = number_format(accuracy = 1)
  ) +
  scale_x_date(
    limits      = x_limits,
    date_breaks = "3 months",
    date_labels = "%b\n%Y",
    expand      = c(0.01, 0.01)
  ) +
  labs(x = "Date", tag = "B") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#eeeeee"),
    axis.title.y     = element_text(colour = "#555555", size = 10),
    axis.title.x     = element_text(colour = "#555555", size = 10,
                                    margin = margin(t = 6)),
    plot.tag         = element_text(face = "bold", size = 13),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

# =============================================================================
# COMBINE AND SAVE
# =============================================================================

library(stringr)

caption_text <- paste0(
  "Dashed vertical lines indicate approximate onset of Alpha (Jan 2021), ",
  "Delta (Jul 2021), and Omicron (Dec 2021) variant waves in Ontario. ",
  "Shaded region in panel B indicates weeks where adjusted cases fell ",
  "below reported cases. WWS = wastewater-based surveillance."
)

fig1 <- panel_a / panel_b +
  plot_layout(heights = c(1.4, 1)) +
  plot_annotation(
    caption = str_wrap(caption_text, width = 120),
    theme = theme(
      plot.caption = element_text(size = 8, colour = "#666666",
                                  hjust = 0, lineheight = 1.3)
    )
  )

ggsave("figure1.png", fig1, width = 10, height = 8,
       dpi = 300, bg = "white")
ggsave("figure1.pdf", fig1, width = 10, height = 8,
       bg = "white")

cat("Saved: figure1.png and figure1.pdf\n")


print(fig1)

