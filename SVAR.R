# ==============================
# 0) Libraries
# ==============================
library(ggplot2)
library(scales)
library(vars)
library(svars)
library(urca)
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(writexl)

# ==============================
# 0.5) Settings
# ==============================

sample_tag <- "full"  # "full", "subsample1", "subsample2"
use_regime_dummies <- TRUE

data_file <- switch(
  sample_tag,
  "full"       = "Full_sample.xlsx",
  "subsample1" = "Sample_2014_2025.xlsx",
  "subsample2" = "Sample_2020_2025.xlsx"
)

# ==============================
# 1) File paths and data
# ==============================

base    <- "C:/Users/khw/Desktop/SVAR_seminar"
in_dir  <- file.path(base, "Data")
fig_dir <- file.path(base, "Figures")
tab_dir <- file.path(base, "Tables")

fp_master <- file.path(in_dir, data_file)
D <- read_excel(fp_master, sheet = 1)

# Date column
if (inherits(D$date, "Date")) {
  # ok
} else if (is.numeric(D$date)) {
  D$date <- as.Date(D$date, origin = "1899-12-30")
} else {
  D$date <- as.Date(D$date)
}

D <- D %>% arrange(date)

# ==============================
# 2) Regime dummies
# ==============================

if (use_regime_dummies) {
  D <- D %>%
    mutate(
      D_covid = if_else(
        date >= as.Date("2020-03-01") & date <= as.Date("2021-12-31"),
        1L, 0L
      ),
      D_hike = if_else(date >= as.Date("2022-07-01"), 1L, 0L)
    )
}

# ==============================
# 3) Stationarity tests
# ==============================

series <- c("opr", "dm3", "pi", "y", "y2", "y5", "y10")

# ADF
adf_list <- lapply(series, function(s) {
  ur.df(na.omit(D[[s]]), type = "drift", lags = 12)
})
names(adf_list) <- series

adf_out_path <- file.path(tab_dir, sprintf("diagnostics_adf_%s.txt", sample_tag))
sink(adf_out_path)
cat("====================================\n")
cat("ADF TEST RESULTS - sample:", sample_tag, "\n")
cat("====================================\n\n")
for (s in series) {
  cat("\nADF test for", s, "\n")
  print(summary(adf_list[[s]]))
}
sink()
message("Saved ADF test results to: ", adf_out_path)

# KPSS
kpss_list <- lapply(series, function(s) {
  ur.kpss(na.omit(D[[s]]), type = "mu")
})
names(kpss_list) <- series

kpss_out_path <- file.path(tab_dir, sprintf("diagnostics_kpss_%s.txt", sample_tag))
sink(kpss_out_path)
cat("====================================\n")
cat("KPSS TEST RESULTS - sample:", sample_tag, "\n")
cat("====================================\n\n")
for (s in series) {
  cat("\nKPSS test for", s, "\n")
  print(summary(kpss_list[[s]]))
}
sink()
message("Saved KPSS test results to: ", kpss_out_path)

# ==============================
# 4) VAR dataset
# ==============================

vars_y <- c("y", "pi", "opr", "dm3", "y2", "y5", "y10")

if (use_regime_dummies) {
  Z <- D %>% select(date, all_of(vars_y), D_covid, D_hike) %>% drop_na()
  X_exog <- Z %>% select(D_covid, D_hike)
} else {
  Z <- D %>% select(date, all_of(vars_y)) %>% drop_na()
  X_exog <- NULL
}

Y_var <- Z %>% select(all_of(vars_y))

# Lag selection
if (use_regime_dummies) {
  sel <- VARselect(Y_var, lag.max = 6, type = "const", exogen = X_exog)$selection
} else {
  sel <- VARselect(Y_var, lag.max = 6, type = "const")$selection
}

lag_candidates <- sel[c("AIC(n)", "HQ(n)", "SC(n)")]
message("Lag selected by criteria:")
print(lag_candidates)

if (all(is.na(lag_candidates))) {
  var_lag <- 2L
  message("Using lag 2.")
} else {
  var_lag <- as.integer(stats::median(lag_candidates, na.rm = TRUE))
  message(sprintf("Chosen VAR lag (median of AIC, HQ, BIC): %s", var_lag))
}

# Estimate VAR (forcing lag=2 based on aic, bic and irfs)
if (use_regime_dummies) {
  v <- VAR(Y_var, p = 2, type = "const", exogen = X_exog)
} else {
  v <- VAR(Y_var, p = 2, type = "const")
}

bic_val_model <- BIC(v)
aic_val_model <- AIC(v)

# ==============================
# 5) Johansen cointegration
# ==============================

joh_data <- as.matrix(Y_var)
K_joh <- max(2L, var_lag)

joh_trace <- ca.jo(joh_data, type = "trace", ecdet = "const", K = K_joh)
joh_eigen <- ca.jo(joh_data, type = "eigen", ecdet = "const", K = K_joh)

joh_out_path <- file.path(tab_dir, sprintf("johansen_results_yields_%s.txt", sample_tag))
sink(joh_out_path)
cat("====================================\n")
cat("JOHANSEN COINTEGRATION TEST RESULTS - sample:", sample_tag, "\n")
cat("====================================\n\n")

cat("--- Trace test (type = 'trace', ecdet = 'const') ---\n\n")
print(summary(joh_trace))

cat("\n\n--- Max eigenvalue test (type = 'eigen', ecdet = 'const') ---\n\n")
print(summary(joh_eigen))
sink()
message("Saved Johansen test results to: ", joh_out_path)

# ==============================
# 6) Structural identification (Cholesky)
# ==============================

Sigma <- crossprod(resid(v)) / nrow(resid(v))
eigen(Sigma)$values

ord <- vars_y
s_chol <- id.chol(v, order_k = ord)

cat("\n=== VAR stability (roots of companion) ===\n")
mod <- roots(v, modulus = TRUE)
print(mod)
if (any(mod >= 1)) message("Warning: at least one root >= 1 (VAR not stable).")

# ==============================
# 7) IRFs: responses to OPR shock
# ==============================

set.seed(123)
ir <- irf(
  v,
  impulse = "opr",
  n.ahead = 12,
  boot = TRUE,
  runs = 2000,
  ci = 0.95,
  ortho = TRUE
)

impulses <- names(ir$irf)

irf_df <- map_dfr(impulses, function(imp) {
  ir_mat  <- ir$irf[[imp]]
  low_mat <- ir$Lower[[imp]]
  up_mat  <- ir$Upper[[imp]]

  tibble(
    impulse  = imp,
    step     = rep(0:(nrow(ir_mat) - 1), times = ncol(ir_mat)),
    response = rep(colnames(ir_mat),      each  = nrow(ir_mat)),
    irf      = as.vector(ir_mat),
    lower95  = as.vector(low_mat),
    upper95  = as.vector(up_mat)
  )
})

write_xlsx(
  irf_df,
  path = file.path(tab_dir, sprintf("irf_results_yields_%s.xlsx", sample_tag))
)

plot_irf <- function(df, title) {
  ggplot(df, aes(x = step, y = irf)) +
    geom_ribbon(aes(ymin = lower95, ymax = upper95),
                fill = "red", colour = NA, alpha = 0.2) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Horizon (months)", y = NULL, title = title) +
    scale_x_continuous(breaks = seq(min(df$step), max(df$step), by = 1)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
}

save_irf <- function(plot, suffix) {
  ggsave(
    filename = file.path(fig_dir, sprintf("irf_opr_%s_%s.pdf", suffix, sample_tag)),
    plot     = plot,
    width    = 8,
    height   = 6
  )
}

# Yields
y2_df  <- irf_df %>% filter(impulse == "opr", response == "y2")
y5_df  <- irf_df %>% filter(impulse == "opr", response == "y5")
y10_df <- irf_df %>% filter(impulse == "opr", response == "y10")

p_y2  <- plot_irf(y2_df,  "Response of 2-year yield to a monetary policy shock")
p_y5  <- plot_irf(y5_df,  "Response of 5-year yield to a monetary policy shock")
p_y10 <- plot_irf(y10_df, "Response of 10-year yield to a monetary policy shock")

save_irf(p_y2,  "y2")
save_irf(p_y5,  "y5")
save_irf(p_y10, "y10")

# ==============================
# 8) FEVD for yields (24-month)
# ==============================

fe <- fevd(v, n.ahead = 24)

extract_fevd <- function(varname) {
  M  <- fe[[varname]]
  df <- as.data.frame(M)
  df$h <- seq_len(nrow(df)) - 1
  df %>%
    pivot_longer(-h, names_to = "shock", values_to = "share") %>%
    mutate(var = varname)
}

fe_df <- bind_rows(
  extract_fevd("y2"),
  extract_fevd("y5"),
  extract_fevd("y10")
)

first_mat   <- fe[[1]]
shock_names <- colnames(first_mat)
if (is.null(shock_names)) shock_names <- paste0("shock.", seq_along(ord))

fe_df$shock <- factor(
  fe_df$shock,
  levels = shock_names,
  labels = ord[seq_len(min(length(ord), length(shock_names)))]
)

p_fevd <- ggplot(fe_df, aes(x = h, y = share, fill = shock)) +
  geom_col(position = "fill", width = 0.9) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = seq(0, 12, 3)) +
  labs(
    title = "Forecast Error Variance Decomposition (FEVD)",
    x = "Horizon (months)",
    y = "Share of forecast error variance",
    fill = "Shock"
  ) +
  facet_wrap(~var, ncol = 1, scales = "free_y") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

fevd_path <- file.path(fig_dir, sprintf("fig_fevd_yields_bars_%s.png", sample_tag))
ggsave(filename = fevd_path, plot = p_fevd, width = 8, height = 6, dpi = 300)
message("FEVD figure saved to: ", fevd_path)

write_xlsx(
  fe_df,
  path = file.path(tab_dir, sprintf("fevd_yields_results_%s.xlsx", sample_tag))
)

# ==============================
# 9) Check for price puzzle (pi, opr)
# ==============================

pi_df  <- irf_df %>% filter(impulse == "opr", response == "pi")
opr_df <- irf_df %>% filter(impulse == "opr", response == "opr")

p_pi  <- plot_irf(pi_df,  "Response of inflation to a monetary policy shock")
p_opr <- plot_irf(opr_df, "Response of opr to a monetary policy shock")

save_irf(p_pi,  "pi")
save_irf(p_opr, "opr")

# ==============================
# 10) Save VAR diagnostics
# ==============================

var_base_out <- file.path(tab_dir, sprintf("var_results_yields_%s.txt", sample_tag))
sink(var_base_out)
cat("====================================\n")
cat("VAR RESULTS (YIELDS) - sample:", sample_tag, "\n")
cat("====================================\n\n")

cat("Selected lag (AIC, BIC, HQ): ", var_lag, "\n\n")
cat("Information criteria for estimated VAR:\n")
cat(sprintf("  AIC: %.4f\n", aic_val_model))
cat(sprintf("  BIC: %.4f\n", bic_val_model))

cat("\n--- VAR Coefficients (summary) ---\n")
print(summary(v))

cat("\n--- Residual covariance matrix ---\n")
print(cov(resid(v)))

cat("\n--- Stability roots (modulus) ---\n")
print(mod)
if (any(mod >= 1)) cat("\nWARNING: At least one root >= 1 → VAR not stable!\n")
sink()
message("Saved VAR results to: ", var_base_out)

# Serial correlation test and IRF object summary
serial.test(v)

cat("\n=== IRF object summary ===\n")
print(str(ir, max.level = 1))
