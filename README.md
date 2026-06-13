## =========================================================
## ONE-SCRIPT RUN (copy-paste):
##   1) Schoenfeld-based design + Monte Carlo power simulation
##   2) Three scenario sets:
##        (A) Extreme HR
##        (B) Unequal allocation (allocation fixed at 0.2)
##        (C) Unequal dropout
##   3) Side-by-side, publication-quality plot (no patchwork, no cairo)
##      - facet strips use plotmath lambda[0] == 0.01 / 0.05 (no boxes)
##      - Extreme HR: x labels "HR = ...", no x-axis title
##      - Unequal allocation: title includes (allocation = 0.2), x labels only HR
##      - Unequal dropout: title includes control/treatment dropout once, x labels only HR
##   4) Save PDF/TIFF/PNG without XQuartz/Cairo
## =========================================================

## ---------------------------
## Packages
## ---------------------------
library(survival)
library(ggplot2)
library(dplyr)

if (!requireNamespace("gridExtra", quietly = TRUE)) {
  install.packages("gridExtra", repos = "https://cloud.r-project.org")
}
library(grid)
library(gridExtra)

## ---------------------------
## Core functions
## ---------------------------

required_events_schoenfeld <- function(HR, alpha, target_power, alloc) {
  z_alpha  <- qnorm(1 - alpha / 2)
  z_beta   <- qnorm(target_power)
  (z_alpha + z_beta)^2 / ((log(HR))^2 * alloc * (1 - alloc))
}

p_event_arm <- function(lambda, lambda_drop, T) {
  lambda * (1 - exp(-(lambda + lambda_drop) * T)) / (lambda + lambda_drop)
}

binom_ci_exact <- function(x, n, conf_level = 0.95) {
  bt <- binom.test(x, n, conf.level = conf_level)
  as.numeric(bt$conf.int)  # c(lower, upper)
}

simulate_trial_once <- function(N, alloc,
                                lambda0, HR,
                                lambda_drop0, lambda_drop1, T) {
  N_treated  <- round(N * alloc)
  N_control  <- N - N_treated
  arm        <- c(rep(0L, N_control), rep(1L, N_treated))
  
  lambda_vec <- ifelse(arm == 0L, lambda0, lambda0 * HR)
  
  t_event <- rexp(N, rate = lambda_vec)
  lambda_drop_vec <- ifelse(arm == 0L, lambda_drop0, lambda_drop1)
  t_drop  <- rexp(N, rate = lambda_drop_vec)
  
  t_obs  <- pmin(t_event, t_drop, T)
  status <- as.integer(t_event <= t_drop & t_event <= T)
  
  fit   <- survdiff(Surv(t_obs, status) ~ arm)
  pval  <- 1 - pchisq(fit$chisq, df = 1)
  
  list(events = sum(status), pval = pval)
}

run_scenario <- function(HR, lambda0, lambda_drop0, lambda_drop1, T,
                         alpha, target_power, alloc, nsim,
                         conf_level = 0.95) {
  
  D_req <- required_events_schoenfeld(
    HR = HR, alpha = alpha, target_power = target_power, alloc = alloc
  )
  
  p0 <- p_event_arm(lambda = lambda0,      lambda_drop = lambda_drop0, T = T)
  p1 <- p_event_arm(lambda = lambda0 * HR, lambda_drop = lambda_drop1, T = T)
  
  p_bar <- (1 - alloc) * p0 + alloc * p1
  N_req <- ceiling(D_req / p_bar)
  
  events_vec <- numeric(nsim)
  sig_vec    <- logical(nsim)
  
  for (b in seq_len(nsim)) {
    sim_b <- simulate_trial_once(
      N = N_req, alloc = alloc,
      lambda0 = lambda0, HR = HR,
      lambda_drop0 = lambda_drop0, lambda_drop1 = lambda_drop1,
      T = T
    )
    events_vec[b] <- sim_b$events
    sig_vec[b]    <- (sim_b$pval < alpha)
  }
  
  x_sig     <- sum(sig_vec)
  emp_power <- x_sig / nsim
  ci <- binom_ci_exact(x = x_sig, n = nsim, conf_level = conf_level)
  
  data.frame(
    HR            = HR,
    lambda0       = lambda0,
    lambda_drop0  = lambda_drop0,
    lambda_drop1  = lambda_drop1,
    T             = T,
    alpha         = alpha,
    target_power  = target_power,
    alloc         = alloc,
    nsim          = nsim,
    D_req         = D_req,
    N_req         = N_req,
    mean_events   = mean(events_vec),
    sd_events     = sd(events_vec),
    emp_power     = emp_power,
    power_ci_low  = ci[1],
    power_ci_high = ci[2],
    stringsAsFactors = FALSE
  )
}

run_scenario_grid <- function(scenarios,
                              alpha, target_power,
                              alloc, nsim,
                              conf_level = 0.95) {
  results_list <- vector("list", nrow(scenarios))
  for (i in seq_len(nrow(scenarios))) {
    s <- scenarios[i, ]
    res_i <- run_scenario(
      HR           = s$HR,
      lambda0      = s$lambda0,
      lambda_drop0 = s$lambda_drop0,
      lambda_drop1 = s$lambda_drop1,
      T            = s$T,
      alpha        = alpha,
      target_power = target_power,
      alloc        = alloc,
      nsim         = nsim,
      conf_level   = conf_level
    )
    res_i$scenario <- s$scenario
    results_list[[i]] <- res_i
  }
  do.call(rbind, results_list)
}

## =========================================================
## Simulation parameters
## =========================================================
seed <- 20251117
set.seed(seed)

alpha        <- 0.05
target_power <- 0.80
nsim         <- 2000
conf_level   <- 0.95
T_follow     <- 3

lambda0_vals <- c(0.01, 0.05)
drop0_fixed  <- 0.001

## =========================================================
## (A) Extreme HR
## =========================================================
alloc_extreme <- 0.50
HR_extreme    <- c(0.1, 2.5)
drop1_extreme <- 0.001

scenarios_extreme <- expand.grid(
  HR           = HR_extreme,
  lambda0      = lambda0_vals,
  lambda_drop0 = drop0_fixed,
  lambda_drop1 = drop1_extreme,
  T            = T_follow,
  stringsAsFactors = FALSE
) %>%
  mutate(scenario = paste0("HR_", HR, "_lam0_", lambda0)) %>%
  select(scenario, HR, lambda0, lambda_drop0, lambda_drop1, T)

res_extreme_hr <- run_scenario_grid(
  scenarios    = scenarios_extreme,
  alpha        = alpha,
  target_power = target_power,
  alloc        = alloc_extreme,
  nsim         = nsim,
  conf_level   = conf_level
)

## =========================================================
## (B) Unequal allocation (allocation fixed at 0.2)
## =========================================================
alloc_unequal <- 0.20
HR_alloc      <- c(0.3, 0.5, 0.7)
drop1_alloc   <- 0.001

scenarios_alloc <- expand.grid(
  HR           = HR_alloc,
  lambda0      = lambda0_vals,
  lambda_drop0 = drop0_fixed,
  lambda_drop1 = drop1_alloc,
  T            = T_follow,
  stringsAsFactors = FALSE
) %>%
  mutate(scenario = paste0("HR_", HR, "_lam0_", lambda0, "_allocation_", alloc_unequal)) %>%
  select(scenario, HR, lambda0, lambda_drop0, lambda_drop1, T)

res_unequal_alloc <- run_scenario_grid(
  scenarios    = scenarios_alloc,
  alpha        = alpha,
  target_power = target_power,
  alloc        = alloc_unequal,
  nsim         = nsim,
  conf_level   = conf_level
)

## =========================================================
## (C) Unequal dropout
## =========================================================
alloc_dropout <- 0.50
HR_dropout    <- c(0.3, 0.5, 0.7)
drop1_dropout <- 0.05

scenarios_dropout <- expand.grid(
  HR           = HR_dropout,
  lambda0      = lambda0_vals,
  lambda_drop0 = drop0_fixed,
  lambda_drop1 = drop1_dropout,
  T            = T_follow,
  stringsAsFactors = FALSE
) %>%
  mutate(scenario = paste0("HR_", HR, "_lam0_", lambda0, "_drop0_", lambda_drop0, "_drop1_", lambda_drop1)) %>%
  select(scenario, HR, lambda0, lambda_drop0, lambda_drop1, T)

res_unequal_dropout <- run_scenario_grid(
  scenarios    = scenarios_dropout,
  alpha        = alpha,
  target_power = target_power,
  alloc        = alloc_dropout,
  nsim         = nsim,
  conf_level   = conf_level
)

## Optional: print results tables
print(res_extreme_hr)
print(res_unequal_alloc)
print(res_unequal_dropout)

## =========================================================
## Plot styling (colorblind-safe; Nature-friendly)
## =========================================================
OKABE_ITO <- c(
  blue   = "#E05454",
  orange = "#C13383",
  green  = "#792CA2",
  grey   = "#443199"
)

theme_natureish <- theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0),
    axis.title      = element_text(face = "bold"),
    axis.text.x     = element_text(size = 9),
    axis.ticks      = element_line(linewidth = 0.4),
    axis.line       = element_line(linewidth = 0.5),
    strip.text      = element_text(face = "bold", size = 10),
    strip.background= element_blank(),
    legend.position = "none",
    plot.margin     = margin(5.5, 8, 5.5, 5.5)
  )

fmt_num <- function(x) format(as.numeric(x), scientific = FALSE, trim = TRUE)

# Facet labeler using plotmath: lambda[0] == 0.01 / 0.05
lab_lam0 <- as_labeller(function(x) {
  paste0("Baseline hazard = ", fmt_num(x))
})

## ---------------------------
## Panel A: Extreme HR
##   x labels: "HR = ..."; no x-axis title
## ---------------------------
dfA <- res_extreme_hr %>%
  mutate(
    HR_f  = factor(HR, levels = sort(unique(HR))),
    HR_lb = factor(paste0("HR = ", HR_f), levels = paste0("HR = ", levels(HR_f)))
  )

pA <- ggplot(dfA, aes(x = HR_lb, y = emp_power)) +
  geom_col(width = 0.72, fill = OKABE_ITO["blue"]) +
  geom_errorbar(aes(ymin = power_ci_low, ymax = power_ci_high),
                width = 0.18, linewidth = 0.5) +
  geom_hline(yintercept = target_power, linetype = "dashed",
             linewidth = 0.5, color = OKABE_ITO["grey"]) +
  facet_wrap(~ lambda0, nrow = 1, labeller = lab_lam0) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Extreme HR", x = NULL, y = "Empirical power") +
  theme_natureish

## ---------------------------
## Panel B: Unequal allocation
##   title includes allocation = 0.2
##   x labels: HR only
## ---------------------------
dfB <- res_unequal_alloc %>%
  mutate(
    HR_f  = factor(HR, levels = sort(unique(HR))),
    HR_lb = factor(paste0("HR = ", HR_f), levels = paste0("HR = ", levels(HR_f)))
  )

alloc_fixed <- unique(res_unequal_alloc$alloc)
alloc_title <- if (length(alloc_fixed) == 1) fmt_num(alloc_fixed) else "varied"

pB <- ggplot(dfB, aes(x = HR_lb, y = emp_power)) +
  geom_col(width = 0.72, fill = OKABE_ITO["orange"]) +
  geom_errorbar(aes(ymin = power_ci_low, ymax = power_ci_high),
                width = 0.18, linewidth = 0.5) +
  geom_hline(yintercept = target_power, linetype = "dashed",
             linewidth = 0.5, color = OKABE_ITO["grey"]) +
  facet_wrap(~ lambda0, nrow = 1, labeller = lab_lam0) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = paste0("Unequal allocation (allocation = ", alloc_title, ")"),
       x = NULL, y = NULL) +
  theme_natureish

## ---------------------------
## Panel C: Unequal dropout
##   title includes control/treatment dropout once
##   x labels: HR only
## ---------------------------
dfC <- res_unequal_dropout %>%
  mutate(
    HR_f  = factor(HR, levels = sort(unique(HR))),
    HR_lb = factor(paste0("HR = ", HR_f), levels = paste0("HR = ", levels(HR_f)))
  )

drop0_u <- unique(res_unequal_dropout$lambda_drop0)
drop1_u <- unique(res_unequal_dropout$lambda_drop1)

drop_title <- if (length(drop0_u) == 1 && length(drop1_u) == 1) {
  paste0("control = ", fmt_num(drop0_u), ", treatment = ", fmt_num(drop1_u))
} else {
  "control/treatment varied"
}

pC <- ggplot(dfC, aes(x = HR_lb, y = emp_power)) +
  geom_col(width = 0.72, fill = OKABE_ITO["green"]) +
  geom_errorbar(aes(ymin = power_ci_low, ymax = power_ci_high),
                width = 0.18, linewidth = 0.5) +
  geom_hline(yintercept = target_power, linetype = "dashed",
             linewidth = 0.5, color = OKABE_ITO["grey"]) +
  facet_wrap(~ lambda0, nrow = 1, labeller = lab_lam0) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = paste0("Unequal dropout (dropout hazards: ", drop_title, ")"),
       x = NULL, y = NULL) +
  theme_natureish

## =========================================================
## Combine side-by-side (gridExtra)
## =========================================================
title_grob <- textGrob(
  "Empirical power under stress-test design scenarios",
  gp = gpar(fontface = "bold", fontsize = 13),
  just = "left", x = 0
)
subtitle_grob <- textGrob(
  "Bars: Monte Carlo power; Error bars: 95% exact binomial CI; Dashed line: target power",
  gp = gpar(fontsize = 10),
  just = "left", x = 0
)
top_grob <- arrangeGrob(title_grob, subtitle_grob, ncol = 1,
                        heights = unit(c(1.2, 1.0), "lines"))

fig <- arrangeGrob(pA, pB, pC, nrow = 1, top = top_grob)

grid.newpage()
grid.draw(fig)

## =========================================================
## Save without Cairo/X11 (macOS-friendly)
## =========================================================

# Vector PDF
pdf("power_side_by_side.pdf", width = 13.5, height = 5.0, onefile = TRUE)
grid.newpage(); grid.draw(fig)
dev.off()

# High-DPI TIFF
tiff("power_side_by_side.tiff",
     width = 13.5, height = 5.0, units = "in",
     res = 600, compression = "lzw")
grid.newpage(); grid.draw(fig)
dev.off()

# Optional PNG (600 dpi)
png("power_side_by_side.png",
    width = 13.5, height = 5.0, units = "in",
    res = 600, type = "quartz")
grid.newpage(); grid.draw(fig)
dev.off()
