rm(list = ls())

library(ggplot2)
library(fitdistrplus)

# Create output folders

dir.create("results", showWarnings = FALSE)
dir.create("images", showWarnings = FALSE)

# Baseline beta-Poisson parameters

alpha <- 0.60
beta <- 1.31e6
tfs <- 1

doses <- 10^seq(1, 9, 0.5)

pMIC_values <- c(0, 0.5, 1, 2.5) / 100
fr_values <- c(0, 0.01, 0.05, 0.10)

fixed_pMIC <- 2 / 100

# Antibiotic-specific parameters

antibiotics <- data.frame(
  antibiotic = c("Cefepime", "Doxycycline"),
  MIC = c(0.5, 0.5),
  Emax = c(58.97, 63.03),
  EC50 = c(0.1137, 0.2184 ),
  fixed_fr = c(0.09, 0.20)
)

# Safe calculation of log(1 - exp(-x))

log1mexpm <- function(a) {
  ans <- numeric(length(a))
  ans[a <= log(2)] <- log(-expm1(-a[a <= log(2)]))
  ans[a > log(2)] <- log1p(-exp(-a[a > log(2)]))
  ans
}

# Simulate r values, apply antibiotic effect, and refit beta distribution

get_ab_beta <- function(alpha, beta, C, Emax, EC50, tfs, seed = 0, nsim = 1000) {
  
  set.seed(seed)
  
  r <- rbeta(nsim, alpha, beta)
  
  mu <- -1 / tfs * log1mexpm(r)
  
  mu_s <- mu + Emax * C / (EC50 + C)
  
  r_s <- -log1mexpm(mu_s * tfs)
  
  r_s <- r_s[is.finite(r_s)]
  r_s <- pmin(pmax(r_s, 1e-16), 1 - 1e-16)
  
  fit <- fitdist(r_s, "beta")
  
  c(
    alpha_s = unname(fit$estimate["shape1"]),
    beta_s = unname(fit$estimate["shape2"])
  )
}

# Beta-Poisson SD-DRM model

calc_bp_sd_drm <- function(fr, C, Emax, EC50, seed = 0) {
  
  pars <- get_ab_beta(
    alpha = alpha,
    beta = beta,
    C = C,
    Emax = Emax,
    EC50 = EC50,
    tfs = tfs,
    seed = seed
  )
  
  alpha_s <- pars["alpha_s"]
  beta_s <- pars["beta_s"]
  
  Ns <- doses * (1 - fr)
  Nr <- doses * fr
  
  riskbp <- 1 -
    (1 + Ns / beta_s)^(-alpha_s) *
    (1 + Nr / beta)^(-alpha)
  
  p_r_only <- 1 - (1 + Nr / beta)^(-alpha)
  
  p_s_only <- ((1 + Nr / beta)^(-alpha)) *
    (1 - (1 + Ns / beta_s)^(-alpha_s))
  
  status <- ifelse(
    p_r_only >= p_s_only,
    "Less likely treatable",
    "More likely treatable")
  
  data.frame(
    dose = log10(doses),
    riskbp = riskbp,
    riskbp_log10 = log10(pmax(riskbp, 1e-12)),
    abfailbp = status,
    alpha_s = alpha_s,
    beta_s = beta_s
  )
}

# Scenario A: varying %MIC

results_C <- do.call(rbind, lapply(1:nrow(antibiotics), function(i) {
  
  ab <- antibiotics[i, ]
  
  do.call(rbind, lapply(pMIC_values, function(pm) {
    
    out <- calc_bp_sd_drm(
      fr = ab$fixed_fr,
      C = pm * ab$MIC,
      Emax = ab$Emax,
      EC50 = ab$EC50,
      seed = 0
    )
    
    out$antibiotic <- ab$antibiotic
    out$pMIC <- factor(pm * 100, levels = c(0, 0.5, 1, 2.5))
    out$fr <- ab$fixed_fr
    out$scenario <- "Varying %MIC"
    
    out
  }))
}))

# Scenario B: varying resistant fraction

results_fr <- do.call(rbind, lapply(1:nrow(antibiotics), function(i) {
  
  ab <- antibiotics[i, ]
  
  do.call(rbind, lapply(fr_values, function(fr) {
    
    out <- calc_bp_sd_drm(
      fr = fr,
      C = fixed_pMIC * ab$MIC,
      Emax = ab$Emax,
      EC50 = ab$EC50,
      seed = 42
    )
    
    out$antibiotic <- ab$antibiotic
    out$pMIC <- fixed_pMIC * 100
    out$fr <- factor(fr, levels = fr_values)
    out$scenario <- "Varying resistant fraction"
    
    out
  }))
}))

# Save CSV outputs

write.csv(
  results_C,
  file.path("results", "BP_SD_DRM_vary_MIC.csv"),
  row.names = FALSE
)

write.csv(
  results_fr,
  file.path("results", "BP_SD_DRM_vary_fr.csv"),
  row.names = FALSE
)

# Plot A: varying %MIC

plot_C <- ggplot(results_C, aes(x = dose, y = riskbp_log10)) +
  geom_line(aes(colour = pMIC), linewidth = 1.2) +
  geom_point(aes(colour = pMIC, shape = abfailbp), size = 2.2, alpha = 0.9) +
  facet_wrap(~ antibiotic) +
  coord_cartesian(ylim = c(-8, 0.5)) +
  labs(
    x = expression(Log[10]~dose),
    y = expression(Log[10]~(P[illness])),
    colour = expression('%'~MIC),
    shape = "Outcome",
    title = "Beta-Poisson SD-DRM: varying %MIC"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

print(plot_C)

ggsave(
  file.path("images", "BP_SD_DRM_vary_MIC.png"),
  plot_C,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)

# Plot B: varying resistant fraction

plot_fr <- ggplot(results_fr, aes(x = dose, y = riskbp_log10)) +
  geom_line(aes(colour = fr), linewidth = 1.2) +
  geom_point(aes(colour = fr, shape = abfailbp), size = 2.2, alpha = 0.9) +
  facet_wrap(~ antibiotic) +
  coord_cartesian(ylim = c(-8, 0.5)) +
  labs(
    x = expression(Log[10]~dose),
    y = expression(Log[10]~(P[illness])),
    colour = expression(f[r]),
    shape = "AB",
    title = "Beta-Poisson SD-DRM: varying resistant fraction"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

print(plot_fr)

ggsave(
  file.path("images", "BP_SD_DRM_vary_fr.png"),
  plot_fr,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)