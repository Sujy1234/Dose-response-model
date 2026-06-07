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


#Sensitivity analysis
# Sensitivity analysis for Vibrio beta-Poisson SD-DRM

library(foreach) 
library(doParallel)

run_pawn_vibrio <- !file.exists(file.path("results", "sens_vibrio.Rda"))

# Risk functions

pawn_bp <- function(x, p) {
  C <- x[1]; fr <- x[2]; dose <- 10^x[3]; Emax_i <- x[4]; EC50_i <- x[5]; alpha_i <- x[6]; beta_i <- x[7]; tfs_i <- x[8]
  pars <- get_ab_beta(alpha_i, beta_i, C, Emax_i, EC50_i, tfs_i, seed = 0)
  Ns <- dose * (1 - fr); Nr <- dose * fr
  1 - (1 + Ns / pars["beta_s"])^(-pars["alpha_s"]) * (1 + Nr / beta_i)^(-alpha_i)
}

pawn_bp_nodose <- function(x, p) {
  C <- x[1]; fr <- x[2]; Emax_i <- x[3]; EC50_i <- x[4]; alpha_i <- x[5]; beta_i <- x[6]; tfs_i <- x[7]
  pars <- get_ab_beta(alpha_i, beta_i, C, Emax_i, EC50_i, tfs_i, seed = 0)
  Ns <- p$dose * (1 - fr); Nr <- p$dose * fr
  1 - (1 + Ns / pars["beta_s"])^(-pars["alpha_s"]) * (1 + Nr / beta_i)^(-alpha_i)
}

# PAWN function

PAWN <- function(model, p, lb, ub, Nu, n, Nc, npts, seed, ncores = 4) {
  M <- length(lb); y_u <- rep(NaN, Nu); KS <- matrix(NaN, M, n); xvals <- matrix(NaN, M, n)
  par_u <- matrix(NaN, Nu, M); par_c <- matrix(NaN, Nc * n * M, M); y_c <- rep(NaN, M * Nc * n); ft <- matrix(NaN, M * n, npts)
  
  set.seed(seed)
  for (ind1 in seq(M)) { par_u[, ind1] <- runif(Nu, lb[ind1], ub[ind1]); par_c[, ind1] <- runif(Nc * n * M, lb[ind1], ub[ind1]) }
  for (ind1 in seq(M)) for (ind2 in seq(n)) {
    left_index <- (ind1 - 1) * Nc * n + (ind2 - 1) * Nc + 1; right_index <- (ind1 - 1) * Nc * n + ind2 * Nc
    xvals[ind1, ind2] <- runif(1, lb[ind1], ub[ind1]); par_c[left_index:right_index, ind1] <- xvals[ind1, ind2]
  }
  
  cl <- makeCluster(ncores); registerDoParallel(cl)
  clusterExport(cl, c("model", "p", "par_u", "par_c", "get_ab_beta", "log1mexpm"), envir = environment())
  clusterEvalQ(cl, library(fitdistrplus))
  outputlist <- foreach(ind = 1:Nu) %dopar% model(par_u[ind, ], p); for (ind in seq(Nu)) y_u[ind] <- outputlist[[ind]]
  outputlist_c <- foreach(ind = 1:nrow(par_c)) %dopar% model(par_c[ind, ], p); for (ind in 1:nrow(par_c)) y_c[ind] <- outputlist_c[[ind]]
  stopCluster(cl)
  
  m1 <- min(c(min(y_c), min(y_u))); m2 <- max(c(max(y_c), max(y_c))); f <- ecdf(y_u)(seq(m1, m2, length.out = npts))
  for (ind1 in seq(M)) for (ind2 in seq(n)) {
    yt <- y_c[((ind1 - 1) * Nc * n + (ind2 - 1) * Nc + 1):((ind1 - 1) * Nc * n + ind2 * Nc)]
    ft[(ind1 - 1) * n + ind2, ] <- ecdf(yt)(seq(m1, m2, length.out = npts))
    KS[ind1, ind2] <- max(abs(ft[(ind1 - 1) * n + ind2, ] - f))
  }
  list(KS = KS, xvals = xvals, ft = ft, par_u = par_u, par_c = par_c, y_u = y_u, y_c = y_c)
}

# Sensitivity settings

n <- 15; Nu <- 100; Nc <- 100; npts <- 100; seed <- 0
crit_c <- c(1.22, 1.36, 1.48, 1.63, 1.73, 1.95); critval <- crit_c[2] * sqrt((Nu + Nc) / (Nu * Nc))
factor <- 0.5; doselb <- 1; doseub <- 9; dval <- 1e6

# Run or load PAWN results

if (run_pawn_vibrio) {
  ab1 <- antibiotics[antibiotics$antibiotic == "Cefepime", ]; ab2 <- antibiotics[antibiotics$antibiotic == "Doxycycline", ]
  
  lb <- c(0, 0, doselb, ab1$Emax-factor*ab1$Emax, ab1$EC50-factor*ab1$EC50, alpha-factor*alpha, beta-factor*beta, tfs*0.5)
  ub <- c(0.025*ab1$MIC, 0.10, doseub, ab1$Emax+factor*ab1$Emax, ab1$EC50+factor*ab1$EC50, alpha+factor*alpha, beta+factor*beta, tfs*1.5)
  res_cef_dose <- PAWN(pawn_bp, list(), lb, ub, Nu, n, Nc, npts, seed, 4)
  
  lb <- c(0, 0, doselb, ab2$Emax-factor*ab2$Emax, ab2$EC50-factor*ab2$EC50, alpha-factor*alpha, beta-factor*beta, tfs*0.5)
  ub <- c(0.025*ab2$MIC, 0.10, doseub, ab2$Emax+factor*ab2$Emax, ab2$EC50+factor*ab2$EC50, alpha+factor*alpha, beta+factor*beta, tfs*1.5)
  res_dox_dose <- PAWN(pawn_bp, list(), lb, ub, Nu, n, Nc, npts, seed, 4)
  
  p <- list(dose = dval)
  lb <- c(0, 0, ab1$Emax-factor*ab1$Emax, ab1$EC50-factor*ab1$EC50, alpha-factor*alpha, beta-factor*beta, tfs*0.5)
  ub <- c(0.025*ab1$MIC, 0.10, ab1$Emax+factor*ab1$Emax, ab1$EC50+factor*ab1$EC50, alpha+factor*alpha, beta+factor*beta, tfs*1.5)
  res_cef_fixed <- PAWN(pawn_bp_nodose, p, lb, ub, Nu, n, Nc, npts, seed, 4)
  
  lb <- c(0, 0, ab2$Emax-factor*ab2$Emax, ab2$EC50-factor*ab2$EC50, alpha-factor*alpha, beta-factor*beta, tfs*0.5)
  ub <- c(0.025*ab2$MIC, 0.10, ab2$Emax+factor*ab2$Emax, ab2$EC50+factor*ab2$EC50, alpha+factor*alpha, beta+factor*beta, tfs*1.5)
  res_dox_fixed <- PAWN(pawn_bp_nodose, p, lb, ub, Nu, n, Nc, npts, seed, 4)
  
  save(res_cef_dose, res_dox_dose, res_cef_fixed, res_dox_fixed, n, critval, Nu, Nc, file = file.path("results", "sens_vibrio.Rda"))
} else load(file.path("results", "sens_vibrio.Rda"))

# Parameter labels

parnames_dose <- c("C", expression(f[r]), "d", expression(E[max]), "EC50", expression(alpha), expression(beta), expression(t[fs]))
parnames_fixed <- c("C", expression(f[r]), expression(E[max]), "EC50", expression(alpha), expression(beta), expression(t[fs]))

# Plot function

# Plot function

plot_ab_sensitivity <- function(res_dose, res_fixed, ab_name, file_name, save = FALSE) {
  if (save) png(file.path("images", file_name), width = 8, height = 4, units = "in", res = 300)
  
  c1 <- "#1b9e7777"; c2 <- "#d95f0277"
  cex_axis <- 0.9; cex_lab <- 1.0; cex_txt <- 1.0; cex_title <- 1.2
  rot <- 30; bwidth <- 1.5; set.seed(0)
  
  layout(matrix(c(1, 2), ncol = 2)); par(oma = c(0, 0, 2, 0))
  
  par(mar = c(4, 4, 2, 1), lwd = bwidth)
  boxplot(t(res_dose$KS), ylim = c(0, 1), ylab = "PAWN index", outline = FALSE,
          xaxt = "n", cex.axis = cex_axis, cex.lab = cex_lab)
  axis(1, at = seq_along(parnames_dose), labels = FALSE, tck = -0.02)
  text(seq_along(parnames_dose), -0.1, labels = parnames_dose, srt = rot, adj = 1, xpd = TRUE, cex = cex_axis)
  points(jitter(rep(seq(length(parnames_dose)), n), factor = 1.5), as.vector(res_dose$KS), bg = c1, pch = 21)
  text(0.7, 0.95, expression(bold("a")~"Varying dose"), adj = 0, cex = cex_txt)
  lines(c(0, 10), c(critval, critval), lty = 2)
  
  par(mar = c(4, 4, 2, 1), lwd = bwidth)
  boxplot(t(res_fixed$KS), ylim = c(0, 1), ylab = "", outline = FALSE,
          xaxt = "n", cex.axis = cex_axis, cex.lab = cex_lab)
  axis(1, at = seq_along(parnames_fixed), labels = FALSE, tck = -0.02)
  text(seq_along(parnames_fixed), -0.1, labels = parnames_fixed, srt = rot, adj = 1, xpd = TRUE, cex = cex_axis)
  points(jitter(rep(seq(length(parnames_fixed)), n), factor = 1.5), as.vector(res_fixed$KS), bg = c2, pch = 21)
  text(0.7, 0.95, expression(bold("b")~"Fixed dose"), adj = 0, cex = cex_txt)
  lines(c(0, 10), c(critval, critval), lty = 2)
  
  mtext(ab_name, outer = TRUE, font = 2, cex = cex_title)
  if (save) dev.off()
}

# Show and save plots

plot_ab_sensitivity(res_cef_dose, res_cef_fixed, "Cefepime", "Cefepime_sensitivity.png", save = FALSE)
plot_ab_sensitivity(res_dox_dose, res_dox_fixed, "Doxycycline", "Doxycycline_sensitivity.png", save = FALSE)

plot_ab_sensitivity(res_cef_dose, res_cef_fixed, "Cefepime", "Cefepime_sensitivity.png", save = TRUE)
plot_ab_sensitivity(res_dox_dose, res_dox_fixed, "Doxycycline", "Doxycycline_sensitivity.png", save = TRUE)