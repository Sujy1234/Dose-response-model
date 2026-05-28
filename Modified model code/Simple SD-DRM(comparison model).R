#load libraries
library(ggplot2)
library(fitdistrplus)
library(gridExtra)
library(egg)
library(grid)

#create output directories for final saved image and the result to be saved as csv
dir.create("results", showWarnings = FALSE)
dir.create("images", showWarnings = FALSE)

# Shared parameters for the antibiotic, borrowed from the dose-response paper
MIC <- 2
Emax <- 51 * 24
EC50 <- 9.93
doses <- 10^seq(2, 8, 0.5)

pMIC_values <- c(0, 0.5, 1, 2.5) / 100
fr_values <- c(0, 0.01, 0.05, 0.10)

fixed_fr <- 0.05
fixed_pMIC <- 1 / 100

# Dataset 1: exponential model, r value was borrowed from the paper
r_exp <- 1.07e-8
tfs_exp <- 1

# Dataset 2: beta-Poisson model, alpha and beta values borrowed from the paper again
alpha_bp <- 0.16
beta_bp <- 1.41e6
tfs_bp <- 2.625

# Safe calculation of log(1 - exp(-x))
log1mexpm <- function(a) {
  ans <- numeric(length(a))
  ans[a <= log(2)] <- log(-expm1(-a[a <= log(2)]))
  ans[a > log(2)] <- log1p(-exp(-a[a > log(2)]))
  ans
}

# Exponential model: convert r to r_s after antibiotic exposure
get_rs <- function(r, C, tfs) {
  mu <- -log1mexpm(r) / tfs
  mu_s <- mu + Emax * C / (EC50 + C)
  -log1mexpm(mu_s * tfs)
}

# Beta-Poisson model: simulate r, apply antibiotic effect, refit beta distribution
get_ab_beta <- function(C, seed = 0, nsim = 1000) {
  set.seed(seed)
  r <- rbeta(nsim, alpha_bp, beta_bp)
  mu <- -log1mexpm(r) / tfs_bp
  mu_s <- mu + Emax * C / (EC50 + C)
  r_s <- -log1mexpm(mu_s * tfs_bp)
  r_s <- pmin(pmax(r_s[is.finite(r_s)], 1e-16), 1 - 1e-16)
  fit <- fitdist(r_s, "beta")
  c(alpha_s = unname(fit$estimate["shape1"]), beta_s = unname(fit$estimate["shape2"]))
}

# Dataset 1: exponential SD-DRM
calc_exp <- function(fr, C) {
  r_s <- get_rs(r_exp, C, tfs_exp)
  Ns <- doses * (1 - fr)
  Nr <- doses * fr
  
  risk <- 1 - exp(-(r_s * Ns + r_exp * Nr))
  p_r <- 1 - exp(-r_exp * Nr)
  p_s <- exp(-r_exp * Nr) * (1 - exp(-r_s * Ns))
  
  data.frame(dose = log10(doses), risk = log10(pmax(risk, 1e-12)),
             outcome = ifelse(p_r >= p_s,   "Lower treatability", "Higher treatability"))
}

# Dataset 2: beta-Poisson SD-DRM
calc_bp <- function(fr, C, seed = 0) {
  pars <- get_ab_beta(C, seed)
  alpha_s <- pars["alpha_s"]
  beta_s <- pars["beta_s"]
  
  Ns <- doses * (1 - fr)
  Nr <- doses * fr
  
  risk <- 1 - (1 + Ns / beta_s)^(-alpha_s) * (1 + Nr / beta_bp)^(-alpha_bp)
  p_r <- 1 - (1 + Nr / beta_bp)^(-alpha_bp)
  p_s <- (1 + Nr / beta_bp)^(-alpha_bp) * (1 - (1 + Ns / beta_s)^(-alpha_s))
  
  data.frame(dose = log10(doses), risk = log10(pmax(risk, 1e-12)),
             outcome = ifelse(p_r >= p_s,   "Lower treatability", "Higher treatability"))
}

# Generate model outputs
make_data <- function(values, model, panel, type) {
  do.call(rbind, lapply(values, function(v) {
    if (type == "C") {
      out <- if (model == "exp") calc_exp(fixed_fr, v * MIC) else calc_bp(fixed_fr, v * MIC, 0)
      out$variable <- factor(v * 100, levels = c(0, 0.5, 1, 2.5))
      out$scenario <- "Varying %MIC"
    } else {
      out <- if (model == "exp") calc_exp(v, fixed_pMIC * MIC) else calc_bp(v, fixed_pMIC * MIC, 42)
      out$variable <- factor(v, levels = fr_values)
      out$scenario <- "Varying resistant fraction"
    }
    
    out$model <- ifelse(model == "exp", "Exponential", "Beta-Poisson")
    out$dataset <- ifelse(model == "exp", "Dataset 1", "Dataset 2")
    out$panel <- panel
    out
  }))
}

df <- rbind(
  make_data(pMIC_values, "exp", "A", "C"),
  make_data(fr_values, "exp", "B", "fr"),
  make_data(pMIC_values, "bp", "C", "C"),
  make_data(fr_values, "bp", "D", "fr")
)

write.csv(df, file.path("results", "SD_DRM_combined_results.csv"), row.names = FALSE)

# Plot settings (similar style and type like the original paper)
sz <- 11
sz_ax <- 14
sz_leg <- 14
sz_lab <- 6
llim <- -7.8

# Reusable paper-style plotting function
make_plot <- function(data, label, xlab, ylab, legend_title, show_x, show_y, show_legend) {
  ggplot(data, aes(dose)) +
    geom_line(aes(y = risk, colour = variable), linewidth = 1.5) +
    geom_point(aes(y = risk, colour = variable, shape = outcome), size = 4, alpha = 0.9) +
    xlab(xlab) + ylab(ylab) +
    labs(colour = legend_title, shape = "AB") +
    theme_bw(base_size = sz) +
    theme(
      legend.position = ifelse(show_legend, "top", "none"),
      legend.box = "vertical",
      axis.line = element_line(linewidth = 1),
      axis.text = element_text(face = "bold", size = sz_ax),
      axis.title = element_text(size = sz_ax),
      axis.text.x = if (show_x) element_text(face = "bold", size = sz_ax) else element_blank(),
      axis.text.y = if (show_y) element_text(face = "bold", size = sz_ax) else element_blank(),
      legend.text = element_text(size = sz_leg),
      legend.title = element_text(size = sz_leg),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA)
    ) +
    guides(shape = guide_legend(direction = "horizontal", order = 2),
           colour = guide_legend(direction = "horizontal", order = 1)) +
    scale_y_continuous(limits = c(llim, 0.5)) +
    annotate("text", x = 2, y = 0, label = label, size = sz_lab, fontface = "bold")
}

plot_specs <- list(
  list("A", "a", "", expression(Log[10]~(P[illness])), expression('%'~MIC), FALSE, TRUE, TRUE),
  list("B", "b", "", "", expression(f[r]), FALSE, FALSE, TRUE),
  list("C", "c", expression(Log[10]~(dose)), expression(Log[10]~(P[illness])), expression('%'~MIC), TRUE, TRUE, FALSE),
  list("D", "d", expression(Log[10]~(dose)), "", expression(f[r]), TRUE, FALSE, FALSE)
)

plots <- lapply(plot_specs, function(s) {
  make_plot(subset(df, panel == s[[1]]), s[[2]], s[[3]], s[[4]], s[[5]], s[[6]], s[[7]], s[[8]])
})

final_plot <- arrangeGrob(
  grobs = lapply(plots, set_panel_size, width = unit(10, "cm"), height = unit(9.5, "cm")),
  ncol = 2
)

grid.newpage()
grid.draw(final_plot)

ggsave(file.path("images", "SD_DRM_four_panel_paper_style.png"),
       final_plot, width = 23.7, height = 25, units = "cm", dpi = 300)


# Sensitivity analysis

library(fitdistrplus)
library(foreach)
library(doParallel)

dir.create("results", showWarnings = FALSE)
dir.create("images", showWarnings = FALSE)
run_pawn <- !file.exists(file.path("results", "sens.Rda"))

# Helper functions

log1mexpm <- function(a) {
  a0 <- log(2); ans <- numeric(length(a))
  ans[a <= a0] <- log(-expm1(-a[a <= a0]))
  ans[a > a0] <- log1p(-exp(-a[a > a0]))
  ans
}

getdat <- function(choice) {
  if (choice == 1) t <- 1
  if (choice == 2) t <- 2.625
  list(t = t, Emax = 1224, Ec50 = 9.93)
}

getr <- function(r, Clist, all_dat) {
  t <- max(all_dat$t); mu <- -1 / t * log1mexpm(r); r_s <- seq(length(Clist))
  for (ind in seq(length(Clist))) {
    mu_s <- mu + all_dat$Emax * Clist[ind] / (all_dat$Ec50 + Clist[ind])
    r_s[ind] <- -log1mexpm(mu_s * t)
  }
  r_s
}

getalphabeta <- function(alpha, beta, Clist, all_dat, seed) {
  set.seed(seed)
  t <- max(all_dat$t); r <- rbeta(1000, alpha, beta); mu <- -1 / t * log1mexpm(r)
  alpha_l <- seq(length(Clist)); beta_l <- seq(length(Clist))
  for (ind in seq(length(Clist))) {
    mu_s <- mu + all_dat$Emax * Clist[ind] / (all_dat$Ec50 + Clist[ind])
    r_s <- -log1mexpm(mu_s * t)
    tmpfit <- fitdist(r_s, "beta")
    alpha_l[ind] <- tmpfit[1]$estimate[1]; beta_l[ind] <- tmpfit[1]$estimate[2]
  }
  list(alpha_l, beta_l)
}

# Risk functions

get_exprisk <- function(x, p) {
  C <- x[1]; fr <- x[2]; dose <- 10^x[3]
  p$Emax <- x[4]; p$Ec50 <- x[5]; r <- x[6]; p$t <- x[7]
  r_s <- getr(r, C, p); Ns <- dose * (1 - fr); Nr <- dose * fr
  1 - exp(-(r_s * Ns + r * Nr))
}

get_bprisk <- function(x, p) {
  C <- x[1]; fr <- x[2]; dose <- 10^x[3]
  p$Emax <- x[4]; p$Ec50 <- x[5]; alpha <- x[6]; beta <- x[7]; p$t <- x[8]
  ab_pars <- getalphabeta(alpha, beta, C, p, 0)
  alpha_s <- ab_pars[[1]]; beta_s <- ab_pars[[2]]
  Ns <- dose * (1 - fr); Nr <- dose * fr
  1 - ((1 + Ns / beta_s)^(-alpha_s)) * ((1 + Nr / beta)^(-alpha))
}

get_exprisk_nodose <- function(x, p) {
  C <- x[1]; fr <- x[2]
  p$Emax <- x[3]; p$Ec50 <- x[4]; r <- x[5]; p$t <- x[6]
  r_s <- getr(r, C, p); Ns <- p$dose * (1 - fr); Nr <- p$dose * fr
  1 - exp(-(r_s * Ns + r * Nr))
}

get_bprisk_nodose <- function(x, p) {
  C <- x[1]; fr <- x[2]
  p$Emax <- x[3]; p$Ec50 <- x[4]; alpha <- x[5]; beta <- x[6]; p$t <- x[7]
  ab_pars <- getalphabeta(alpha, beta, C, p, 0)
  alpha_s <- ab_pars[[1]]; beta_s <- ab_pars[[2]]
  Ns <- p$dose * (1 - fr); Nr <- p$dose * fr
  1 - ((1 + Ns / beta_s)^(-alpha_s)) * ((1 + Nr / beta)^(-alpha))
}

# PAWN sensitivity analysis function

PAWN <- function(model, p, lb, ub, Nu, n, Nc, npts, seed, ncores = 4) {
  M <- length(lb)
  y_u <- rep(NaN, Nu); KS <- matrix(NaN, M, n); xvals <- matrix(NaN, M, n)
  par_u <- matrix(NaN, Nu, M); par_c <- matrix(NaN, Nc * n * M, M)
  y_c <- rep(NaN, M * Nc * n); ft <- matrix(NaN, M * n, npts)
  
  set.seed(seed)
  for (ind1 in seq(M)) {
    par_u[, ind1] <- runif(Nu, lb[ind1], ub[ind1])
    par_c[, ind1] <- runif(Nc * n * M, lb[ind1], ub[ind1])
  }
  for (ind1 in seq(M)) {
    for (ind2 in seq(n)) {
      left_index <- (ind1 - 1) * Nc * n + (ind2 - 1) * Nc + 1
      right_index <- (ind1 - 1) * Nc * n + ind2 * Nc
      xvals[ind1, ind2] <- runif(1, lb[ind1], ub[ind1])
      par_c[left_index:right_index, ind1] <- xvals[ind1, ind2]
    }
  }
  
  cl <- makeCluster(ncores); registerDoParallel(cl)
  clusterExport(cl, c("model", "p", "par_u", "par_c", "getr", "getalphabeta",
                      "log1mexpm", "get_exprisk", "get_bprisk",
                      "get_exprisk_nodose", "get_bprisk_nodose"), envir = environment())
  clusterEvalQ(cl, library(fitdistrplus))
  
  outputlist <- foreach(ind = 1:Nu) %dopar% model(par_u[ind, ], p)
  for (ind in seq(Nu)) y_u[ind] <- outputlist[[ind]]
  
  outputlist_c <- foreach(ind = 1:nrow(par_c)) %dopar% model(par_c[ind, ], p)
  for (ind in 1:nrow(par_c)) y_c[ind] <- outputlist_c[[ind]]
  
  stopCluster(cl)
  
  m1 <- min(c(min(y_c), min(y_u))); m2 <- max(c(max(y_c), max(y_c)))
  f <- ecdf(y_u)(seq(m1, m2, length.out = npts))
  
  for (ind1 in seq(M)) {
    for (ind2 in seq(n)) {
      yt <- y_c[((ind1 - 1) * Nc * n + (ind2 - 1) * Nc + 1):((ind1 - 1) * Nc * n + ind2 * Nc)]
      ft[(ind1 - 1) * n + ind2, ] <- ecdf(yt)(seq(m1, m2, length.out = npts))
      KS[ind1, ind2] <- max(abs(ft[(ind1 - 1) * n + ind2, ] - f))
    }
  }
  
  list(KS = KS, xvals = xvals, ft = ft, par_u = par_u, par_c = par_c, y_u = y_u, y_c = y_c)
}

# Sensitivity settings and fitted parameters

n <- 15; Nu <- 100; Nc <- 100; npts <- 100; seed <- 0
crit_c <- c(1.22, 1.36, 1.48, 1.63, 1.73, 1.95)
critval <- crit_c[2] * sqrt((Nu + Nc) / (Nu * Nc))

r <- 1.066597e-08
alpha <- 1.623144e-01; beta <- 1.414959e+06
Emax <- 1224; Ec50 <- 9.93
doselb <- 1; doseub <- 4; factor <- 0.5; dval <- 1e4

# Run or load PAWN results

if (run_pawn) {
  
  # Exponential model with varying dose
  lb <- c(0.00, 0.0, doselb, Emax-factor*Emax, Ec50-factor*Ec50, r-factor*r, 1.5)
  ub <- c(0.05, 0.1, doseub, Emax+factor*Emax, Ec50+factor*Ec50, r+factor*r, 2.5)
  res_exp <- PAWN(get_exprisk, getdat(1), lb, ub, Nu, n, Nc, npts, seed, 4)
  
  # Beta-Poisson model with varying dose
  lb <- c(0.00, 0.0, doselb, Emax-factor*Emax, Ec50-factor*Ec50, alpha-factor*alpha, beta-factor*beta, 2)
  ub <- c(0.05, 0.1, doseub, Emax+factor*Emax, Ec50+factor*Ec50, alpha+factor*alpha, beta+factor*beta, 3)
  res_beta <- PAWN(get_bprisk, getdat(2), lb, ub, Nu, n, Nc, npts, seed, 4)
  
  # Exponential model with fixed dose
  lb <- c(0.00, 0.0, Emax-factor*Emax, Ec50-factor*Ec50, r-factor*r, 1.5)
  ub <- c(0.05, 0.1, Emax+factor*Emax, Ec50+factor*Ec50, r+factor*r, 2.5)
  p <- getdat(1); p$dose <- dval
  resnd_exp <- PAWN(get_exprisk_nodose, p, lb, ub, Nu, n, Nc, npts, seed, 4)
  
  # Beta-Poisson model with fixed dose
  lb <- c(0.00, 0.0, Emax-factor*Emax, Ec50-factor*Ec50, alpha-factor*alpha, beta-factor*beta, 2)
  ub <- c(0.05, 0.1, Emax+factor*Emax, Ec50+factor*Ec50, alpha+factor*alpha, beta+factor*beta, 3)
  p <- getdat(2); p$dose <- dval
  resnd_beta <- PAWN(get_bprisk_nodose, p, lb, ub, Nu, n, Nc, npts, seed, 4)
  
  save(res_exp, resnd_exp, res_beta, resnd_beta, n, critval, Nu, Nc,
       file = file.path("results", "sens.Rda"))
  
} else {
  load(file.path("results", "sens.Rda"))
}

# Parameter labels

parnames1 <- c("C", expression(f[r]), "d", expression(E[max]), "EC50", "r", expression(t[fs]))
parnames2 <- c("C", expression(f[r]), "d", expression(E[max]), "EC50", expression(alpha), expression(beta), expression(t[fs]))
parnames3 <- c("C", expression(f[r]), expression(E[max]), "EC50", "r", expression(t[fs]))
parnames4 <- c("C", expression(f[r]), expression(E[max]), "EC50", expression(alpha), expression(beta), expression(t[fs]))

# Plot sensitivity figure

plot_sensitivity <- function(save = FALSE) {
  if (save) png(file.path("images", "sensitivity.png"), width = 8, height = 8, units = "in", res = 300)
  
  lwidth <- 5; ht <- 0.5; wd <- 0.5; bwidth <- 1.5
  c1 <- "#1b9e7777"; c2 <- "#d95f0277"
  cex_plot <- 1.3; cex_txt <- 1.5; rot <- 30
  set.seed(0)
  
  m <- matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE)
  layout(mat = m, heights = c(ht, 1 - ht), widths = c(wd, 1 - wd))
  
  # Panel A: exponential model with varying dose
  par(mar = c(4, 4, 0.5, 0), lwd = bwidth)
  boxplot(t(res_exp$KS), ylim = c(0, 1), ylab = "PAWN index", outline = FALSE,
          cex.axis = cex_plot, cex.sub = cex_plot, cex.lab = cex_plot, xaxt = "n")
  axis(1, labels = FALSE)
  text(x = seq_along(parnames1), y = -0.1, srt = rot, adj = 1, xpd = TRUE, labels = parnames1, cex = cex_plot)
  points(jitter(rep(seq(length(parnames1)), n), factor = 1.5), as.vector(res_exp$KS), bg = c1, pch = 21)
  text(0.7, 0.95, expression(bold("a")), cex = cex_txt)
  lines(c(0, 10), c(critval, critval), lty = 2)
  
  # Panel B: beta-Poisson model with varying dose
  par(mar = c(4, 4, 0.5, 0.2), lwd = bwidth)
  boxplot(t(res_beta$KS), ylim = c(0, 1), ylab = "", outline = FALSE,
          cex.axis = cex_plot, cex.sub = cex_plot, cex.lab = cex_plot, xaxt = "n", xlab = "")
  axis(1, at = seq_along(parnames2), labels = FALSE)
  text(x = seq_along(parnames2), y = -0.1, srt = rot, adj = 1, xpd = TRUE, labels = parnames2, cex = cex_plot)
  points(jitter(rep(seq(length(parnames2)), n), factor = 1.5), as.vector(res_beta$KS), bg = c1, pch = 21)
  text(0.7, 0.95, expression(bold("b")), cex = cex_txt)
  lines(c(0, 10), c(critval, critval), lty = 2)
  
  # Panel C: exponential model with fixed dose
  par(mar = c(4, 4, 0.5, 0), lwd = bwidth)
  boxplot(t(resnd_exp$KS), ylim = c(0, 1), ylab = "PAWN index", outline = FALSE,
          cex.axis = cex_plot, cex.sub = cex_plot, cex.lab = cex_plot, xaxt = "n")
  axis(1, labels = FALSE)
  text(x = seq_along(parnames3), y = -0.1, srt = rot, adj = 1, xpd = TRUE, labels = parnames3, cex = cex_plot)
  points(jitter(rep(seq(length(parnames3)), n), factor = 1.5), as.vector(resnd_exp$KS), bg = c2, pch = 21)
  text(0.7, 0.95, expression(bold("c")), cex = cex_txt)
  lines(c(0, 10), c(critval, critval), lty = 2)
  
  # Panel D: beta-Poisson model with fixed dose
  par(mar = c(4, 4, 0.5, 0.2), lwd = bwidth)
  boxplot(t(resnd_beta$KS), ylim = c(0, 1), outline = FALSE,
          cex.axis = cex_plot, cex.sub = cex_plot, cex.lab = cex_plot, xaxt = "n")
  axis(1, labels = FALSE)
  text(x = seq_along(parnames4), y = -0.1, srt = rot, adj = 1, xpd = TRUE, labels = parnames4, cex = cex_plot)
  points(jitter(rep(seq(length(parnames4)), n), factor = 1.5), as.vector(resnd_beta$KS), bg = c2, pch = 21)
  text(0.7, 0.95, expression(bold("d")), cex = cex_txt)
  lines(c(0, 10), c(critval, critval), lty = 2)
  
  if (save) dev.off()
}

# Show and save plot

plot_sensitivity(save = FALSE)
plot_sensitivity(save = TRUE)