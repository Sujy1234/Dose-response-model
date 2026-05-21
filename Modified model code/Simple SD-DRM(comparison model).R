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