# Estimate pharmacodynamic parameters from 0–6 h time-kill data
# Sigmoid Emax model with bootstrap confidence intervals

# Working directory
setwd("/Users/gooner/Antibiotic Emax and EC50/Final emax model")

# Packages
library(dplyr)

# Inputs
dat <- read.csv("Antibiotic data with control.csv", check.names = FALSE)
time_window <- c(Cefepime = 6, Doxycycline = 6)

n_boot <- 1000
set.seed(123)

all_boot <- list()

# Function
get_pd <- function(abx) {
  
  # Select early time-kill data
  early <- filter(dat, Antibiotic == abx, Time_h <= time_window[abx])
  
  # Estimate bacterial growth rate in control group
  E0_fit <- lm(logCFU ~ Time_h,
               data = filter(early, Concentration_mg_L == 0))
  
  E0 <- coef(E0_fit)[2]
  E0_r2 <- summary(E0_fit)$r.squared
  
  cat("\n", abx, "\n")
  cat("Control growth rate, E0 =",
      round(E0, 4),
      "log10 CFU/mL/h; R² =",
      round(E0_r2, 4), "\n")
  
  # Estimate concentration-specific net growth/kill rates
  rates <- early %>%
    filter(Concentration_mg_L > 0) %>%
    group_by(Concentration_mg_L) %>%
    summarise(
      E_C = coef(lm(logCFU ~ Time_h))[2],
      r_squared = summary(lm(logCFU ~ Time_h))$r.squared,
      .groups = "drop"
    )
  
  print(rates)
  
  # Fit sigmoid Emax model
  fit <- nls(
    E_C ~ E0 - (Emax * Concentration_mg_L^Hill) /
      (EC50^Hill + Concentration_mg_L^Hill),
    data = rates,
    start = list(
      Emax = E0 - min(rates$E_C),
      EC50 = 0.2,
      Hill = 1
    ),
    algorithm = "port",
    lower = c(Emax = 0, EC50 = 0.001, Hill = 0.1),
    upper = c(Emax = 5, EC50 = 10, Hill = 10)
  )
  
  p <- coef(fit)
  
  # Bootstrap confidence intervals
  boot <- replicate(n_boot, {
    
    d <- slice_sample(rates, n = nrow(rates), replace = TRUE)
    
    tryCatch({
      coef(
        nls(
          E_C ~ E0 - (Emax * Concentration_mg_L^Hill) /
            (EC50^Hill + Concentration_mg_L^Hill),
          data = d,
          start = as.list(p),
          algorithm = "port",
          lower = c(Emax = 0, EC50 = 0.001, Hill = 0.1),
          upper = c(Emax = 5, EC50 = 10, Hill = 10)
        )
      )
    }, error = function(e) c(Emax = NA, EC50 = NA, Hill = NA))
  })
  
  boot <- na.omit(as.data.frame(t(boot)))
  boot$Antibiotic <- abx
  
  all_boot[[abx]] <<- boot
  
  ci <- apply(
    boot[, c("Emax", "EC50", "Hill")],
    2,
    quantile,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  
  # Return final SD-DRM parameters
  data.frame(
    Antibiotic = abx,
    
    Emax_ln_per_day = p["Emax"] * 24 * log(10),
    Emax_LCL_ln_per_day = ci[1, "Emax"] * 24 * log(10),
    Emax_UCL_ln_per_day = ci[2, "Emax"] * 24 * log(10),
    
    EC50_mg_L = p["EC50"],
    EC50_LCL_mg_L = ci[1, "EC50"],
    EC50_UCL_mg_L = ci[2, "EC50"],
    
    Hill = p["Hill"],
    Hill_LCL = ci[1, "Hill"],
    Hill_UCL = ci[2, "Hill"],
    
    Successful_bootstrap_runs = nrow(boot)
  )
}

# Run analysis
results <- bind_rows(lapply(names(time_window), get_pd))
boot_results <- bind_rows(all_boot)

# Display results
print(results)

# Save results
write.csv(results, "PD_PK_estimates_with_CI.csv", row.names = FALSE)
write.csv(boot_results, "PD_PK_bootstrap_parameter_distributions.csv", row.names = FALSE)