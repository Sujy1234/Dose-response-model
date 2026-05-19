getdat = function(choice) {
  if (choice == 1) {
    # pooled 39, 40 (DuPont 1971)
    dose = c(1e4, 1e4, 1e6, 1e6, 1e8, 1e8)
    nill = c(0, 0, 0, 1, 5, 3)
    ntot = c(5, 5, 5, 9, 8, 5)
    t = 1;
    all_dat = list(dose = dose, nill = nill, ntot = ntot, t = t, N = length(ntot))
  }
  else if (choice == 2) {
    # pooled, 153, 157, 159, 214, 216, 217
    # 214, 216, 217 (Bieber, 1998) -> tmax = 2
    # (Levine 1978) -> tmax = 2.625
    # (Levine 1973)
    dose = c(1e6, 1e6, 1e8, 5e8, 2.5e9, 1e10, 1e10, 1e10, 1e10, 2e10, 2.3e10)
    nill = c(0, 1, 1, 3, 6, 9, 9, 3, 5, 2, 14)
    ntot = c(4, 5, 5, 5, 6, 10, 14, 5, 5, 2, 19)
    t = 2.625;
    all_dat = list(dose = dose, nill = nill, ntot = ntot, t = t, N = length(ntot))
  }
  all_dat$t[all_dat$t > 20] = 0;
  all_dat$upperlim = 12;
  all_dat$lowerlim = 0;
  Emax = 51 * 24.0;
  # per day
  Ec50 = 9.93;
  # mg/L
  all_dat$Emax = Emax;
  all_dat$Ec50 = Ec50;
  return(all_dat)
}

log1mexpm <- function(a) {
  # https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
  # A numerically accurate way to compute log(1-exp(-x))
  a0 = log(2)
  less_than_mask = a <= a0
  greater_than_mask = a > a0
  ans = seq(length(a))
  ans[less_than_mask] = log(-expm1(-a[less_than_mask]));
  ans[greater_than_mask] = log1p(-exp(-a[greater_than_mask]));
  return(ans)
}

exp_deviance <- function(x, parameter_list) {
  # Return the deviance of the exponential model.
  target = 0
  dose = parameter_list$dose;
  nill = parameter_list$nill;
  ntot = parameter_list$ntot;
  r = 10 ^ x;
  devi = 0;
  for (n in 1:parameter_list$N) {
    temp = 1 - exp(-r * dose[n]);
    if (nill[n] == 0) {
      devi = devi + -2 * ntot[n] * dose[n] * (-r);
    }
    else if (nill[n] == ntot[n]) {
      devi = devi + -2 * ntot[n] * log1mexpm(r * dose[n]);
    }
    else {
      devi = devi + -2 * (nill[n] * log(temp * ntot[n] / (nill[n])) + (ntot[n] - nill[n]) * log((1 - temp) * ntot[n] / (ntot[n] - nill[n])));
    }
  }
  return(devi)
}


betap_deviance <- function(x, parameter_list) {
  # Return the deviance of the beta Poisson model.
  target = 0
  dose = parameter_list$dose;
  nill = parameter_list$nill;
  ntot = parameter_list$ntot;
  alpha = 10 ^ x[1];
  beta = 10 ^ x[2];
  devi = 0;
  for (n in 1:parameter_list$N) {
    temp = 1 - (1 + dose[n] / beta) ^ (-alpha);
    if (nill[n] == 0) {
      devi = devi + -2 * (ntot[n] * log(1 - temp));
    }
    else if (nill[n] == ntot[n]) {
      devi = devi + -2 * (ntot[n] * log(temp));
    }
    else {
      devi = devi + -2 * (nill[n] * log(temp * ntot[n] / (nill[n])) + (ntot[n] - nill[n]) * log((1 - temp) * ntot[n] / (ntot[n] - nill[n])));
    }
  }
  return(devi)
}


getr <- function(r, Clist, all_dat) {
  # For a given r, list of antibiotic concentrations and list of t, calculate 
  # r_{s,ab}(C) for each C in the list of concentrations.
  t = max(all_dat$t) # Get tfs
  mu = -1 / t * log1mexpm(r);
  # Calculate death rate, mu
  # Calculate mu_s and then r_{s,ab}
  r_s = seq(length(Clist));
  for (ind in seq(length(Clist))) {
    # for a given concentration
    mu_s = mu + all_dat$Emax * Clist[ind] / (all_dat$Ec50 + Clist[ind]);
    # day^-1
    r_s[ind] = -log1mexpm(mu_s * t);
  }
  return(r_s)
}


getalphabeta <- function(alpha, beta, Clist, all_dat, seed) {
  # For a given alpha, beta, list of antibiotic concentrations, list of t and
  # random seed, calculate alpha_{s,ab}(C) and beta_{s,ab}(C) for each C in
  # the list of concentrations.
  require(fitdistrplus)
  set.seed(seed)
  t = max(all_dat$t) # Get tfs
  # Sample 1000 values of r from a beta distribution.
  r = rbeta(1000, alpha, beta);
  # r[r<1e-16] = 1e-16;
  # Calculate the mu for each r.
  mu = -1 / t * log1mexpm(r);
  r_s = seq(length(r))
  alpha_l = seq(length(Clist))
  beta_l = seq(length(Clist))
  # Calculate mu_s and then fit a beta distribution to compute 
  # alpha_{s,AB} and beta_{s,AB}.
  for (ind in seq(length(Clist))) {
    # for a given concentration
    mu_s = mu + all_dat$Emax * Clist[ind] / (all_dat$Ec50 + Clist[ind]);
    # day^-1
    r_s = -log1mexpm(mu_s * max(all_dat$t));
    tmpfit = fitdist(r_s, 'beta');
    alpha_l[ind] = tmpfit[1]$estimate[1];
    beta_l[ind] = tmpfit[1]$estimate[2];
  }
  return(list(alpha_l, beta_l))
}


get_exprisk <- function(x, p) {
  # Return the risk for exponential model. Takes arguments in a format 
  # convenient for sensitivity analysis with varying dose.
  C = x[1];
  fr = x[2];
  dose = 10 ^ x[3];
  p$Emax = x[4];
  p$Ec50 = x[5];
  r = x[6];
  p$t = x[7];
  r_s = getr(r, C, p);
  Ns = dose * (1 - fr);
  Nr = dose * fr;
  risk = 1 - exp(-(r_s * Ns + r * Nr));
  return(risk)
}


get_bprisk <- function(x, p) {
  # Return the risk for beta Poisson model. Takes arguments in a format 
  # convenient for sensitivity analysis with varying dose.
  C = x[1];
  fr = x[2];
  dose = 10 ^ x[3];
  p$Emax = x[4];
  p$Ec50 = x[5];
  alpha = x[6];
  beta = x[7];
  p$t = x[8];
  ab_pars = getalphabeta(alpha, beta, C, p, 0);
  alpha_s = ab_pars[[1]];
  beta_s = ab_pars[[2]]
  Ns = dose * (1 - fr);
  Nr = dose * fr;
  risk = 1 - ((1 + Ns / beta_s) ^ (-alpha_s)) * ((1 + Nr / beta) ^ (-alpha));
}


get_exprisk_nodose <- function(x, p) {
  # Return the risk for exponential model. Takes arguments in a format 
  # convenient for sensitivity analysis with fixed dose.
  C = x[1];
  fr = x[2];
  p$Emax = x[3];
  p$Ec50 = x[4];
  r = x[5];
  p$t = x[6];
  dose = p$dose;
  r_s = getr(r, C, p);
  Ns = dose * (1 - fr);
  Nr = dose * fr;
  risk = 1 - exp(-(r_s * Ns + r * Nr));
  return(risk)
}


get_bprisk_nodose <- function(x, p) {
  # Return the risk for beta Poisson model. Takes arguments in a format 
  # convenient for sensitivity analysis with fixed dose.
  C = x[1];
  fr = x[2];
  p$Emax = x[3];
  p$Ec50 = x[4];
  alpha = x[5];
  beta = x[6];
  p$t = x[7];
  dose = p$dose;
  ab_pars = getalphabeta(alpha, beta, C, p, 0);
  alpha_s = ab_pars[[1]];
  beta_s = ab_pars[[2]]
  Ns = dose * (1 - fr);
  Nr = dose * fr;
  risk = 1 - ((1 + Ns / beta_s) ^ (-alpha_s)) * ((1 + Nr / beta) ^ (-alpha));
  return(risk)
}


get_plot_dataframes <- function(dataset_no, res) {
  # For a given dataset number and the results dataframe, 
  # return the data frames used to make the effect figure.
  all_dat = getdat(dataset_no)
  # Extract best fitting parameters
  r = res$r[dataset_no]
  alpha = res$alpha[dataset_no]
  beta = res$beta[dataset_no]
  doses = 10 ^ seq(2, 8, 0.5)
  ndose = length(doses);
  MIC = 2;
  #mg/L or microgram/mL

  pMIC = c(0, 0.5, 1, 2.5) / 100
  nconc = length(pMIC);
  df_dose = seq(nconc * ndose);
  df_C = seq(nconc * ndose);
  df_abfail = seq(nconc * ndose);
  df_risk = seq(nconc * ndose);
  df_riskbp = seq(nconc * ndose);
  df_abfailbp = seq(nconc * ndose);
  fr = 0.05;

  for (ind2 in seq(nconc)) {
    # for a given concentration
    leftind = (ind2 - 1) * ndose + 1
    rightind = ind2 * ndose
    # print(c(ind1, ind2, leftind, rightind))
    df_C[leftind:rightind] = pMIC[ind2] * 100
    df_dose[leftind:rightind] = doses
    r_s = getr(r, pMIC[ind2] * MIC, all_dat)
    res = getalphabeta(alpha, beta, pMIC[ind2] * MIC, all_dat, 0);
    alpha_s = res[[1]];
    beta_s = res[[2]];
    Nr = doses * fr;
    Ns = doses * (1 - fr);
    df_risk[leftind:rightind] = 1 - exp(-(r_s * Ns + r * Nr));
    df_riskbp[leftind:rightind] = 1 - (1 + Nr / beta) ^ (-alpha) * (1 + Ns / beta_s) ^ (-alpha_s)
    df_abfail[leftind:rightind] = ifelse(1 - exp(-r * Nr) >= exp(-r * Nr) * (1 - exp(-r_s * Ns)), "Untreatable", "Treatable");
    df_abfailbp[leftind:rightind] = ifelse(1 - (1 + Nr / beta) ^ (-alpha) >= ((1 + Nr / beta) ^ (-alpha)) * (1 - (1 + Ns / beta_s) ^ (-alpha_s)), "Untreatable", "Treatable");
  }
  df1 = data.frame(C = factor(df_C), risk = log10(df_risk), dose = log10(df_dose), abfail = df_abfail,
                   riskbp = log10(df_riskbp), abfailbp = df_abfailbp)

  C = 1 / 100 * MIC;
  # mg/L
  set.seed(42)
  r_s = getr(r, C, all_dat)
  res = getalphabeta(alpha, beta, C, all_dat, seed = 42);
  res
  alpha_s = res[[1]];
  beta_s = res[[2]];

  frlist = c(0, 0.01, 0.05, 0.1);
  nfr = length(frlist);
  df_fr = seq(nfr * ndose);
  df_dose = seq(nfr * ndose)
  df_risk = seq(nfr * ndose);
  df_abfail = seq(nfr * ndose)
  df_riskbp = seq(nfr * ndose);
  df_abfailbp = seq(nfr * ndose)
  for (ind2 in seq(nfr)) {
    # for a given fr
    leftind = (ind2 - 1) * ndose + 1
    rightind = ind2 * ndose
    df_fr[leftind:rightind] = frlist[ind2]
    df_dose[leftind:rightind] = doses
    Nr = doses * frlist[ind2];
    Ns = doses * (1 - frlist[ind2]);
    df_risk[leftind:rightind] = 1 - exp(-(r_s * Ns + r * Nr));
    # 
    df_abfail[leftind:rightind] = ifelse(1 - exp(-r * Nr) >= exp(-r * Nr) * (1 - exp(-r_s * Ns)), "Untreatable", "Treatable");
    df_riskbp[leftind:rightind] = 1 - (1 + Nr / beta) ^ (-alpha) * (1 + Ns / beta_s) ^ (-alpha_s)
    df_abfailbp[leftind:rightind] = ifelse(1 - (1 + Nr / beta) ^ (-alpha) >= ((1 + Nr / beta) ^ (-alpha)) * (1 - (1 + Ns / beta_s) ^ (-alpha_s)), "Untreatable", "Treatable");
  }
  df2 = data.frame(fr = factor(df_fr), risk = log10(df_risk), dose = log10(df_dose), abfail = df_abfail,
                   riskbp = log10(df_riskbp), abfailbp = df_abfailbp)
  return(list(df1, df2))
}

get_beta_verification <- function(dataset_no, res) {
  # Return the elements necessary to make the plot to visually demonstrate
  # the accuracy of the beta estimation procedure.
  all_dat = getdat(dataset_no)
  r = res$r[dataset_no]
  alpha = res$alpha[dataset_no]
  beta = res$beta[dataset_no]
  dose = 10 ^ seq(log10(min(all_dat$dose)), log10(max(all_dat$dose)), 0.02)
  MIC = 2;
  #mg/L or microgram/mL
  dose = 10 ^ seq(2, 8, 0.5)
  # Clist = c(0, 0.01,0.02,0.05)
  pMIC = c(0, 0.5, 1, 2.5) / 100
  t = max(all_dat$t)
  r = rbeta(1000, alpha, beta);
  # r[r<1e-16] = 1e-16;
  mu = -1 / t * log1mexpm(r);


  ndose = length(dose);
  nconc = length(pMIC);
  alpha_s = seq(nconc);
  beta_s = seq(nconc);
  r_s_list = seq(nconc * length(r));
  dim(r_s_list) = c(nconc, length(r))
  mu_s_list = seq(nconc * length(r));
  dim(mu_s_list) = c(nconc, length(r))
  for (ind2 in seq(nconc)) {
    # for a given concentration
    C = pMIC[ind2] * MIC
    set.seed(0)

    mu_s = mu + all_dat$Emax * C / (all_dat$Ec50 + C);
    # day^-1
    r_s = -log1mexpm(mu_s * max(all_dat$t));
    tmpfit = fitdist(r_s, 'beta');

    # res = getalphabeta(alpha, beta, pMIC[ind2] * MIC, all_dat, 0);
    mu_s_list[ind2,] = mu_s
    r_s_list[ind2,] = r_s
    alpha_s[ind2] = tmpfit[1]$estimate[1];
    beta_s[ind2] = tmpfit[1]$estimate[2];
  }

  return(list(alpha_s = alpha_s, beta_s = beta_s, r = r, mu = mu,
              mu_s = mu_s_list, r_s = r_s_list, alpha = alpha, beta = beta,
              pMIC = pMIC, MIC = MIC))
}