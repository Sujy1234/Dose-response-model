library(foreach)
library(doParallel)

PAWN <- function (model, p, lb, ub, Nu, n, Nc, npts, seed, ncores){
  #PAWN Run PAWN for Global Sensitivity Analysis of a supplied model
  #   [KS,xvals,y_u, y_c, par_u, par_c, ft] = PAWN(model, p, lb, ub, ...
  # Nu, n, Nc, npts, seed)
  #
  #     model : A function that takes a vector of parameters x and a 
  #             structure p as input to provide the model output as a scalar.
  #             The vector holds parameters we need sensitivity of and p
  #             holds all other parameters required to run the model.
  #     lb : A vector (1xM) of the lower bounds for each parameter
  #     ub : A vector (1xM) of the upper bounds for each parameter
  #     Nu : Number of samples of the unconditioned parameter space
  #     n : Number of conditioning values in each dimension
  #     Nc : Number of samples of the conditioned parameter space
  #     npts : Number of points to use in kernel density estimation
  #     seed : Random number seed
  #
  # ## Refernces
  # [1]: Pianosi, F., Wagener, T., 2015. A simple and efficient method for 
  # global sensitivity analysis based on cumulative distribution functions. 
  # Environ. Model. Softw. 67, 1?11. doi:10.1016/j.envsoft.2015.01.004
  
  
  # Initializations
  M = length(lb); # Number of parameters
  y_u = rep(NaN, Nu); # Ouput of unconditioned simulations
  y_c = n * Nc * M; # Output of conditioned simulations
  KS = rep(NaN, M*n); dim(KS) = c(M,n)# Kolmogorov-Smirnov statistic
  xvals = rep(NaN, M*n); dim(xvals) = c(M,n)# Container for conditioned samples
  ft = rep(NaN, M*n*npts); dim(ft) = c(M*n, npts)# CDF container
  
  set.seed(seed); # Set random seed
  # Containers for parameters
  # par_u = bsxfun(@plus, lb, bsxfun(@times, rand(Nu, M), (ub-lb)));
  par_u = rep(NaN, Nu*M); dim(par_u) = c(Nu, M);
  par_c = rep(NaN, M*M*Nc*n); dim(par_c) = c(Nc*n*M, M);
  xvals = rep(NaN, n*M); dim(xvals) = c(M, n);
  y_u = rep(NaN, Nu); y_c = rep(NaN, M*Nc*n);
  ft = rep(NaN, M*n*npts); dim(ft) = c(M*n, npts);
  for (ind1 in seq(M)){
    par_u[,ind1] = runif(Nu, lb[ind1], ub[ind1]);
    par_c[,ind1] = runif(Nc*n*M, lb[ind1], ub[ind1]);
  }
  for (ind1 in seq(M)){
    for (ind2 in seq(n)){
      left_index = (ind1-1)*Nc*n + (ind2-1)*Nc + 1
      right_index = (ind1-1)*Nc*n + (ind2)*Nc
      xvals[ind1, ind2] = runif(1, lb[ind1], ub[ind1])
      # print(c(left_index, right_index))
      par_c[left_index:right_index,ind1] = xvals[ind1, ind2];
    }
  }
  
  # Start parallel model evaluation
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  parallel::clusterExport(cl, "getr")
  parallel::clusterExport(cl, "getalphabeta")
  parallel::clusterExport(cl, "log1mexpm")
  # Evaluate model output of unconditioned samples
  outputlist <- foreach(ind=1:Nu) %dopar% {
    model(par_u[ind,], p);
  }
  
  for (ind in seq(Nu)) {
    y_u[ind] = outputlist[[ind]]
  }
  
  # Evaluate model output of conditioned samples
  outputlist_c <- foreach(ind=1:nrow(par_c)) %dopar% {
    model(par_c[ind,], p);
  }
  for (ind in 1:nrow(par_c)){
    y_c[ind] = outputlist_c[[ind]]
  }
  
  # Stop cluster
  stopCluster(cl)
  
  # Find bounds of the model outputs
  m1 = min(c(min(y_c), min(y_u)));
  m2 = max(c(max(y_c), max(y_c)));
  
  # Evaluate the CDF with kernel density for unconditioned samples
  emp_cdf = ecdf(y_u)
  f = emp_cdf(seq(m1, m2, length.out = npts))
  
  # Evaluate the CDF with kernel density for conditioned samples and use
  # that to find the KS statistic (Eqn 4 in the paper). 
  for (ind1 in seq(M)){
    for (ind2 in seq(n)){
      # Temporarily store the current conditioned samples
      yt = y_c[((ind1-1)*Nc*n+(ind2-1)*Nc+1):((ind1-1)*Nc*n+(ind2*Nc))];
      emp_cdf = ecdf(yt)
      ft[(ind1-1)*n+ind2,] = emp_cdf(seq(m1, m2, length.out = npts))
      KS[ind1,ind2] = max(abs(ft[(ind1-1)*n+ind2,]-f)); # Eqn 4
    }
  }
  res = list()
  res$KS = KS; res$xvals = xvals; res$ft = ft;
  res$par_u = par_u; res$par_c = par_c;
  res$y_u = y_u; res$y_c = y_c
  return (res)
}
