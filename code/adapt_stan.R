## Bayesian Regularized SEM Software: Analysis with Stan
## Author: Sara van Erp

library(rstan)
options(mc.cores = 3)
rstan_options(auto_write = TRUE)
library(shinystan)

set.seed(09052023)

source("./code/stan_conv_fun.R")

## Data: already standardized
load("./data/testDatSD_adapt.dat")

## Run Stan with ridge and regularized horseshoe prior
run.stan <- function(dat, prior, n.chains = 3, sample = 10000, algorithm = c("mcmc", "vb")){
  
  if(prior == "ridge"){
    standat <- list(N = nrow(testDatSD),
                    P = ncol(testDatSD),
                    Q = 3,
                    y = testDatSD,
                    nC = 72)
    stanmod <- stan_model("./code/stanmodels/adapt_ridge.stan")
    parsel <- c("psi", "L_main_C", "L_cross_C", "phi_C", "sd0")
  } else if(prior == "reghs_def"){
    standat <- list(N = nrow(testDatSD),
                    P = ncol(testDatSD),
                    Q = 3,
                    y = testDatSD,
                    nC = 72,
                    scale_global = 1,
                    nu_global = 1,
                    nu_local = 1,
                    slab_scale = 2.5,
                    slab_df = 3)
    stanmod <- stan_model("./code/stanmodels/adapt_reghs.stan")
    parsel <- c("psi", "L_main_C", "L_cross_C", "phi_C")
  } else {
    # Base global scale on prior guess, which results in a value much lower than 1
    p0 = 7 # estimated number of relevant cross-loadings (based on estimated cross-loadings > 0.20 in initial EFA)
    p = 73 # total number of cross loadings 
    scale_global = p0/(p-p0)*1/sqrt(nrow(testDatSD))
    
    standat <- list(N = nrow(testDatSD),
                    P = ncol(testDatSD),
                    Q = 3,
                    y = testDatSD,
                    nC = 72,
                    scale_global = scale_global,
                    nu_global = 1,
                    nu_local = 1,
                    slab_scale = 2.5,
                    slab_df = 3)
    stanmod <- stan_model("./code/stanmodels/adapt_reghs.stan")
    parsel <- c("psi", "L_main_C", "L_cross_C", "phi_C")
  }
  
  if(algorithm == "mcmc"){
    fit <- sampling(stanmod, data = standat, chains = n.chains, iter = sample, pars = parsel)
  } else {
    fit <- vb(stanmod, data = standat, pars = parsel)
  }
  
  save(fit, file = paste0("./results/fit_adapt_stan_", prior, "_", algorithm, ".RData"))
  
}

run.stan(testDatSD, prior = "ridge", algorithm = "mcmc")
run.stan(testDatSD, prior = "ridge", algorithm = "vb")

run.stan(testDatSD, prior = "reghs_def", algorithm = "mcmc")
run.stan(testDatSD, prior = "reghs_def", algorithm = "vb")

run.stan(testDatSD, prior = "reghs_p0", algorithm = "mcmc")
run.stan(testDatSD, prior = "reghs_p0", algorithm = "vb")

## Check convergence
cutoffs <- list("rhat" = 1.01, "n_eff" = 1000, "ess_bulk" = 1000, "ess_tail" = 1000, "div" = 0)

# full Bayesian models
stanls.mcmc <- list(NA)
prior <- c("ridge", "reghs_def", "reghs_p0")
for(i in 1:length(prior)){
  algorithm <- "mcmc"
  load(paste0("./results/fit_adapt_stan_", prior[i], "_", algorithm, ".RData"))
  stanls.mcmc[[i]] <- fit
}

names(stanls.mcmc) <- c("fit_ridge_mcmc", "fit_reghs_def_mcmc", "fit_reghs_p0_mcmc")

stan_conv_fun(stanls.mcmc, cutoffs, algorithm = "mcmc") # effective N is quite low especially for the variance of factor 3

# Approximate Bayesian models
stanls.vb <- list(NA)
for(i in 1:length(prior)){
  algorithm <- "vb"
  load(paste0("./results/fit_adapt_stan_", prior[i], "_", algorithm, ".RData"))
  stanls.vb[[i]] <- fit
}

names(stanls.vb) <- c("fit_ridge_vb", "fit_reghs_def_vb", "fit_reghs_p0_vb")

stan_conv_fun(stanls.vb, cutoffs, algorithm = "vb") 

## Extract results
prior <- c("ridge", "reghs_def", "reghs_p0")
algorithm <- c("mcmc", "vb")
eg <- expand.grid(prior, algorithm)
df <- list(NA)
for(i in 1:nrow(eg)){
  prior <- eg$Var1[i]
  algorithm <- eg$Var2[i]
  load(paste0("./results/fit_adapt_stan_", prior, "_", algorithm, ".RData"))
  
  # extract estimates and hpd
  summ <- summary(fit)$summary
  summ <- summ[grep("L_main_C|L_cross_C|phi_C", rownames(summ)), ] # not needed after rerun
  est <- summ[, "mean"]
  lower <- summ[, "2.5%"]
  upper <- summ[, "97.5%"]
  # TODO: check if parnms order is still correct after rerun
  parnms <- c(paste0("F1=~y", 1:18), paste0("F2=~y", 19:24), paste0("F3=~y", 25:36),
              paste0("F1=~y", 19:36), paste0("F2=~y", 1:18), paste0("F2=~y", 25:36),
              paste0("F3=~y", 1:24), "F1~~F1", "F1~~F2", "F1~~F3", "F2~~F1", "F2~~F2",
              "F2~~F3", "F3~~F1", "F3~~F2", "F3~~F3")
  df[[i]] <- cbind.data.frame("par" = parnms, "postmean" = est, "lower" = lower,
                              "upper" = upper, "prior" = prior, "algorithm" = algorithm)
}

out <- do.call(rbind, df)
save(out, file = "./results/output_stan.RData")


