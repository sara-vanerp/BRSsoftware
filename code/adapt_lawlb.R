## Bayesian Regularized SEM Software: Analysis with LAWLB
## Author: Sara van Erp

library(LAWBL)
library(coda)
library(bayestestR)

set.seed(09052023)

# See also the quick start guide for the package here: https://cran.r-project.org/web/packages/LAWBL/vignettes/LAWBL.html

## Data: already standardized
load("./data/testDatSD_adapt.dat")

## Create design matrix
Q <- matrix(-1, nrow = ncol(testDatSD), ncol = 3) # -1 for unspecified (regularized) and 1 for specified (main loadings)
f1 <- which(colnames(testDatSD) %in% paste0("y", 1:18))
f2 <- which(colnames(testDatSD) %in% paste0("y", 19:24))
f3 <- which(colnames(testDatSD) %in% paste0("y", 25:36))
Q[f1, 1] <- Q[f2, 2] <- Q[f3, 3] <- 1

## Run the model with the lasso
fit <- pcfa(dat = testDatSD,
            Q = Q,
            LD = FALSE,
            cati = NULL,
            PPMC = FALSE,
            burn = 2000,
            iter = 5000,
            update = 1000,
            verbose = TRUE,
            rseed = 09052023)
save(fit, file = "./results/fit_adapt_lawlb_lasso.RData")

## Ideally, one would run multiple chains using the random seed and combine draws

## Run the model with the adaptive lasso
fit <- pcfa(dat = testDatSD,
            Q = Q,
            LD = FALSE,
            cati = NULL,
            PPMC = FALSE,
            alas = TRUE,
            burn = 2000,
            iter = 5000,
            update = 1000,
            verbose = TRUE,
            rseed = 09052023)
save(fit, file = "./results/fit_adapt_lawlb_alas.RData")

## Check convergence
# LA = posterior draws loadings (main and cross)
# Omega = Factor scores per person on each factor
# PSX = posterior draws item variances
# PHI = posterior draws factor covariance matrix
# gammal = ?
# gammas = ?
# Eigen = posterior draws eigenvalues
load("./results/fit_adapt_lawlb_lasso.RData")
plot_lawbl(fit)
plot_lawbl(fit, what = "density")
plot_lawbl(fit, what = "EPSR")

par(mfrow = c(2, 2))
traceplot(as.mcmc(fit$LA[, 1:108]))
densplot(as.mcmc(fit$PHI))

load("./results/fit_adapt_lawlb_alas.RData")
plot_lawbl(fit)
plot_lawbl(fit, what = "density")
plot_lawbl(fit, what = "EPSR")

par(mfrow = c(2, 2))
traceplot(as.mcmc(fit$LA[, 1:108]))
densplot(as.mcmc(fit$PHI))

## Extract results
res.fun <- function(prior){
  load(paste0("./results/fit_adapt_lawlb_", prior, ".RData"))
  summary(fit, what = "phi") # only summarizes significant loadings
  
  # loadings
  parL <- c(paste0("F1=~y", 1:36),
           paste0("F2=~y", 1:36),
           paste0("F3=~y", 1:36))
  postmeanL <- apply(fit$LA, 2, mean)
  ci_hdiL <- apply(fit$LA, 2, function(x) bayestestR::ci(x, ci = 0.95, method = "HDI"))
  ciL <- do.call(rbind.data.frame, ci_hdiL)
  lowerL <- ciL$CI_low
  upperL <- ciL$CI_high
  
  outL <- cbind.data.frame("par" = parL, 
                           "postmean" = postmeanL, 
                           "lower" = lowerL,
                           "upper" = upperL, 
                           "prior" = prior, 
                           "algorithm" = "LAWLB")
  
  # factor correlations
  parF <- c("F1~~F2", "F1~~F3", "F2~~F3")
  postmeanF <- apply(fit$PHI, 2, mean) 
  ci_hdiF <- apply(fit$PHI, 2, function(x) bayestestR::ci(x, ci = 0.95, method = "HDI"))
  ciF <- do.call(rbind.data.frame, ci_hdiF)
  lowerF <- ciF$CI_low
  upperF <- ciF$CI_high
  
  outF <- cbind.data.frame("par" = parF, 
                           "postmean" = postmeanF, 
                           "lower" = lowerF,
                           "upper" = upperF, 
                           "prior" = prior, 
                           "algorithm" = "LAWLB")
  
  out <- rbind.data.frame(outL, outF)
  
  return(out)
}

out.lasso <- res.fun("lasso")
out.alas <- res.fun("alas")
out <- rbind.data.frame(out.lasso, out.alas)

save(out, file = "./results/output_lawlb.RData")

