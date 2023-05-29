## Bayesian Regularized SEM Software: Analysis with blavaan
## Author: Sara van Erp

library(blavaan)
options(mc.cores = 3)
future::plan("multicore")
library(ggplot2)

source("./code/stan_conv_fun.R")

set.seed(09052023)

## Data: already standardized
load("./data/testDatSD_adapt.dat")

## Run model with blavaan and different prior variances
run.blavaan <- function(dat, priorvar, n.chains = 3, sample = 5000, algorithm = c("mcmc", "vb")){
  priorsd <- sqrt(priorvar)
  priorCrossL <- paste0("normal(0,", priorsd, ")")
  priorMainL <- "normal(0,10)" # default
  
  # specify model with all cross-loadings
  modAllCL <- paste(c(paste0("F1 =~ prior(\"", priorMainL, "\")*", c(paste0("y", 1:18))," \n"),
                      paste0("F1 =~ prior(\"", priorCrossL, "\")*", c(paste0("y", 19:36))," \n"),
                      paste0("F2 =~ prior(\"", priorMainL, "\")*", c(paste0("y", 19:24))," \n"),
                      paste0("F2 =~ prior(\"", priorCrossL, "\")*", c(paste0("y", 1:18))," \n"),
                      paste0("F2 =~ prior(\"", priorCrossL, "\")*", c(paste0("y", 25:36))," \n"),
                      paste0("F3 =~ prior(\"", priorMainL, "\")*", c(paste0("y", 25:36))," \n"),
                      paste0("F3 =~ prior(\"", priorCrossL, "\")*", c(paste0("y", 1:24)))
  ))
  
  fit <- bcfa(modAllCL, data = dat, n.chains = n.chains, sample = sample, target = algorithm)
  save(fit, file = paste0("./results/fit_adapt_blavaan_", priorvar, "_", algorithm, ".RData"))
}

run.blavaan(testDatSD, priorvar = 0.001, algorithm = "mcmc")
run.blavaan(testDatSD, priorvar = 0.001, algorithm = "vb")

run.blavaan(testDatSD, priorvar = 0.01, algorithm = "mcmc")
run.blavaan(testDatSD, priorvar = 0.01, algorithm = "vb")

run.blavaan(testDatSD, priorvar = 0.1, algorithm = "mcmc")
run.blavaan(testDatSD, priorvar = 0.1, algorithm = "vb")

## Check convergence 
cutoffs <- list("rhat" = 1.01, "n_eff" = 1000, "ess_bulk" = 1000, "ess_tail" = 1000, "div" = 0)

# full Bayesian models
stanls.mcmc <- list(NA)
priorvar = c(0.001, 0.01, 0.1)
for(i in 1:length(priorvar)){
  algorithm <- "mcmc"
  load(paste0("./results/fit_adapt_blavaan_", priorvar[i], "_", algorithm, ".RData"))
  stanobj <- blavInspect(fit, "mcobj")
  stanls.mcmc[[i]] <- stanobj
}

names(stanls.mcmc) <- c("fit_0.001_mcmc", "fit_0.01_mcmc", "fit_0.1_mcmc")

stan_conv_fun(stanls.mcmc, cutoffs, algorithm = "mcmc") # convergence seems fine

# Approximate Bayesian models
stanls.vb <- list(NA)
for(i in 1:length(priorvar)){
  algorithm <- "vb"
  load(paste0("./results/fit_adapt_blavaan_", priorvar[i], "_", algorithm, ".RData"))
  stanobj <- blavInspect(fit, "mcobj")
  stanls.vb[[i]] <- stanobj
}

names(stanls.vb) <- c("fit_0.001_vb", "fit_0.01_vb", "fit_0.1_vb")

stan_conv_fun(stanls.vb, cutoffs, algorithm = "vb") # convergence seems fine

## Extract results
priorvar <- c(0.001, 0.01, 0.1)
algorithm <- c("mcmc", "vb")
eg <- expand.grid(priorvar, algorithm)
df <- list(NA)
for(i in 1:nrow(eg)){
  priorvar <- eg$Var1[i]
  algorithm <- eg$Var2[i]
  load(paste0("./results/fit_adapt_blavaan_", priorvar, "_", algorithm, ".RData"))
  
  # extract estimates and hpd
  est <- blavInspect(fit, what = "postmean")
  hpd <- blavInspect(fit, what = "hpd") # 95% hpd
  df[[i]] <- cbind.data.frame("par" = rownames(hpd), "postmean" = est, hpd, "prior" = priorvar, "algorithm" = algorithm)
}

out <- do.call(rbind, df)
save(out, file = "./results/output_blavaan.RData")

## Check fit indices
load("./results/fit_adapt_blavaan_0.1_mcmc.RData")
fit1 <- fit
load("./results/fit_adapt_blavaan_0.01_mcmc.RData")
fit01 <- fit
load("./results/fit_adapt_blavaan_0.001_mcmc.RData")
fit001 <- fit

# posterior predictive model checking
ppmc1 <- ppmc(fit1)
ppmc01 <- ppmc(fit01)
ppmc001 <- ppmc(fit001)
ppmc <- list(ppmc1, ppmc01, ppmc001)
names(ppmc) <- c("var0 = 0.1", "var0 = 0.01", "var0 = 0.001")
save(ppmc, file = "./results/ppmc_blavaan.RData")

lapply(ppmc, function(x){
  summary(x, central.tendency = c("mean","median","mode"), prob = .90)
})

# specify null model for incremental fit indices
mod0 <- paste0(paste0("y", 1:36), "~~", paste0("y", 1:36), " \n")
fit0 <- bcfa(mod0, data = testDatSD, sample = 5000)

# compute fit indices
FI1 <- blavFitIndices(fit1, baseline.model = fit0)
FI01 <- blavFitIndices(fit01, baseline.model = fit0)
FI001 <- blavFitIndices(fit001, baseline.model = fit0)
fitind <- list(FI1, FI01, FI001)
names(fitind) <- c("var0 = 0.1", "var0 = 0.01", "var0 = 0.001")
save(fitind, file = "./results/fit_indices_blavaan.RData")

lapply(fitind, function(x) x)

out <- lapply(fitind, function(x){
  summary(x, central.tendency = c("mean","median","mode"), prob = .90)
})

out[[1]]$prior <- "0.1"
out[[1]]$measure <- rownames(out[[1]])
out[[2]]$prior <- "0.01"
out[[2]]$measure <- rownames(out[[2]])
out[[3]]$prior <- "0.001"
out[[3]]$measure <- rownames(out[[3]])

fitind.df <- rbind.data.frame(out[[1]], out[[2]], out[[3]])

pos <- position_dodge(width = 0.2)

# For easier comparison, group BRMSEA with BMc
# Interestingly, BMc is way off (1 = good fit), but the others are quite ok
sel <- fitind.df[(which(fitind.df$measure %in% c("BRMSEA", "BMc"))), ]
ggplot(sel, aes(x = EAP, y = measure, color = prior)) +
  geom_point(position = pos) +
  geom_linerange(aes(xmin = lower, xmax = upper), position = pos) +
  xlab("Mean estimate") + ylab("") + theme_bw(base_size = 15) +
  theme(legend.title = element_blank())

sel <- fitind.df[(-c(which(fitind.df$measure %in% c("BRMSEA", "BMc")))), ]
ggplot(sel, aes(x = EAP, y = measure, color = prior)) +
  geom_point(position = pos) +
  geom_linerange(aes(xmin = lower, xmax = upper), position = pos) +
  xlab("Mean estimate") + ylab("") + theme_bw(base_size = 15) +
  theme(legend.title = element_blank())
