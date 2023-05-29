## Bayesian Regularized SEM Software: Combine all results
## Author: Sara van Erp

library(psych)
library(ggplot2)
library(bayesplot)

load("./results/output_blavaan.RData")
res1 <- out
levels(res1$algorithm) <- list(mcmc_blav = "mcmc", vb_blav = "vb")
load("./results/output_lawlb.RData")
res2 <- out
load("./results/output_stan.RData")
res3 <- out
load("./results/output_lslx.RData")
res4 <- out

# add loadings EFA for comparison
# original EFA was with 8 factors on the training data, so need to rerun but this leads to some loadings belonging to different factors
load("./data/testDatSD_adapt.dat")
fitEFA <- fa(r = testDatSD,
             nfactors = 3,
             rotate = "oblimin",
             fm = "pa")
print(fitEFA, sort = TRUE)
L <- as.data.frame(unclass(fitEFA$loadings))
L1 <- data.frame(par = paste0("F1=~", rownames(L)), postmean = L$PA1,
                 lower = NA, upper = NA, prior = "none", algorithm = "EFA")
L2 <- data.frame(par = paste0("F2=~", rownames(L)), postmean = L$PA2,
                 lower = NA, upper = NA, prior = "none", algorithm = "EFA")
L3 <- data.frame(par = paste0("F3=~", rownames(L)), postmean = L$PA3,
                 lower = NA, upper = NA, prior = "none", algorithm = "EFA")

# combine
res <- rbind.data.frame(res1, res2, res3, res4, L1, L2, L3) 
out <- res

# compare blavaan algorithms
sel1 <- which(res$algorithm == "mcmc_blav")
sel2 <- which(res$algorithm == "vb_blav")
out <- res[c(sel1, sel2), ]

# compare Stan algorithms
sel3 <- which(res$algorithm == "mcmc" & res$prior == "ridge")
sel4 <- which(res$algorithm == "vb" & res$prior == "ridge")
out <- res[c(sel3, sel4), ]

sel5 <- which(res$algorithm == "mcmc" & res$prior != "ridge")
sel6 <- which(res$algorithm == "vb" & res$prior != "ridge")
out <- res[c(sel5, sel6), ]

# compare classical approaches
sel7 <- which(res$algorithm == "lslx")
out <- res[c(sel7), ]

# compare priors in blavaan
sel8 <- which(res$algorithm == "mcmc_blav")
out <- res[c(sel8, sel3), ]

# compare regularized horseshoe priors
sel9 <- which(res$prior == "reghs_def" & res$algorithm == "mcmc")
sel10 <- which(res$prior == "reghs_p0" & res$algorithm == "mcmc")
out <- res[c(sel9, sel10), ]
  
# compare lasso specifications
sel11 <- which(res$algorithm == "LAWLB")
sel12 <- which(res$algorithm == "lslx" & res$prior == "lasso")
out <- res[c(sel11, sel12), ]

# compare multiple priors and classical methods
sel13 <- which(res$algorithm == "LAWLB" & res$prior == "lasso")
sel14 <- which(res$prior == "reghs_def" & res$algorithm == "mcmc")
sel15 <- which(res$prior == "ridge" & res$algorithm == "mcmc")
sel16 <- which(res$algorithm == "lslx" & res$prior == "mcp")
sel17 <- which(res$prior == "none" & res$algorithm == "EFA")

out <- res[c(sel13, sel14, sel15, sel16, sel17), ]

out$analysis <- apply(out, 1, function(x) paste(x["prior"], x["algorithm"], sep = "_"))

pos <- position_dodge(width = 0.2)

# cross-loadings
cross_F1 <- paste0("F1=~y", 19:23, "$")
cross_F1 <- paste0("F1=~y", 23:27, "$")
cross_F1 <- paste0("F1=~y", 27:31, "$")
cross_F1 <- paste0("F1=~y", 31:36, "$")
outsel <- out[grep(paste(cross_F1, collapse = "|"), out$par), ]
ggplot(outsel, aes(x = postmean, y = par, color = analysis)) +
  geom_point(position = pos) +
  geom_linerange(aes(xmin = lower, xmax = upper), position = pos) +
  xlab("Mean estimate") + ylab("") + theme_bw(base_size = 15) +
  theme(legend.title = element_blank())

cross_F2 <- paste0("F2=~y", c(1:5), "$")
cross_F2 <- paste0("F2=~y", c(6:10), "$")
cross_F2 <- paste0("F2=~y", c(11:15), "$")
cross_F2 <- paste0("F2=~y", c(16:18,25:26), "$")
cross_F2 <- paste0("F2=~y", c(27:31), "$")
cross_F2 <- paste0("F2=~y", c(31:36), "$")
outsel <- out[grep(paste(cross_F2, collapse = "|"), out$par), ]
ggplot(outsel, aes(x = postmean, y = par, color = analysis)) +
  geom_point(position = pos) +
  geom_linerange(aes(xmin = lower, xmax = upper), position = pos) +
  xlab("Mean estimate") + ylab("") + theme_bw(base_size = 15) +
  theme(legend.title = element_blank())

cross_F3 <- paste0("F3=~y", 1:5, "$")
cross_F3 <- paste0("F3=~y", 6:10, "$")
cross_F3 <- paste0("F3=~y", 11:15, "$")
cross_F3 <- paste0("F3=~y", 16:20, "$")
cross_F3 <- paste0("F3=~y", 21:24, "$")
outsel <- out[grep(paste(cross_F3, collapse = "|"), out$par), ]
ggplot(outsel, aes(x = postmean, y = par, color = analysis)) +
  geom_point(position = pos) +
  geom_linerange(aes(xmin = lower, xmax = upper), position = pos) +
  xlab("Mean estimate") + ylab("") + theme_bw(base_size = 15) +
  theme(legend.title = element_blank())

# nuisance parameters
main_F1 <- paste0("F1=~y", 1:9, "$")
outsel <- out[grep(paste(main_F1, collapse = "|"), out$par), ]
ggplot(outsel, aes(x = postmean, y = par, color = analysis)) +
  geom_point(position = pos) +
  geom_linerange(aes(xmin = lower, xmax = upper), position = pos) +
  xlab("Mean estimate") + ylab("") + theme_bw(base_size = 15) +
  theme(legend.title = element_blank())

main_F1 <- paste0("F1=~y", 10:18, "$")
outsel <- out[grep(paste(main_F1, collapse = "|"), out$par), ]
ggplot(outsel, aes(x = postmean, y = par, color = analysis)) +
  geom_point(position = pos) +
  geom_linerange(aes(xmin = lower, xmax = upper), position = pos) +
  xlab("Mean estimate") + ylab("") + theme_bw(base_size = 15) +
  theme(legend.title = element_blank())

main_F2 <- paste0("F2=~y", 19:24, "$")
outsel <- out[grep(paste(main_F2, collapse = "|"), out$par), ]
ggplot(outsel, aes(x = postmean, y = par, color = analysis)) +
  geom_point(position = pos) +
  geom_linerange(aes(xmin = lower, xmax = upper), position = pos) +
  xlab("Mean estimate") + ylab("") + theme_bw(base_size = 15) +
  theme(legend.title = element_blank())


main_F3 <- paste0("F3=~y", 25:36, "$")
outsel <- out[grep(paste(main_F3, collapse = "|"), out$par), ]
ggplot(outsel, aes(x = postmean, y = par, color = analysis)) +
  geom_point(position = pos) +
  geom_linerange(aes(xmin = lower, xmax = upper), position = pos) +
  xlab("Mean estimate") + ylab("") + theme_bw(base_size = 15) +
  theme(legend.title = element_blank())

outsel <- out[grep("F1~~F2|F1~~F3|F2~~F3", out$par), ]
ggplot(outsel, aes(x = postmean, y = par, color = analysis)) +
  geom_point(position = pos) +
  geom_linerange(aes(xmin = lower, xmax = upper), position = pos) +
  xlab("Mean estimate") + ylab("") + theme_bw(base_size = 15) +
  theme(legend.title = element_blank())

# Check AVBB21 (y36) which has a relatively high loading on both F1 and F3 according to the initial EFA.
outsel <- out[grep("=~y36", out$par), ]
ggplot(outsel, aes(x = postmean, y = par, color = analysis)) +
  geom_point(position = pos) +
  geom_linerange(aes(xmin = lower, xmax = upper), position = pos) +
  xlab("Mean estimate") + ylab("") + theme_bw(base_size = 15) +
  theme(legend.title = element_blank())

# Plots for paper

# Show per prior one parameter that is illustrative for the general pattern
# comparison reghs algorithms
out <- res[c(which(res$prior == "reghs_def")), ]
out$analysis <- factor(apply(out, 1, function(x) paste(x["prior"], x["algorithm"], sep = "_")))
levels(out$analysis) <- list("Reg. hs with VB" = "reghs_def_vb",
                              "Reg. hs with MCMC" = "reghs_def_mcmc")
outsel1 <- out[grep("F2=~y16", out$par), ]

# comparison ridge specifications
out <- res[c(which(res$algorithm %in% c("mcmc_blav", "mcmc"))), ]
out <- out[-c(which(out$prior %in% c("reghs_p0", "reghs_def"))), ]
out$analysis <- factor(apply(out, 1, function(x) paste(x["prior"], x["algorithm"], sep = "_")))
levels(out$analysis) <- list("Ridge var = 0.001" = "0.001_mcmc_blav",
                             "Ridge var = 0.01" = "0.01_mcmc_blav",
                             "Ridge var = 0.1" = "0.1_mcmc_blav",
                             "Ridge var = est." = "ridge_mcmc")
outsel2 <- out[grep("F3=~y22", out$par), ]

# comparison lasso specifications
out <- res[c(which(res$algorithm %in% c("LAWLB", "lslx"))), ]
out <- out[-c(which(out$prior %in% c("mcp", "en"))), ]
out$analysis <- factor(apply(out, 1, function(x) paste(x["prior"], x["algorithm"], sep = "_")))
levels(out$analysis) <- list("Adaptive lasso" = "alas_LAWLB",
                             "Bayesian lasso" = "lasso_LAWLB",
                             "Classical lasso" = "lasso_lslx")
outsel3 <- out[grep("F1=~y30", out$par), ]

# combine and plot
outsel <- rbind.data.frame(outsel1, outsel2, outsel3)

png("./results/selected_crossloadings.png", width = 1500, height = 1000)
ggplot(outsel, aes(y = postmean, x = analysis)) +
  geom_point(position = pos, size = 5) +
  geom_linerange(aes(ymin = lower, ymax = upper), position = pos, size = 2) +
  facet_grid(~par, scales = "free_x") +
  ylab("Mean estimate and 95% CI") + xlab("") + theme_bw(base_size = 40) +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# Show interval estimates for a selection of priors and classical approaches for multiple cross-loadings
cross_F1 <- paste0("F1=~y", 27:33, "$")
out <- res[c(sel13, sel14, sel15, sel16, sel17), ]
out$analysis <- factor(apply(out, 1, function(x) paste(x["prior"], x["algorithm"], sep = "_")))
levels(out$analysis) <- list("Minimax concave penalty" = "mcp_lslx",
                             "Bayesian lasso" = "lasso_LAWLB",
                             "EFA estimate" = "none_EFA",
                             "Reg. hs with MCMC" = "reghs_def_mcmc",
                             "Ridge var = est." = "ridge_mcmc")
outsel <- out[grep(paste(cross_F1, collapse = "|"), out$par), ]
png("./results/crossloadings.png", width = 1700, height = 1000)
ggplot(outsel, aes(y = postmean, x = par, color = analysis)) +
  geom_point(position = pos, size = 5) +
  geom_hline(aes(yintercept = 0), colour = "darkgrey") +
  geom_linerange(aes(ymin = lower, ymax = upper), position = pos, size = 2) +
  ylab("Mean estimate and 95% CI") + xlab("") + theme_bw(base_size = 40) +
  theme(legend.title = element_blank(), legend.position = "bottom")
dev.off()

# Influence on nuisance parameters
# factor correlations with same selection of priors as above + mcmc_blav options
out <- res[c(sel13, sel14, sel15, sel16, sel17, sel1), ]
out$analysis <- factor(apply(out, 1, function(x) paste(x["prior"], x["algorithm"], sep = "_")))
out <- out[-c(which(out$analysis == "0.01_mcmc_blav")), ]
levels(out$analysis) <- list("Minimax concave penalty" = "mcp_lslx",
                             "Bayesian lasso" = "lasso_LAWLB",
                             "EFA estimate" = "none_EFA",
                             "Reg. hs with MCMC" = "reghs_def_mcmc",
                             "Ridge var = est." = "ridge_mcmc",
                             "Ridge var = 0.001" = "0.001_mcmc_blav",
                             "Ridge var = 0.1" = "0.1_mcmc_blav")
outsel <- out[grep("F1~~F2|F1~~F3|F2~~F3", out$par), ]
png("./results/factorcorrelations.png", width = 1700, height = 1000)
ggplot(outsel, aes(y = postmean, x = par, color = analysis)) +
  geom_point(position = pos, size = 5) +
  geom_hline(aes(yintercept = 0), colour = "darkgrey") +
  geom_linerange(aes(ymin = lower, ymax = upper), position = pos, size = 2) +
  ylab("Mean estimate and 95% CI") + xlab("") + theme_bw(base_size = 40) +
  theme(legend.title = element_blank(), legend.position = "bottom")
dev.off()
