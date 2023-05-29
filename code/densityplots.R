## Density plots shrinkage priors
## Author: Sara van Erp

library(ggplot2)
library(gridExtra)
library(tidyr)
library(rmutil) # for double-exponential 
library(LaplacesDemon) # for horseshoe density 

set.seed(12052023)

##### Sample from prior densities -----
ndraws <- 1e+05

## ridge 
ridge <- rnorm(ndraws, mean=0, sd=0.5)
ridge2 <- rep(NA, ndraws)
for(i in 1:ndraws){
  sd0 <- rhalfcauchy(1, scale=1)
  ridge2[i] <- rnorm(1, 0, sd0)
}

## lasso 
lasso <- rmutil::rlaplace(ndraws, m=0, s=0.5)

## horseshoe
hs <- rep(NA, ndraws)
for(i in 1:ndraws){
  lambda <- rhalfcauchy(1, scale=1)
  tau <- rhalfcauchy(1, scale=lambda)
  hs[i] <- rnorm(1, 0, tau)
}

##### Plot -----
## Create plot data
df.comb <- data.frame(ridge, ridge2, lasso, hs)
df.long <- gather(df.comb, Prior, value) # long format
df.long$Prior <- factor(df.long$Prior)
levels(df.long$Prior) <- list("Ridge fix SD"="ridge", "Ridge est SD"="ridge2", "Lasso"="lasso", "Horseshoe"="hs")
df.long$asymp <- rep(NA, nrow(df.long)) 
df.long[which(df.long$Prior=="Horseshoe"), "asymp"] <- 0

# xlim removes values outside the range, coord_cartesian does not, but behaves strange sometimes
# solution: restrict range with xlim first & then use coord-cartesian  
sub1 <- subset(df.long, Prior %in% c("Ridge fix SD", "Lasso", "Horseshoe"))
p1 <- ggplot(sub1, aes(x = value, linetype = Prior)) +
  stat_density(geom = "line", position = "identity") +
  geom_vline(aes(xintercept = asymp), colour = "grey") +
  scale_linetype_manual(values = c(1, 2, 3)) +
  xlim(-10,10) + coord_cartesian(xlim=c(-5,5)) + ylim(0, 1.5) +
  theme_bw(base_size = 15) + xlab("") + ylab("") +
  theme(axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank()) 

png("./results/densities_sel.png", width = 400, height = 400)
p1
dev.off()

# Compare ridge with estimated SD to horseshoe
sub2 <- subset(df.long, Prior %in% c("Ridge fix SD", "Ridge est SD"))
p2 <- ggplot(sub2, aes(x = value, linetype = Prior)) +
  stat_density(geom = "line", position = "identity") +
  geom_vline(aes(xintercept = asymp), colour = "grey") +
  scale_linetype_manual(values = c(1, 2, 3)) +
  xlim(-10,10) + coord_cartesian(xlim=c(-5,5)) + ylim(0, 1.5) +
  theme_bw(base_size = 15) + xlab("") + ylab("") +
  theme(axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank()) 

png("./results/comparison_ridge.png", width = 400, height = 400)
p2
dev.off()

# Combine all
p3 <- ggplot(df.long, aes(x = value, linetype = Prior, color = Prior)) +
  stat_density(geom = "line", position = "identity") +
  geom_vline(aes(xintercept = asymp), colour = "lightgrey") +
  scale_linetype_manual(values = c(1, 1, 2, 3)) +
  scale_color_manual(values = c("black", "grey", "black", "black")) +
  xlim(-10,10) + coord_cartesian(xlim=c(-5,5)) + ylim(0, 1.5) +
  theme_bw(base_size = 15) + xlab("") + ylab("") +
  theme(axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank()) 

png("./results/densities.png", width = 800, height = 800)
p3
dev.off()
