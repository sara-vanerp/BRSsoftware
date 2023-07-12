## Prior-likelihood-posterior plot to visualize behavior shrinkage priors
library(ggplot2)

# data
n <- 1
md <- 0.15
sdd <- 0.1
x <- rnorm(n, mean = md, sd = sdd)

# prior
m0 <- 0
sd01 <- 0.01
sd02 <- 0.05
sd03 <- 0.15

# posterior
sd0 <- sd03 # change per prior
post.mean <- ((n/(sdd^2) + 1/(sd0^2))^-1)*(n/(sdd^2) * mean(x) + 1/(sd0^2)*m0)
post.sd <- sqrt((n/(sdd^2) + 1/(sd0^2))^-1)

# plot
g <- ggplot(data = data.frame(x = c(-0.5, 0.5)), aes(x)) + 
  stat_function(fun = dnorm, n = 200, args = list(mean = m0, sd = sd0), aes(linetype = "Prior"), linewidth = 1) +
  stat_function(fun = dnorm, n = 200, args = list(mean = md, sd = sdd), aes(linetype = "Likelihood"), linewidth = 1) +
  stat_function(fun = dnorm, n = 200, args = list(mean = post.mean, sd = post.sd), aes(linetype = "Posterior"), linewidth = 1) +
  scale_linetype_manual("", values = c(3, 1, 2)) +
  labs(x = "", y = "") +
  theme_bw(base_size = 25) + 
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        legend.position = "none",
        legend.key.size = unit(1.5, 'cm'))

# save
png("./results/shrinkage_prior3.png", width = 800, height = 800)
g + ggtitle(paste0("Prior standard deviation ", sd0))
dev.off()

