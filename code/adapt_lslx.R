## Bayesian Regularized SEM Software: Analysis with lslx
## Author: Sara van Erp

library(lslx)
library(lavaan)

## Data: already standardized
load("./data/testDatSD_adapt.dat")

## Run analysis with lslx
# based on this tutorial: https://cran.r-project.org/web/packages/lslx/vignettes/factor-analysis.html
# default settings are used for delta and lambda grid

# specify model
modAllCL <- paste(c(paste("F1 =~ y", 1:18, " \n"), # main loadings: freely estimated
                    paste("pen() * F1 =~ y", 19:36," \n"), # cross-loadings: penalized
                    paste("F2 =~ y", 19:24," \n"),
                    paste("pen() * F2 =~ y", 1:18," \n"),
                    paste("pen() * F2 =~ y", 25:36," \n"),
                    paste("F3 =~ y", 25:36," \n"),
                    paste("pen() * F3 =~ y", 1:24," \n"),
                    "F1  ~~  1 * F1 \n", # fix latent variances for identification
                    "F2  ~~ 1 * F2 \n",
                    "F3  ~~ 1 * F3")
)

# fit the model using lasso penalty
lslx_lasso <- plsem(model = modAllCL, 
                    data = testDatSD,
                    penalty_method = "lasso")
save(lslx_lasso, file = "./results/fit_adapt_lslx_lasso.RData")

# fit the model using elastic net penalty
lslx_en <- plsem(model = modAllCL, 
                 data = testDatSD,
                 penalty_method = "elastic_net")
save(lslx_en, file = "./results/fit_adapt_lslx_en.RData")

# fit the model using mcp penalty
lslx_mcp <- plsem(model = modAllCL, 
                 data = testDatSD,
                 penalty_method = "mcp")
save(lslx_mcp, file = "./results/fit_adapt_lslx_mcp.RData")

## Extract results
res_fun <- function(lslx_fit, selector = "bic", penalty){
  est <- lslx_fit$extract_coefficient(selector = selector)
  ci <- lslx_fit$test_coefficient(selector, alpha_level = 0.05)
  lower <- ci$lower
  upper <- ci$upper
  
  # adapt names to be comparable
  spec <- lslx_fit$extract_specification()
  spec_noint <- spec[!spec$right == 1, ] # remove intercepts
  spec_load <- spec_noint[grep("<-F", spec_noint$relation), ]
  spec_load$par <- paste0(spec_load$right, "=~", spec_load$left)
  spec_var <- spec_noint[grep("<->", spec_noint$relation), ]
  spec_var$par <- paste0(spec_var$right, "~~", spec_var$left)
  spec <- rbind.data.frame(spec_load, spec_var)
  spec$oldPar <- rownames(spec)
  df <- spec[, c("oldPar", "par")]
  
  # combine
  df$postmean <- est[match(df$oldPar, names(est))] # is no posterior mean, but colnames need to be the same across methods
  df$lower <- lower[match(df$oldPar, rownames(ci))]
  df$upper <- upper[match(df$oldPar, rownames(ci))]
  df$prior <- penalty
  df$algorithm <- "lslx"
  df <- df[, -c(grep("oldPar", colnames(df)))]
  return(df)
}

# lasso
penalty <- "lasso"
load(paste0("./results/fit_adapt_lslx_", penalty, ".RData"))
df1 <- res_fun(lslx_fit = lslx_lasso, selector = "bic", penalty = penalty)

# elastic net
penalty <- "en"
load(paste0("./results/fit_adapt_lslx_", penalty, ".RData"))
df2 <- res_fun(lslx_fit = lslx_en, selector = "bic", penalty = penalty)

# mcp
penalty <- "mcp"
load(paste0("./results/fit_adapt_lslx_", penalty, ".RData"))
df3 <- res_fun(lslx_fit = lslx_mcp, selector = "bic", penalty = penalty)

# combine
out <- rbind.data.frame(df1, df2, df3)
save(out, file = "./results/output_lslx.RData")


