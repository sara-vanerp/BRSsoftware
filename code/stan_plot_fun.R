# Function to combine results from one or more stanfit objects
# Provide a list of stan fitobjects and one or more parameters 
# Note: a selection of parameters to plot must be provided, to avoid models having different
# numbers of parameters, which messes up the naming of the models

#stanls <- list("fit1"= stanfit, "fit2" = stanvi)

#p <- stan_plot_fun(stanls, pars = c("L_cross_C|L_main_C"))
#plot(p)

stan_plot_fun <- function(stan_fits, pars, dodge = 0.2){
  out = lapply(stan_fits, function(x){
  
    post_draws = as.matrix(x)
    # TODO: use summary within rstan instead of bayesplot to avoid unnecessary dependency
    # add option to specify quantiles and type of point estimate (+ option to plot both)
    mcmc_int = bayesplot::mcmc_intervals_data(post_draws)
  })

  combined = do.call(rbind, out)
  combined_sel = combined[grep(pars, combined$parameter), ]
  combined_sel$model <- rep(names(stan_fits), each = nrow(combined_sel)/2)
  
  pos = ggplot2::position_dodge(width = dodge)
  
  plot_ls = ggplot2::ggplot(combined_sel, aes(x = m, y = parameter, color = model)) +
    geom_point(position = pos) +
    geom_linerange(aes(xmin = ll, xmax = hh), position = pos) +
    xlab("Mean estimate")
  
  return(plot_ls)

}



