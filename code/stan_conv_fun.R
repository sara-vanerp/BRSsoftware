# See also this post: https://mc-stan.org/misc/warnings.html

# Function to automatically check convergence for one or more stanfit objects
# Returns a list with summarized convergence information per stanfit object
# Note: divergent transitions and n_eff are not relevant when algorithm = vb

#stanls <- list("fit1"= stanfit, "fit2" = stanvi)
#cutoffs <- list("rhat" = 1.01, "n_eff" = 1000, "ess_bulk" = 1000, "ess_tail" = 1000, "div" = 0)
#stan_conv_fun(stanls, cutoffs)

stan_conv_fun <- function(stan_fits, cutoffs, algorithm){
  out = lapply(stan_fits, function(x){

    summ_draws = posterior::summarise_draws(x)
    
    rhat = summ_draws$rhat
    names(rhat) = summ_draws$variable
    rhat_out = rhat[which(rhat > cutoffs$rhat)]
    
    n_eff = summary(x)$summary[, "n_eff"]
    n_eff_out = n_eff[which(n_eff < cutoffs$n_eff)]
    
    ess_bulk = summ_draws$ess_bulk
    names(ess_bulk) = summ_draws$variable
    ess_bulk_out = ess_bulk[which(ess_bulk < cutoffs$ess_bulk)]
    
    ess_tail = summ_draws$ess_tail
    names(ess_tail) = summ_draws$variable
    ess_tail_out = ess_tail[which(ess_tail < cutoffs$ess_tail)]
    
    if(algorithm == "mcmc"){
      sp = rstan::get_sampler_params(x, inc_warmup=FALSE)
      div = sum(sapply(sp, function(x) sum(x[, "divergent__"])))
      div_out = ifelse(div <= cutoffs$div, FALSE, TRUE)
    } else if(algorithm == "vb"){
      div = NA
      div_out = NA
    }
    
    df_res = cbind.data.frame("max_rhat" = max(rhat, na.rm = TRUE),
                              "min_n_eff" = min(n_eff, na.rm = TRUE),
                              "min_ess_bulk" = min(ess_bulk, na.rm = TRUE),
                              "min_ess_tail" = min(ess_tail, na.rm = TRUE),
                              "max_div" = max(div, na.rm = TRUE))
    
    ls_cutoff = list("rhat_exceed_cutoff" = rhat_out,
                     "n_eff_exceed_cutoff" = n_eff_out,
                     "ess_bulk_exceed_cutoff" = ess_bulk_out,
                     "ess_tail_exceed_cutoff" = ess_tail_out,
                     "div_exceed_cutoff" = div_out)
    
    res = list("conv_crit" = df_res, ls_cutoff)

    return(res)
    
  })
  
  names(out) = names(stan_fits)
  return(out)
}






