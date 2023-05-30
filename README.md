# BRSsoftware
A repository for the project "Bayesian regularized SEM: Current capabilities and constraints" in which an overview of literature and software on Bayesian regularized SEM is provided, with an illustration of various software packages. The preprint is available here: https://psyarxiv.com/92vh8

All code for running the empirical example and recreating the figures in the manuscript is available in the "code" directory.

"densitplots.R" and "visualization_shrinkage.R" recreate two general figures.

"stan_conv_fun.R" and "stan_plot_fun.R" contain helper functions to extract and plot results from stanfit objects.

For the empirical example, first run "adapt_dataprep.R", then run all the analyses before combining the results and recreating the figures (using "adapt_results.R").

For the manual implementation in Stan, the models are available in the subdirectory "stanmodels". Note that these models are only valid for the factormodel considered in the manuscript and require adaptation for other models. 

Please contact Sara van Erp (sara.vanerp@gmail.com) if you have any questions about this code.
