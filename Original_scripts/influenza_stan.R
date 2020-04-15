library(deSolve)
library(dplyr)
library(rstan)
library(outbreaks)

rstan_options(auto_write = TRUE)

options(mc.cores = parallel::detectCores())

onset <- influenza_england_1978_school$date
cases <- influenza_england_1978_school$in_bed
N = length(onset)
pop = 763
sample_time=1:N

flu_data = list(n_obs =N,
                n_theta = 2,
                n_difeq = 3,
                n_pop = pop,
                y = cases,
                t0 = 0,
                ts = sample_time
            )

parameters = c("y_hat", "y_init", "theta", "R_0")

n_chains = 5
n_warmups=500
n_iter = 100500
n_thin = 50
set.seed(1234)

ini = function(){
  list(theta=c(runif(1,0,5), runif(1,0.2,0.4)),
       S0=runif(1, (pop-3)/pop, (pop-1)/pop))
}


nuts_fit = stan(file = "infectious.stan",
                data = flu_data,
                pars = parameters,
                init=ini,
                chains=n_chains,
                warmup = n_warmups,
                iter = n_iter,
                thin=n_thin,
                seed=13219 #period for saving samples
              )

print(nuts_fit)
nuts_fit_summary <- summary(nuts_fit, pars = c("lp_", "theta[1]", 'theta[2]', 'y_init[1]', 'R_0'))$summary
print(nuts_fit_summary, scientific=FALSE, digits=2)

#Obtain the generated samples:
posts <- rstan::extract(nuts_fit)
sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)

check_divergences(nuts_fit) 
check_treedepth(nuts_fit)

library(bayesplot)
posterior <- as.array(nuts_fit)
mcmc_trace(posterior_1, pars=c("lp_", "theta[1]", "theta[2]", "y_init[1]", "R_0"))

pairs(nuts_fit_2, pars=c("theta[1]", "theta[2]", "y_init[1]"), labels=c("beta","gamma", "s(0)"), cex.labels = 1.5, font.labels = 9, condition = "accept_stat__")

ini_vb = function(){
  list(params =c(runif(1,1.85, 1.92), runif(1, 0.47, 0.49)), S0 = runif(1, (pop-2)/pop, (pop-1)/pop))
}

fit_vb=vb(mod, data=flu_data, pars= parameters, init=ini_vb, iter = 10000, tol_rel_obj = 0.001, seed=16735679)

print(vb_fit)
vb_fit_summary <- summary(vb_fit, pars = c("theta[1]", "theta[2]", "y_init[1]", "R_0"))$summary
print(vb_fit_summary, scientific = FALSE, digits=2)

#Extract the approximate samples"
posts_vb <- rstan::extract(vb_fit)

