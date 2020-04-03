# Copyright 2020 Google LLC
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# https://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# Prior-posterior coverage test for the SEIR model.
# See http://modernstatisticalworkflow.blogspot.com/2017/04/an-easy-way-to-simulate-fake-data-from.html for
# details on the technique and motivation.

library(rstan)
library(tidyverse)

options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)

source("corona_settings.R")

# Prepare to generate fake data, using draws from the priors.
n_simulated_obs <- 25
n_forecast <- 2
population <- 1000000
sim_data <- list(n_obs=n_simulated_obs,
                  n_forecast=n_forecast,
                  n_theta=3,
                  n_difeq=4,
                  n_pop=population,
                  y=rep(1, n_simulated_obs),
                  t0=0,
                  ts=1:2,
                  tf=1:n_forecast,
                  fit_model=0)
n_simulations <- 20
n_sim_warmup <- 1000
simulated_model <- stan("sir.stan", data=sim_data, init=initial_values_function(n_simulated_obs),
                        chains=1, warmup=n_sim_warmup, iter=n_simulations + n_sim_warmup)
simulated_draws <- rstan::extract(simulated_model)

# Initialise empty storage for rank statistics and effective sample sizes.
create_empty_df <- function() {
  df <- as.data.frame(simulated_model)[0,]
  for (i in 1:n_simulations) {
    df <- df %>% rbind(as.data.frame(simulated_model)[1,] %>% mutate_all(function(c) {NA}))
  }
  return(df)
}
rank_statistics <- create_empty_df()
eff <- create_empty_df()

# For each draw from the prior, fit the model and summarise the posterior draws.
for (sim_i in 1:n_simulations) {
  fake_data <- sim_data
  fake_data$y <- simulated_draws$y_sim[sim_i,]
  fitted_model <- stan("sir.stan", data=fake_data, init=initial_values_function(n_simulated_obs),
                       chains=n_chains, iter=n_warmup_draws + n_posterior_draws_per_chain,
                       warmup=n_warmup_draws)
  
  # Find the percentile that each actual parameter value occupies in the posterior draws.
  fitted_values <- as.data.frame(fitted_model)
  true_values <- as.data.frame(simulated_model)[sim_i,]
  new_rank_statistics <- as.data.frame(simulated_model)[1,] %>% mutate_all(function(c) {NA})
  for (param_i in seq_along(fitted_values)) {
    param_name <- names(fitted_values)[param_i]
    param_draws <- fitted_values[[param_name]]
    new_rank_statistics[[param_name]] <- mean(true_values[[param_name]] > param_draws)
  }
  rank_statistics[sim_i, ] <- new_rank_statistics
  
  # Collate the effective sample sizes for each parameter.
  fitted_eff <- summary(fitted_model)$summary[,"n_eff"]
  for (param_i in seq_along(fitted_values)) {
    param_name <- names(fitted_values)[param_i]
    eff[sim_i, param_name] <- fitted_eff[param_name]
  }
}

# If the model is well identified, then the rank statistics will be uniformly distributed.
ks_stats <- data.frame(
  parameter=character(),
  ks_p_value=numeric()
)
for (param_i in seq_along(rank_statistics)) {
  ks_stats <- rbind(ks_stats,
                    data.frame(
                      parameter=names(rank_statistics)[param_i],
                      ks_p_value=ks.test(rank_statistics[, param_i], punif)$p.value
                    ))
}
ks_warning_threshold <- 0.1
expected_rejections <- ks_warning_threshold * length(names(fitted_values))
poorly_estimated_params <- ks_stats %>% arrange(ks_p_value) %>% filter(ks_p_value < ks_warning_threshold)
print(paste("There were", nrow(poorly_estimated_params), "poorly estimated parameters. (Expected: ",
            round(expected_rejections), ")"))
if (nrow(poorly_estimated_params) > expected_rejections) {
  poorly_estimated_params %>% head(10) %>% as.tbl() %>% print()
  stop("TEST FAILURE: Too many unidentified parameters")
}

# Check that all parameters are likely to be sampled efficiently.
low_effs <- data.frame(
  parameter=character(),
  low_eff=numeric()
)
for (param_i in seq_along(eff)) {
  # Replace total failures with 0.
  eff_values <- eff[, param_i]
  eff_values[is.na(eff_values)] <- 0
  low_effs <- rbind(low_effs,
                    data.frame(
                      parameter=names(eff)[param_i],
                      low_eff=quantile(eff_values, probs=(0.25))
                    ))
}
eff_warning_threshold <- 0.25 * n_posterior_draws
poorly_mixed_params <- low_effs %>% arrange(low_eff) %>% filter(low_eff < eff_warning_threshold)
print(paste("There were", nrow(poorly_mixed_params), "inefficiently mixed parameters."))
if (nrow(poorly_mixed_params) > 25) {
  poorly_mixed_params %>% head(10) %>% as.tbl() %>% print()
  stop("TEST FAILURE: parameter mixing is not good enough")
}

