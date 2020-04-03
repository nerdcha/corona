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


# Initial values and sampler settings for Stan, extracted here so they
# can be shared between production and testing.

n_chains <- 6
n_posterior_draws_per_chain <- 1000
n_posterior_draws <- n_chains * n_posterior_draws_per_chain
n_warmup_draws <- 1000

# A function-factory that returns an initial value function.
# (Some initial values depend on the number of obs., which varies between
# production and tests.)
initial_values_function <- function(n_obs) {
  function() {
    list(theta=c(runif(1, 0.2, 0.35), runif(1, 0.2, 0.4), runif(1, 0.2, 0.4)), 
         I0=runif(1, 90, 150),
         E0=runif(1, 50, 1000),
         R0=runif(1, 0, 20),
         phi=runif(1, 0.3, 0.9),
         sigma=runif(1, 0.002, 0.02),
         log_beta_t_deviation=runif(n_obs, -0.05, 0.05)
    )
  }
}
