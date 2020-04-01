// Copyright 2020 Google LLC
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
//     https://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// A Stan version of an SEIR model with a time-varying contact parameter.
// See https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SEIR_model -- here, the beta parameter
// is allowed to vary via an AR(1).
// The initial version was based on code generously provided by https://github.com/anastasiachtz/COMMAND_stan to
// accompany Chatzilena1a et al, "Contemporary statistical inference for infectious disease models using Stan", Epidemics (2019),
// available at https://arxiv.org/pdf/1903.00423.pdf.

functions {
  // Deterministic ODE of the SEIR model. Function signature matches rk45 solver's expectation.
  real[] SIR(
    real t,       // Time period; unused.
    real[] y,     // System state {susceptible,exposed,infected,recovered}.
    real[] theta, // Parameters.
    real[] x_r,   // Continuous-valued data; unused.
    int[] x_i     // Integer-valued data; unused.
  ) {
    real dy_dt[4];
    dy_dt[1] = -theta[1] * y[1] * y[3];
    dy_dt[2] = theta[1] * y[1] * y[3] - theta[2] * y[2];
    dy_dt[3] = theta[2] * y[2] - theta[3] * y[3];
    dy_dt[4] = theta[3] * y[3];
    return dy_dt;
  }
}

data {
  int<lower=1> n_obs;       // Number of days observed.
  int<lower=0> n_forecast;  // Number of days to project.
  int<lower=1> n_theta;     // Number of model parameters.
  int<lower=1> n_difeq;     // Number of differential equations.
  int<lower=1> n_pop;       // Population.
  int y[n_obs];             // Data: total number of infected individuals each day
  real t0;      // Initial time tick for the ODE solver. Must be provided in the data block.
  real ts[2];   // Time ticks for the ODE solver. Must be provided in the data block.
  real tf[n_forecast];  // Time ticks for the forward projections.
}
  
transformed data {
  // Covariates for the ODE solver. Not used here, but must be provided.
  real x_r[0];
  // Integer covariates for the ODE solver. Not used here, but must be provided.
  int x_i[0];
}
  
parameters {
  real<lower=0> theta[n_theta];   // ODE model parameters.
  real<lower=0> I0;               // Initial number of infected individuals
  real<lower=0> E0;               // Initial number of exposed individuals.
  real<lower=0, upper=1> phi;     // Contact-rate AR(1) parameter.
  real<lower=0> sigma;            // Std dev of innovations in the contact rate.
  real beta_t[n_obs];             // Time-varying contact rate.
}
  
transformed parameters{
  real y_init[n_difeq];       // Initial conditions for both fractions of S and I
  real y_hat[n_obs, n_difeq];
  real dy_dt[n_obs, n_difeq];
  
  y_init[1] = 1.0 - I0/n_pop - E0/n_pop;
  y_init[2] = E0/n_pop;
  y_init[3] = I0/n_pop;
  y_init[4] = 0;
  
  // Use a deterministic ODE solver to go, in continuous time, from time `t` to `t+1`,
  // conditional on the contact-rate parameter beta_t.
  dy_dt[1,] = SIR(0, y_init, theta, x_r, x_i);
  for (i in 1:n_difeq) {
    y_hat[1,i] = y_init[i] + dy_dt[1,i];
  }
  for (t in 2:n_obs) {
    dy_dt[t,] = SIR(0, y_hat[t-1,], {beta_t[t], theta[2], theta[3]}, x_r, x_i);
    for (i in 1:n_difeq) {
      y_hat[t,i] = y_hat[t-1,i] + dy_dt[t,i];
    }
  }
}
  
model {
  // Poisson rate parameter at each time point.
  real lambda[n_obs];
  
  // Prior on steady-state contact rate: given 7-day infectious period, R0 ~ 2.
  theta[1] ~ normal(2.0/7.0, 0.05);
  // Prior on mean incubation period: 5 days.
  theta[2] ~ normal(1.0/5.0, 1.0);
  // Prior on infectious period: 7 days.
  theta[3] ~ normal(1.0/7.0, 0.01);
  // Prior on t0 incubation number: 300 ± 200. Note that data starts with 100 measured infections.
  E0 ~ normal(300, 100);
  // Prior on initial infections number: 100 ± 5.
  I0 ~ normal(100, 2.5);
  // Prior on changes to the contact rate: up to 5% of the steady-state rate.
  sigma ~ normal(0, 0.015);
  // Prior on regression to the mean contact rate.
  phi ~ beta(0.5, 0.5);
  
  // Define the AR(1) process for the contact-rate parameter.
  beta_t[1] ~ normal(theta[1], sigma);
  for (i in 2:n_obs) {
    beta_t[i] ~ normal(phi*beta_t[i-1] + (1-phi)*theta[1], sigma);
  }
  
  // Define the observation process.
  for (i in 1:n_obs){
    // Public datasets report cumulative confirmed cases, which is I+R in this model.
    lambda[i] = (y_hat[i,3] + y_hat[i,4]) * n_pop;
  }
  y ~ poisson(lambda);
}
  
generated quantities {
  // The solution from the deterministic ODE solver.
  real y_forc[n_forecast, n_difeq];
  // The initial condition for the forward projection = the end state of the data.
  real y_init_forc[n_difeq];
  y_init_forc[1] = y_hat[n_obs, 1];
  y_init_forc[2] = y_hat[n_obs, 2];
  y_init_forc[3] = y_hat[n_obs, 3];
  y_init_forc[4] = y_hat[n_obs, 4];
  
  // Generate forward projections deterministically, conditional on the final value of the contact-rate parameter.
  y_forc = integrate_ode_rk45(SIR, y_init_forc, t0, tf, {beta_t[n_obs], theta[2], theta[3]}, x_r, x_i);
}
