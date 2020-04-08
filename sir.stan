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
    real[] y,     // System state {susceptible,exposed,infected,recovered,untested-infectious,untested-recovered}.
    real[] theta, // Parameters.
    real[] x_r,   // Continuous-valued data; unused.
    int[] x_i     // Integer-valued data; unused.
    ) {
      real dy_dt[6];
      dy_dt[1] = -theta[1] * y[1] * (y[3] + y[5]);
      dy_dt[2] = theta[1] * y[1] * (y[3] + y[5]) - theta[2] * y[2] - theta[4] * y[2];
      dy_dt[3] = theta[2] * y[2] - theta[3] * y[3];
      dy_dt[4] = theta[3] * y[3];
      dy_dt[5] = theta[4] * y[2] - theta[3] * y[5];
      dy_dt[6] = theta[3] * y[5];
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
  int<lower=0, upper=1> fit_model;  // If 0, only draws from the priors will be used.
}

transformed data {
  // Covariates for the ODE solver. Not used here, but must be provided.
  real x_r[0];
  // Integer covariates for the ODE solver. Not used here, but must be provided.
  int x_i[0];
}

parameters {
  real<lower=0> R0;               // The unconditional basic reproduction number.
  real<lower=0> average_incubation_days;  // An ODE parameter: {a^-1} in the Wikipedia notation.
  real<lower=0> average_infectious_days;  // An ODE parameter: {gamma^-1} in the Wikipedia notation.
  real<lower=0, upper=1> untested_share; // Share of infectious individuals who don't appear in data.
  real<lower=0> I0;               // Initial number of infected individuals
  real<lower=0> E0;               // Initial number of exposed individuals.
  real<lower=0> Rec0;               // Initial number of recovered individuals.
  real<lower=0, upper=1> phi;     // Contact-rate AR(1) parameter.
  real sigma;                     // Std dev of innovations in the contact rate.
  real<lower=0> log_beta_t_deviation[n_obs];             // Time-varying contact rate shock.
}

transformed parameters{
  real y_init[n_difeq];       // Initial conditions for both fractions of S and I
  real y_hat[n_obs, n_difeq];
  real dy_dt[n_obs, n_difeq];
  real gamma;
  real a;
  real beta;
  
  a = 1.0/average_incubation_days;
  gamma = 1.0/average_infectious_days;
  beta = R0 * gamma;
  
  y_init[1] = 1.0 - I0/n_pop - E0/n_pop;
  y_init[2] = E0/n_pop;
  y_init[3] = (1.0 - untested_share) * I0/n_pop;
  y_init[4] = (1.0 - untested_share) * Rec0/n_pop;
  y_init[5] = untested_share * I0/n_pop;
  y_init[6] = untested_share * Rec0/n_pop;
  
  // Use a deterministic ODE solver to go, in continuous time, from time `t` to `t+1`,
  // conditional on the contact-rate parameter = exp(log_beta_t_deviation)*theta[1].
  dy_dt[1,] = SIR(0, y_init, {beta, a * (1.0 - untested_share), gamma, a * untested_share}, x_r, x_i);
  for (i in 1:n_difeq) {
    y_hat[1,i] = y_init[i] + dy_dt[1,i];
  }
  for (t in 2:n_obs) {
    dy_dt[t,] = SIR(0, y_hat[t-1,],
                    {exp(-log_beta_t_deviation[t])*beta, a  * (1.0 - untested_share), gamma, a * untested_share},
                    x_r, x_i);
    for (i in 1:n_difeq) {
      y_hat[t,i] = y_hat[t-1,i] + dy_dt[t,i];
    }
  }
}

model {
  // Poisson rate parameter at each time point.
  real lambda[n_obs];
  
  // Prior weights on ODE parameters are informed by literature on the novel coronavirus outbreak
  // in China.
  R0 ~ normal(3.5, 0.5);
  average_incubation_days ~ normal(5.5, 1.5);
  // Recovery reports from VIC and WA lag case reports by ~7 days and ~10 days respectively.
  average_infectious_days ~ normal(7.0, 0.75);
  
  // Heavy-tailed priors on initial conditions.
  E0 ~ student_t(3, 0, 100);
  I0 ~ student_t(3, 0, 100);
  Rec0 ~ student_t(3, 0, 100);
  
  // Prior on log changes to the contact rate.
  sigma ~ student_t(3, 0, 0.02);
  // Prior on regression to the mean contact rate.
  phi ~ beta(20.0, 2.0);
  
  // Prior on time-invariant testing policy/prevalence.
  untested_share ~ beta(2.0, 2.0);
  
  // Define the AR(1) process for the percent deviation in the contact-rate parameter.
  log_beta_t_deviation[1] ~ normal(0, sigma);
  for (i in 2:n_obs) {
    log_beta_t_deviation[i] ~ normal(phi*log_beta_t_deviation[i-1], sigma);
  }
  
  // Define the observation process.
  if (fit_model == 1) {
    for (i in 1:n_obs){
      // Public datasets report cumulative confirmed cases, which is I+R in this model.
      lambda[i] = (y_hat[i,3] + y_hat[i,4]) * n_pop;
    }
    y ~ poisson(lambda);
  }
}

generated quantities {
  // The solution from the deterministic ODE solver.
  real y_forc[n_forecast, n_difeq];
  // The initial condition for the forward projection = the end state of the data.
  real y_init_forc[n_difeq];
  // For model-checking: simulated count data.
  real y_sim[n_obs];
  // For model-checking: simulated Poisson rates.
  real lambda_sim[n_obs];
  
  // Generate forward projections deterministically, conditional on the final value of the contact-rate parameter.
  y_init_forc[1] = y_hat[n_obs, 1];
  y_init_forc[2] = y_hat[n_obs, 2];
  y_init_forc[3] = y_hat[n_obs, 3];
  y_init_forc[4] = y_hat[n_obs, 4];
  y_init_forc[5] = y_hat[n_obs, 5];
  y_init_forc[6] = y_hat[n_obs, 6];
  y_forc = integrate_ode_rk45(SIR, y_init_forc, t0, tf,
                              {exp(-log_beta_t_deviation[n_obs])*beta, a * (1.0 - untested_share), gamma, a * untested_share},
                              x_r, x_i);
  
  if (fit_model == 0) {
    // If doing model-checking, generate simulated count data.
    for (i in 1:n_obs){
      // Public datasets report cumulative confirmed cases, which is I+R in this model.
      lambda_sim[i] = (y_hat[i,3] + y_hat[i,4]) * n_pop;
    }
    if (min(lambda_sim) < 0) {
      print({beta, a, gamma});
      print(sigma);
      print(log_beta_t_deviation);
      print(y_hat);
    }
    y_sim = poisson_rng(lambda_sim);
  } else {
    // Fill the simulated values with noise to prevent the model-fitting code throwing an error
    // due to NaN effective sample sizes.
    for (i in 1:n_obs) {
      lambda_sim[i] = uniform_rng(0.0, n_pop);
    }
    y_sim = poisson_rng(lambda_sim);
  }
}
