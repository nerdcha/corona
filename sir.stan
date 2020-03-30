// Based on https://github.com/anastasiachtz/COMMAND_stan/blob/master/SingleStrainStan.Rmd
// documented in https://arxiv.org/pdf/1903.00423.pdf

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
  int<lower=0> n_forecast;   // Number of days to project.
  int<lower=1> n_theta;     // Number of model parameters.
  int<lower=1> n_difeq;     // Number of differential equations.
  int<lower=1> n_pop;       // Population.
  int y[n_obs];           // Data: total number of infected individuals each day
  real t0;      // Initial time tick for the ODE solver. Must be provided in the data block.
  real ts[n_obs];  // Time ticks for the ODE solver. Must be provided in the data block.
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
  real<lower=0> I0;      // initial number of infected individuals
  real<lower=0> E0;      // Initial number of exposed individuals
}
  
transformed parameters{
  real y_hat[n_obs, n_difeq]; // solution from the ODE solver
  real y_init[n_difeq];       // initial conditions for both fractions of S and I

  y_init[1] = 1.0 - I0/n_pop - E0/n_pop;
  y_init[2] = E0/n_pop;
  y_init[3] = I0/n_pop;
  y_init[4] = 0;
  y_hat = integrate_ode_rk45(SIR, y_init, t0, ts, theta, x_r, x_i);
}
  
model {
  real lambda[n_obs];              // Poisson rate parameter at each time point.
  
  // Priors.
  theta[1] ~ lognormal(0, 1);
  theta[2] ~ gamma(20.0, 20.0/5.0);  // Mean asymptomatic period = 5 days.
  theta[3] ~ gamma(20.0, 20.0/7.0);  // Mean infectious period = 7 days.
  E0 ~ normal(300, 100);
  I0 ~ normal(100, 2);
  
  // Likelihood
  for (i in 1:n_obs){
    // Public datasets report cumulative confirmed cases, which is I+R in this model.
    lambda[i] = (y_hat[i,3] + y_hat[i,4]) * n_pop;
  }
  y ~ poisson(lambda);
}
  
generated quantities {
  real y_forc[n_forecast, n_difeq]; // solution from the ODE solver
  real y_init_forc[n_difeq];  // Initial condition for the forward projection.
  real R0;
  
  R0 = theta[1]/theta[3];
  
  y_init_forc[1] = y_hat[n_obs, 1];
  y_init_forc[2] = y_hat[n_obs, 2];
  y_init_forc[3] = y_hat[n_obs, 3];
  y_init_forc[4] = y_hat[n_obs, 4];
  y_forc = integrate_ode_rk45(SIR, y_init_forc, t0, tf, theta, x_r, x_i);
}
