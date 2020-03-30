library(rstan)
library(tidyverse)
options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)

# https://www.abs.gov.au/AUSSTATS/abs@.nsf/Latestproducts/3101.0Main%20Features3Sep%202019?opendocument&tabname=Summary&prodno=3101.0&issue=Sep%202019&num=&view=
population <- 25464100

australia_cases <- read_csv("https://raw.githubusercontent.com/pappubahry/AU_COVID19/master/time_series_cases.csv")
# Growth rates are very different <100 cases, so start from that point.
obs_dataframe <- australia_cases %>% filter(Total >= 100)
obs <- obs_dataframe %>% pull(Total)
n_obs <- length(obs)

corona_data <- list(n_obs=n_obs,
                    n_theta=2,
                    n_difeq=3,
                    n_pop=population,
                    y=obs,
                    t0=0,
                    ts=1:n_obs)
initial_values <- function() {
  list(theta=c(runif(1,0,5), runif(1,0.2,0.4)), 
       S0=runif(1,(population-130)/population,(population-95)/population)
  )
}
# model <- vb(stan_model("sir.stan"), data=corona_data, init=initial_values)
model <- stan("sir.stan", data=corona_data, init=initial_values)
fitted_values <- rstan::extract(model)
fitted_I <- fitted_values$y_hat[,,2] * population
fitted_cases <- (fitted_values$y_hat[,,2] + fitted_values$y_hat[,,3]) * population
median_cases <- apply(fitted_cases, 2, median)
low_cases <-  apply(fitted_cases, 2, quantile, probs=c(0.025))
high_cases = apply(fitted_cases, 2, quantile, probs=c(0.975))

fitted_plot <- ggplot(data.frame(day=1:n_obs, obs, fitted_cases, median_cases, low_cases, high_cases)) + aes(x=day) +
  geom_ribbon(aes(ymin=low_cases, ymax=high_cases), fill = "light green", alpha = 0.6) +
  geom_line(aes(y=median_cases, color = "Fitted median"), size = 1.3) +
  geom_point(shape=19, size=1, (aes(y=obs, color="Data"))) +
  scale_y_log10()

actuals_plot <- ggplot(data.frame(day=obs_dataframe %>% pull(Date),
                                  median=apply(fitted_I, 2, median),
                                  low=apply(fitted_I, 2, quantile, probs=c(0.025)),
                                  hi=apply(fitted_I, 2, quantile, probs=c(0.975)))) +
  aes(x=day) + geom_ribbon(aes(ymin=low, ymax=hi), fill="light blue", alpha=0.5) +
  geom_line(aes(y=median), colour="blue", size=1.5) +
  scale_y_log10() + ylab("") + xlab("")





