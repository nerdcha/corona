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
n_forecast <- 150

corona_data <- list(n_obs=n_obs,
                    n_forecast=n_forecast,
                    n_theta=3,
                    n_difeq=4,
                    n_pop=population,
                    y=obs,
                    t0=0,
                    ts=1:2,
                    tf=1:n_forecast)
initial_values <- function() {
  list(theta=c(runif(1, 0, 5), runif(1, 0.2, 0.4), runif(1, 0.2, 0.4)), 
       I0=runif(1, 90, 150),
       E0=runif(1, 50, 1000)
  )
}
# model <- vb(stan_model("sir.stan"), data=corona_data, init=initial_values)
model <- stan("sir.stan", data=corona_data, init=initial_values, chains=2)
fitted_values <- rstan::extract(model)
fitted_I <- fitted_values$y_hat[,,3] * population
fitted_cases <- (fitted_values$y_hat[,,3] + fitted_values$y_hat[,,4]) * population
median_cases <- apply(fitted_cases, 2, median)
low_cases <-  apply(fitted_cases, 2, quantile, probs=c(0.025))
high_cases = apply(fitted_cases, 2, quantile, probs=c(0.975))

fitted_plot <- ggplot(data.frame(day=1:n_obs, obs, median_cases, low_cases, high_cases)) + aes(x=day) +
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

forecast_plot_data <- data.frame(day=seq(from=obs_dataframe$Date[1], length.out=n_obs + n_forecast, by="day"),
                                 exposed=c(apply(fitted_values$y_hat[,,2] * population, 2, median), rep(NA, n_forecast)),
                                 fitted=c(apply(fitted_values$y_hat[,,3] * population, 2, median), rep(NA, n_forecast)),
                                 projected=c(rep(NA, n_obs), apply(fitted_values$y_forc[,,3]*population, 2, median)),
                                 hi=c(rep(NA, n_obs), apply(fitted_values$y_forc[,,3]*population, 2, quantile, probs=c(0.975))),
                                 lo=c(rep(NA, n_obs), apply(fitted_values$y_forc[,,3]*population, 2, quantile, probs=c(0.025))))
forecast_plot <- ggplot(forecast_plot_data) + aes(x=day) +
  geom_ribbon(aes(ymin=lo, ymax=hi), fill="purple", alpha=0.5) +
  geom_line(aes(y=exposed), colour="red") +
  geom_line(aes(y=fitted), colour="blue") + geom_line(aes(y=projected), colour="purple") +
  scale_y_log10()

forecast_cases_plot_data <- data.frame(day=seq(from=obs_dataframe$Date[1], length.out=n_obs + n_forecast, by="day"),
                                       obs=c(obs, rep(NA, n_forecast)),
                                       fitted=c(median_cases, rep(NA, n_forecast)),
                                       projected=c(rep(NA, n_obs), apply((fitted_values$y_forc[,,3] + fitted_values$y_forc[,,4])*population, 2, median)),
                                       hi=c(rep(NA, n_obs), apply((fitted_values$y_forc[,,3] + fitted_values$y_forc[,,4])*population, 2, quantile, probs=c(0.975))),
                                       lo=c(rep(NA, n_obs), apply((fitted_values$y_forc[,,3] + fitted_values$y_forc[,,4])*population, 2, quantile, probs=c(0.025))))
forecast_cases_plot <- ggplot(forecast_cases_plot_data) + aes(x=day) +
  geom_ribbon(aes(ymin=lo, ymax=hi), fill="purple", alpha=0.5) +
  geom_line(aes(y=fitted), colour="blue") + geom_line(aes(y=projected), colour="purple") +
  geom_point(shape=19, size=1, aes(y=obs)) +
  scale_y_log10()

print(forecast_cases_plot)
print(forecast_plot)

