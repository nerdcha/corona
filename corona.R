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


# Fits an SEIR model to Australian observations of Covid-19.

library(rstan)
library(tidyverse)
library(ggthemes)
library(scales)

options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)

source("corona_settings.R")

output_dir <- paste(getwd(), "pics", sep="/")

# https://www.abs.gov.au/AUSSTATS/abs@.nsf/Latestproducts/3101.0Main%20Features3Sep%202019?opendocument&tabname=Summary&prodno=3101.0&issue=Sep%202019&num=&view=
population <- 25464100

australia_cases <- read_csv("https://raw.githubusercontent.com/pappubahry/AU_COVID19/master/time_series_cases.csv")
# Growth rates are very different <100 cases, so start from that point. (That insight is the basis of
# the Cowgill Chart; see https://blog.grattan.edu.au/2020/03/australian-governments-can-choose-to-slow-the-spread-of-coronavirus-but-they-must-act-immediately/)
obs_dataframe <- australia_cases %>% filter(Total >= 100)
obs <- obs_dataframe %>% pull(Total)
n_obs <- length(obs)
n_forecast <- 90

corona_data <- list(n_obs=n_obs,
                    n_forecast=n_forecast,
                    n_theta=3,
                    n_difeq=4,
                    n_pop=population,
                    y=obs,
                    t0=0,
                    ts=1:2,
                    tf=1:n_forecast,
                    fit_model=1)


model <- stan("sir.stan", data=corona_data, init=initial_values_function(n_obs),
              chains=n_chains, iter=n_warmup_draws + n_posterior_draws_per_chain,
              warmup=n_warmup_draws)
fitted_values <- rstan::extract(model)

# Given a TxD array and a string name, return a data frame of T-length vectors suitable for a fanchart.
# Truncates negative values to a small positive number.
get_fanchart_series <- function(values, name) {
  result <- list()
  small_positive_number <- 1
  result[[paste0(name, "_median")]] <- apply(values, 2, median)
  result[[paste0(name, "_lo5")]] <- apply(values, 2, quantile, probs=c(0.025)) %>% ifelse(. < 0, small_positive_number, .)
  result[[paste0(name, "_lo50")]] <- apply(values, 2, quantile, probs=c(0.25)) %>% ifelse(. < 0, small_positive_number, .)
  result[[paste0(name, "_hi50")]] <- apply(values, 2, quantile, probs=c(0.75)) %>% ifelse(. < 0, small_positive_number, .)
  result[[paste0(name, "_hi5")]] <- apply(values, 2, quantile, probs=c(0.975)) %>% ifelse(. < 0, small_positive_number, .)
  return(as.data.frame(result))
}


# Generate a plot of datapoints and in-sample fit.
fitted_plot_df <- obs_dataframe %>% select(Date, Total) %>%
  cbind(get_fanchart_series(population*(fitted_values$y_hat[,,3] + fitted_values$y_hat[,,4]), "fitted"))
fitted_plot <- ggplot(fitted_plot_df) + aes(x=Date) +
  geom_ribbon(aes(ymin=fitted_lo5, ymax=fitted_hi5), fill="lightblue", alpha=0.3) +
  geom_ribbon(aes(ymin=fitted_lo50, ymax=fitted_hi50), fill="lightblue", alpha=0.3) +
  geom_line(aes(y=fitted_median), colour="darkblue") +
  geom_point(aes(y=Total), colour="red", size=1.5) +
  scale_y_log10() +
  ggtitle("Australia—Covid19 cases", subtitle="Datapoints and in-sample model fit") +
  xlab("") + ylab("") +
  theme_solarized() + theme(
    title=element_text(face="bold", colour="gray20"),
    plot.subtitle=element_text(face="plain", colour="gray50"))
ggsave(fitted_plot, filename=paste(output_dir, "fitted.png", sep="/"), width=5, height=4)

# Generate a plot of Exposed and Infected out to forecast end, dashed line for forecast start.
fitted_exposed <- get_fanchart_series(population*(fitted_values$y_hat[,,2]), "exposed")
# Replace the Exposed projections with NA, to make the chart less confusing.
out_of_sample_exposed <- matrix(NA, n_forecast, ncol(fitted_exposed)) %>% as.data.frame()
names(out_of_sample_exposed) <- names(fitted_exposed)
forecast_plot_df <- fitted_exposed %>% rbind(out_of_sample_exposed) %>%
  cbind(get_fanchart_series(population*(cbind(fitted_values$y_hat[,,3], fitted_values$y_forc[,,3])), "infectious")) %>%
  cbind(data.frame(date=seq(from=obs_dataframe$Date[1], length.out=(n_obs + n_forecast), by="day")))
forecast_plot <- ggplot(forecast_plot_df) + aes(x=date) +
  geom_ribbon(aes(ymin=exposed_lo5, ymax=exposed_hi5), fill="red", alpha=0.1) +
  geom_ribbon(aes(ymin=exposed_lo50, ymax=exposed_hi50), fill="red", alpha=0.1) +
  geom_line(aes(y=exposed_median), colour="red") +
  geom_ribbon(aes(ymin=infectious_lo5, ymax=infectious_hi5), fill="blue", alpha=0.3) +
  geom_ribbon(aes(ymin=infectious_lo50, ymax=infectious_hi50), fill="blue", alpha=0.3) +
  geom_line(aes(y=infectious_median), colour="blue") +
  geom_vline(xintercept=obs_dataframe$Date[nrow(obs_dataframe)], linetype=2) +
  scale_y_log10(labels=comma) + coord_cartesian(ylim=c(50, 5e6)) +
  ggtitle("Covid19 in Australia—Projections", subtitle="Red: incubating; Blue: infectious") +
  labs(caption="Forecasts assume shutdowns stay at current level.") +
  xlab("") + ylab("") +
  theme_solarized() + theme(
    title=element_text(face="bold", colour="gray20"),
    plot.subtitle=element_text(face="plain", colour="gray50"),
    plot.caption=element_text(face="italic", colour="gray20"))
ggsave(forecast_plot, filename=paste(output_dir, "forecast.png", sep="/"), width=5, height=4)
write_csv(forecast_plot_df, path=paste(output_dir, "forecast.csv", sep="/"))

# Generate a plot of contact-rate over time.
contact_rate_df <- obs_dataframe %>% select(Date, Total) %>%
  cbind(get_fanchart_series(exp(fitted_values$log_beta_t_deviation) * fitted_values$theta[,1]
                            / fitted_values$theta[,3], "R_t"))
contact_rate_plot <- ggplot(contact_rate_df) + aes(x=Date) +
  annotate(geom="rect", ymin=0, ymax=1,
           xmin=contact_rate_df$Date[1] - 1,
           xmax=contact_rate_df$Date[n_obs] + 1,
           fill="royalblue1", alpha = 0.1) + 
  geom_ribbon(aes(ymin=R_t_lo5, ymax=R_t_hi5), fill="tomato3", alpha=0.3) +
  geom_ribbon(aes(ymin=R_t_lo50, ymax=R_t_hi50), fill="tomato3", alpha=0.3) +
  geom_line(aes(y=R_t_median), colour="tomato1", size=1.5) +
  ggtitle("Australia—Covid19 infection rate") +
  xlab("") + ylab("") +
  theme_solarized() + theme(
    title=element_text(face="bold", colour="gray20"),
    plot.subtitle=element_text(face="plain", colour="gray50"),
    plot.caption=element_text(face="italic", colour="gray20"))
ggsave(contact_rate_plot, filename=paste(output_dir, "contact_rate.png", sep="/"), width=5, height=4)

