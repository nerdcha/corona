# A Coronavirus forecasting model

This repo contains a Bayesian model that can be used to analyse and forecast the COVID-19 pandemic in a local area.

## Methodology

This is an [SEIR](https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SEIR_model) model with a time-varying contact rate (beta). The beta estimate is constrained to start near the prior distribution of R0, then is allowed to wander as lockdowns and restrictions are put in place. Priors are placed on the dynamic parameters based on prior literature; the unobserved beta component uses somewhat-informative priors which keep the model stable.

To generate charts for Australian data, run [corona.R](./corona.R), which uses the Stan model specified in [sir.stan](./sir.stan).

Comments, feedback, questions and suggestions are all welcome: feel free to open an Issue or a Pull Request here, ping me on [Twitter](https://twitter.com/jamie_hall), or email me an email (`jamie1212@gmail.com`).

## Current results

![Chart of projections for Australia](./pics/forecast.png?raw=true "Forecasts")

Here is the most recently-generated forecast for Australia. Note that standard caveats and modelling humility apply here. One way to understand this chart is "Taking the reported data at face value, and making reasonable assumptions about how the virus might behave, what might the future look like?"

## Contributing

This repo follows Google's standard open-source conditions; see the guidelines in `CONTRIBUTING.md` for more details. Please note: this is not an officially supported Google product.