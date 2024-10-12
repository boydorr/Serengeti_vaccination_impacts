rm(list=ls())

library(rstan) 
library(loo)
library(brms)


## Load data
#______________________

data_vill <- read.csv("Output/incidence_coverage_model_data_village.csv")
data <- readRDS("Output/power_mean_model_village_data.rds")
colnames(data$X)

# Variations of the data for different models
data_standardised <- data
data_standardised$X[,2:ncol(data_standardised$X)] <- scale(data_standardised$X[,2:ncol(data_standardised$X)])
data_wo_distant_cases <- data
data_wo_distant_cases$X <- data$X[,which(!colnames(data$X) %in% c("log_case_rate_neighbours_last2monthMean","log_case_rate_notNeighbours_last2monthMean"))]
data_wo_distant_cases$K <- data_wo_distant_cases$K-2
data_wo_distant_cases_standardised <- data_wo_distant_cases
data_wo_distant_cases_standardised$X[,2:ncol(data_wo_distant_cases_standardised$X)] <- scale(data_wo_distant_cases_standardised$X[,2:ncol(data_wo_distant_cases_standardised$X)])
data_wo_cases <- data
data_wo_cases$X <- data$X[,which(!colnames(data$X) %in% c("log_case_rate_last2monthMean","log_case_rate_neighbours_last2monthMean","log_case_rate_notNeighbours_last2monthMean"))]
data_wo_cases$K <- data_wo_cases$K-3
data_wo_cases_standardised <- data_wo_cases
data_wo_cases_standardised$X[,2:ncol(data_wo_cases_standardised$X)] <- scale(data_wo_cases_standardised$X[,2:ncol(data_wo_cases_standardised$X)])




# Fit models
#______________________

# Full model
model_vill_month_stan <- stan("stan/power_mean_model_village.stan", data = data, iter = 3000,
                              warmup = 1500, thin = 1, chains = 4, verbose = TRUE, cores=4, seed=5,
                              control = list(adapt_delta = 0.95, max_treedepth = 10),include = TRUE, 
                              pars = c("beta","p","phi","Intercept","sigma_village","gamma_t"))
plot(model_vill_month_stan,pars=c("Intercept","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","beta[7]","beta[8]","beta[9]","p","phi"))
traceplot(model_vill_month_stan,inc_warmup=F,pars=c("Intercept","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","beta[7]","beta[8]","beta[9]","p","phi"))
summary(model_vill_month_stan,pars=c("Intercept","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","beta[7]","beta[8]","beta[9]","p","phi"),probs = c(0.025,0.5, 0.975))$summary
saveRDS(model_vill_month_stan,"output/stan_models/incidence_from_vax_model_village_stan.rds")

# Without distant cases
model_vill_month_stan_wo_distant_cases <- stan("stan/power_mean_model_village.stan", data = data_wo_distant_cases, iter = 3000,
                                               warmup = 1500, thin = 1, chains = 4, verbose = TRUE, cores=4, seed=5,
                                               control = list(adapt_delta = 0.95, max_treedepth = 10),include = TRUE, 
                                               pars = c("beta","p","phi","Intercept","sigma_village","gamma_t"))
plot(model_vill_month_stan_wo_distant_cases,pars=c("Intercept","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","beta[7]","p","phi"))
traceplot(model_vill_month_stan_wo_distant_cases,inc_warmup=F,pars=c("Intercept","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","beta[7]","p","phi"))
summary(model_vill_month_stan_wo_distant_cases,pars=c("Intercept","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","beta[7]","p","phi"),probs = c(0.025,0.5, 0.975))$summary
saveRDS(model_vill_month_stan_wo_distant_cases,"output/stan_models/incidence_from_vax_model_village_stan_wo_distant_cases.rds")

# Without cases
model_vill_month_stan_wo_cases <- stan("stan/power_mean_model_village.stan", data = data_wo_cases, iter = 3000,
                                       warmup = 1500, thin = 1, chains = 4, verbose = TRUE, cores=4, seed=5,
                                       control = list(adapt_delta = 0.95, max_treedepth = 10),include = TRUE, 
                                       pars = c("beta","p","phi","Intercept","sigma_village","gamma_t"))
plot(model_vill_month_stan_wo_cases,pars=c("Intercept","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","p","phi"))
traceplot(model_vill_month_stan_wo_cases,inc_warmup=F,pars=c("Intercept","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","p","phi"))
summary(model_vill_month_stan_wo_cases,pars=c("Intercept","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","p","phi"),probs = c(0.025,0.5, 0.975))$summary
saveRDS(model_vill_month_stan_wo_cases,"output/stan_models/incidence_from_vax_model_village_stan_wo_cases.rds")


# Standardised
#-----------

# Full model
model_vill_month_stan_standardised <- stan("stan/power_mean_model_village_standardised.stan", data = data_standardised, iter = 3000,
                                           warmup = 1500, thin = 1, chains = 4, verbose = TRUE, cores=4, seed=5,
                                           control = list(adapt_delta = 0.95, max_treedepth = 10),include = TRUE, 
                                           pars = c("beta","p","phi","Intercept","sigma_village","gamma_t"))
plot(model_vill_month_stan_standardised,pars=c("Intercept","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","beta[7]","beta[8]","beta[9]","p","phi"))
traceplot(model_vill_month_stan_standardised,inc_warmup=F,pars=c("Intercept","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","beta[7]","beta[8]","beta[9]","p","phi"))
summary(model_vill_month_stan_standardised,pars=c("Intercept","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","beta[7]","beta[8]","beta[9]","p","phi"),probs = c(0.025,0.5, 0.975))$summary
saveRDS(model_vill_month_stan_standardised,"output/stan_models/incidence_from_vax_model_village_stan_standardised.rds")

# Without distant cases
model_vill_month_stan_wo_distant_cases_standardised <- stan("stan/power_mean_model_village_standardised.stan", data = data_wo_distant_cases_standardised, iter = 3000,
                                                            warmup = 1500, thin = 1, chains = 4, verbose = TRUE, cores=4, seed=5,
                                                            control = list(adapt_delta = 0.95, max_treedepth = 10),include = TRUE, 
                                                            pars = c("beta","p","phi","Intercept","sigma_village","gamma_t"))
plot(model_vill_month_stan_wo_distant_cases_standardised,pars=c("Intercept","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","beta[7]","p","phi"))
traceplot(model_vill_month_stan_wo_distant_cases_standardised,inc_warmup=F,pars=c("Intercept","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","beta[7]","p","phi"))
summary(model_vill_month_stan_wo_distant_cases_standardised,pars=c("Intercept","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","beta[7]","p","phi"),probs = c(0.025,0.5, 0.975))$summary
saveRDS(model_vill_month_stan_wo_distant_cases_standardised,"output/stan_models/incidence_from_vax_model_village_stan_wo_distant_cases_standardised.rds")

# Without cases
model_vill_month_stan_wo_cases_standardised <- stan("stan/power_mean_model_village_standardised.stan", data = data_wo_cases_standardised, iter = 3000,
                                                    warmup = 1500, thin = 1, chains = 4, verbose = TRUE, cores=4, seed=5,
                                                    control = list(adapt_delta = 0.95, max_treedepth = 10),include = TRUE, 
                                                    pars = c("beta","p","phi","Intercept","sigma_village","gamma_t"))
plot(model_vill_month_stan_wo_cases_standardised,pars=c("Intercept","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","p","phi"))
traceplot(model_vill_month_stan_wo_cases_standardised,inc_warmup=F,pars=c("Intercept","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","p","phi"))
summary(model_vill_month_stan_wo_cases_standardised,pars=c("Intercept","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","p","phi"),probs = c(0.025,0.5, 0.975))$summary
saveRDS(model_vill_month_stan_wo_cases_standardised,"output/stan_models/incidence_from_vax_model_village_stan_wo_cases_standardised.rds")



