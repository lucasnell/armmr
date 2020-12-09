#==========
#========== Preliminaries
#==========

library(tidyverse)
library(armmr)
library(rstan)

options(mc.cores = parallel::detectCores()-4)

data <- read_csv("data_fit.csv")





#==========
#========== Exmaple 1: variation by year in lycosids
#==========

# prepare data and examine
data_lyco <- data %>%
  filter(taxon == "lyco") %>%
  mutate(plot_f = as_factor(plot),
         trans_f = as_factor(trans),
         year_f = as_factor(year))
data_lyco %>%
  ggplot(aes(year, y, group = plot_f))+
  geom_line()

# fit model using `armmr`
# include the variable that changes within years as both a fixed effect and a random effect
mod <- armm(formula = y ~ year_f +
              (1 | plot_f) +
              (year_f | year_f),
            time_form = ~  time | plot_f,
            ar_form = ~ 1,
            obs_error = F,
            family = "ss",
            data = data_lyco,
            x_scale = FALSE,
            y_scale = NULL,
            hmc = T,
            change = T,
            rstan_control = list(iter = 2000, chains = 4, seed = 3e3,
                                 control = list(adapt_delta = 0.95)))

# extract stan_data
stan_data_timehack <- mod$stan_data

# define binary indicator for whether to include the fixed effect component
stan_data_timehack$coef_b <- c(1, rep(0, 9))

# define a vecotr mapping random effect to standard deviations
stan_data_timehack$sig_beta_map <- c(1, rep(2, 9))


# fit model using `rstan` and revised stan file
chains <- 4
iter <- 2000
adapt_delta <- 0.9
max_treedepth <- 10
fit_timehack <- stan(file = "armm_ss_timehack.stan", 
                     data = stan_data_timehack, 
                     seed=2e3,
            chains = chains, iter = iter,
            control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth))

# summarize
fit_th_sum <- rstan::summary(fit_timehack, probs = c(0.16,0.50,0.84))$summary %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    as_tibble() %>%
    rename(var = rowname,
           lo = `16%`,
           mi = `50%`,
           hi = `84%`) %>%
    select(var, lo, mi, hi, n_eff, Rhat)

# examine random effects standard deviations
fit_th_sum %>%
  filter(str_detect(fit_th_sum$var, "sig_beta")) 

# examine fixed effects
fit_th_sum %>%
  filter(str_detect(fit_th_sum$var, "alpha")) 

# pull out year effect
trial <- fit_th_sum %>%
  filter(str_detect(fit_th_sum$var, "beta"), 
         !str_detect(fit_th_sum$var, "sig")) %>%
  mutate(id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         coef = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3]))) %>%
  select(-var, -id) %>%
  filter(coef > 1) %>%
  unique() %>%
  mutate(year = 2009:2017)

# plot
data_lyco %>%
  ggplot(aes(year, y, group = plot_f))+
  geom_line(size = 0.2, alpha = 0.5)+
  geom_ribbon(data = trial,
            aes(x = year,
                ymin = lo,
                ymax = hi),
            alpha = 0.2,
            fill = "dodgerblue",
            inherit.aes = F)+
  geom_line(data = trial,
            aes(x = year,
                y = mi),
            size = 1.5,
            color = "dodgerblue",
            inherit.aes = F)+
  theme_bw()





#==========
#========== Exmaple 2: variation in year x taxon (direct parameterization)
#==========

# prepare data and examines
data_2<- data %>%
  filter(taxon %in% c("lyco","gnap")) %>%
  mutate(plot_f = as_factor(plot),
         trans_f = as_factor(trans),
         taxon_f = as_factor(taxon),
         year_taxon = as_factor(paste(taxon, year, sep = "_")))
data_2 %>%
  ggplot(aes(year, y, group = plot_f))+
  facet_wrap(~taxon)+
  geom_line()
data_2 %>%
  ggplot(aes(year, y, color = year_taxon))+
  facet_wrap(~taxon)+
  geom_point()

# fit model using `armmr`
# include the variable that changes within years as both a fixed effect and a random effect
mod <- armm(formula = y ~ taxon + year_taxon +
              (1 | plot_f) +
              (year_taxon | year_taxon),
            time_form = ~  time | plot_f + taxon_f,
            ar_form = ~ taxon_f,
            obs_error = F,
            family = "ss",
            data = data_2,
            x_scale = FALSE,
            y_scale = NULL,
            hmc = T,
            change = T,
            rstan_control = list(iter = 1, chains = 1, seed = 3e3,
                                 control = list(adapt_delta = 0.95)))

# extract stan_data
stan_data_timehack <- mod$stan_data

# define binary indicator for whether to include the fixed effect component
stan_data_timehack$coef_b <- c(1, 1, rep(0, 19))

# define a vecotr mapping random effect to standard deviations
stan_data_timehack$sig_beta_map <- c(1, 
                                     rep(2, 9),
                                     rep(3, 10))

# fit model using `rstan` and revised stan file
chains <- 4
iter <- 2000
adapt_delta <- 0.9
max_treedepth <- 10
fit_timehack <- stan(file = "armm_ss_timehack.stan",
                     data = stan_data_timehack,
                     seed=2e3,
                     chains = chains, iter = iter,
                     control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth))

# summarize
fit_th_sum <- rstan::summary(fit_timehack, probs = c(0.16,0.50,0.84))$summary %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  rename(var = rowname,
         lo = `16%`,
         mi = `50%`,
         hi = `84%`) %>%
  select(var, lo, mi, hi, n_eff, Rhat)

# examine random effects standard deviations
fit_th_sum %>%
  filter(str_detect(fit_th_sum$var, "sig_beta")) 

# examine fixed effects
fit_th_sum %>%
  filter(str_detect(fit_th_sum$var, "alpha")) 

# pull out predicted values (summarize by year and taxon)
trial <- fit_th_sum %>%
  filter(str_detect(fit_th_sum$var, "y_pred")) %>%
  mutate(id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  full_join(data_2 %>%
              mutate(id = row_number()) %>%
              select(id, taxon, year)) %>%
  group_by(taxon, year) %>%
  summarize(mi = mean(mi))
data_2 %>%
  ggplot(aes(year, y, group = plot_f))+
  facet_wrap(~taxon)+
  geom_line(size = 0.2, alpha = 0.5)+
  geom_line(data = trial,
            aes(x = year,
                y = mi),
            size = 1.5,
            color = "dodgerblue",
            inherit.aes = F)+
  theme_bw()

# year x taxon effect 
coef <- fit_th_sum %>%
  filter(str_detect(fit_th_sum$var, "beta"),
         !str_detect(fit_th_sum$var, "sig")) %>%
  mutate(id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         coef = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])))  %>%
  filter(coef > 2) %>%
  select(-var, -id, -Rhat, -n_eff) %>%
  group_by(coef) %>%
  summarize(mi = mean(mi)) %>%
  mutate(taxon = ifelse(coef < 12, "gnap","lyco")) %>%
  group_by(taxon) %>%
  mutate(year = row_number())
coef %>% 
  ggplot(aes(year, mi))+
  facet_wrap(~taxon)+
  geom_line()






#==========
#========== Exmaple 3: variation in year x taxon (contrasts parameterization)
#==========

# prepare data and examines
data_2<- data %>%
  filter(taxon %in% c("lyco","gnap")) %>%
  mutate(plot_f = as_factor(plot),
         trans_f = as_factor(trans),
         taxon_f = as_factor(taxon),
         dummy = 1,
         year_f = as_factor(year),
         year_taxon = ifelse(taxon == "lyco", "null", year_f),
         year_taxon = as_factor(year_taxon))
data_2 %>%
  ggplot(aes(year, y, group = plot_f))+
  facet_wrap(~taxon)+
  geom_line()
data_2 %>%
  ggplot(aes(year, y, color = year_taxon))+
  facet_wrap(~taxon)+
  geom_point()

# fit model using `armmr`
# include the variable that changes within years as both a fixed effect and a random effect
mod <- armm(formula = y ~ taxon + year_f + year_taxon +
              (1 | plot_f) +
              (year_f | year_f) +
              (year_taxon | year_taxon),
            time_form = ~  time | plot_f + taxon_f,
            ar_form = ~ taxon_f,
            obs_error = F,
            family = "ss",
            data = data_2,
            x_scale = FALSE,
            y_scale = NULL,
            hmc = T,
            change = T,
            rstan_control = list(iter = 1, chains = 1, seed = 3e3,
                                 control = list(adapt_delta = 0.95)))

# extract stan_data
stan_data_timehack <- mod$stan_data

# define binary indicator for whether to include the fixed effect component
stan_data_timehack$coef_b <- c(1, 1, rep(0, 19))

# define a vecotr mapping random effect to standard deviations
stan_data_timehack$sig_beta_map <- c(1, 
                                     rep(2, 9),
                                     rep(3, 10))

# fit model using `rstan` and revised stan file
chains <- 4
iter <- 2000
adapt_delta <- 0.9
max_treedepth <- 10
fit_timehack <- stan(file = "armm_ss_timehack.stan",
                     data = stan_data_timehack,
                     seed=2e3,
                     chains = chains, iter = iter,
                     control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth))

# summarize
fit_th_sum <- rstan::summary(fit_timehack, probs = c(0.16,0.50,0.84))$summary %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  rename(var = rowname,
         lo = `16%`,
         mi = `50%`,
         hi = `84%`) %>%
  select(var, lo, mi, hi, n_eff, Rhat)

# examine random effects standard deviations
fit_th_sum %>%
  filter(str_detect(fit_th_sum$var, "sig_beta")) 

# examine fixed effects
fit_th_sum %>%
  filter(str_detect(fit_th_sum$var, "alpha")) 

# pull out predicted values (summarize by year and taxon)
trial <- fit_th_sum %>%
  filter(str_detect(fit_th_sum$var, "y_pred")) %>%
  mutate(id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  full_join(data_2 %>%
              mutate(id = row_number()) %>%
              select(id, taxon, year)) %>%
  group_by(taxon, year) %>%
  summarize(mi = mean(mi))
data_2 %>%
  ggplot(aes(year, y, group = plot_f))+
  facet_wrap(~taxon)+
  geom_line(size = 0.2, alpha = 0.5)+
  geom_line(data = trial,
            aes(x = year,
                y = mi),
            size = 1.5,
            color = "dodgerblue",
            inherit.aes = F)+
  theme_bw()

# year x taxon effect 
coef <- fit_th_sum %>%
  filter(str_detect(fit_th_sum$var, "beta"),
         !str_detect(fit_th_sum$var, "sig")) %>%
  mutate(id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         coef = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])))  %>%
  filter(coef > 2) %>%
  select(-var, -id, -Rhat, -n_eff) %>%
  group_by(coef) %>%
  summarize(mi = mean(mi)) %>%
  mutate(taxon = ifelse(coef < 12, "lyco","gnap")) %>%
  group_by(taxon) %>%
  mutate(year = row_number())
coef %>% 
  ggplot(aes(year, mi))+
  facet_wrap(~taxon)+
  geom_line()

