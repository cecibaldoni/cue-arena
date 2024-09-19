library(brms)
library(tidyverse)
library(bayesplot)
library(rstan)
library(tidybayes)
library(ggdist)
library(bayestestR)
library(rstanarm)
library(cmdstanr)

# Load & Libraries ####

path_data <- read.csv("~/data/cue_learning/distance/master_distance.csv") %>% 
  mutate(unique_trial_ID = as.factor(paste(season, ID, trial, sep = "_")))

path_data <- path_data %>%
  left_join(speed_summary, by = "unique_trial_ID")

path_data <- path_data %>% 
  mutate(trial = str_remove(trial, "^T"),
         season = as.factor(season),
         ID = as.factor(ID)) %>% 
  mutate(season_status = case_when(season == "summer" ~ "summer_wild",
                            season == "winter" & grepl("^2021", ID) ~ "winter_wild",
                            season == "winter" & grepl("^2020", ID) ~ "winter_captive",
                            season == "winter" & grepl("^2024", ID) ~ "winter_wild",
                            season == "spring" & grepl("^2024", ID) ~ "spring_wild",
                            season == "spring" & grepl("^2020", ID) ~ "spring_captive"),
         status = case_when(season == "summer" ~ "wild",
                            season =="winter" & grepl("^2021", ID) ~ "wild",
                            season =="winter" & grepl("^2020", ID) ~ "captive",
                            season == "winter" & grepl("^2024", ID) ~ "wild",
                            season == "spring" & grepl("^2024", ID) ~ "wild",
                            season == "spring" & grepl("^2020", ID) ~ "captive"),
    status = as.factor(status),
    season_status = as.factor(season_status),
    ID = as.factor(ID)) %>%
  filter(speed_trip_to <40)

path_data$trial <- factor(path_data$trial, levels = 1:10)
path_data$season <- factor(path_data$season, levels = c("summer", "winter", "spring"))
path_data$ID <- as.factor(path_data$ID)

#EDA ####
summary(path_data)

sapply(path_data, function(x) sum(is.na(x)))

library(corrplot)
cor_matrix <- cor(path_data[, sapply(path_data, is.numeric)], use = "complete.obs")
corrplot(cor_matrix, method = "circle")

ggplot(path_data, aes(x = status, y = smooth_straightness_to)) +
  geom_boxplot() +
  labs(title = "Walked To by Status")

ggplot(path_data, aes(x = smooth_straightness_to)) +
  geom_density(fill = "blue", alpha = 0.4) +
  facet_wrap(~status)

ggplot(path_data) +
  geom_jitter(mapping = aes(x = trial, y = smooth_straightness_to, color = status, group = status), alpha =0.3) +
  geom_smooth(mapping = aes(x = trial, y = smooth_straightness_to, color = status, group = status)) +
  facet_wrap(~status) +
  theme_bw()

ggplot(path_data) +
  geom_jitter(mapping = aes(x = trial, y = speed_trip_to, color = status, group = status), alpha =0.3) +
  geom_smooth(mapping = aes(x = trial, y = speed_trip_to, color = status, group = status), method = "lm") +
  facet_wrap(~status) +
  theme_bw()
ggplot(path_data) +
  geom_jitter(mapping = aes(x = trial, y = speed_trip_back, color = status, group = status), alpha =0.3) +
  geom_smooth(mapping = aes(x = trial, y = speed_trip_back, color = status, group = status), method = "lm") +
  facet_wrap(~status) +
  theme_bw()

ggplot(path_data) +
  geom_jitter(mapping = aes(x = trial, y = walked_to, color = status, group = status), alpha =0.3) +
  geom_smooth(mapping = aes(x = trial, y = walked_to, color = status, group = status)) +
  facet_wrap(~status, scales = "free_y") +
  theme_bw()

ggplot(path_data) +
  geom_jitter(mapping = aes(x = trial, y = smooth_straightness_back, color = status, group = status), alpha =0.3) +
  geom_smooth(mapping = aes(x = trial, y = smooth_straightness_back, color = status, group = status)) +
  facet_wrap(~status) +
  theme_bw()
ggplot(path_data) +
  geom_jitter(mapping = aes(x = trial, y = walked_back, color = status, group = status), alpha =0.3) +
  geom_smooth(mapping = aes(x = trial, y = walked_back, color = status, group = status)) +
  facet_wrap(~status, scales = "free_y") +
  theme_bw()

ggplot(path_data) +
  geom_jitter(mapping = aes(x = speed_trip_to, y = smooth_straightness_to, color = status, group = status)) +
  geom_smooth(mapping = aes(x = speed_trip_to, y = smooth_straightness_to, color = status, group = status)) +
  theme_bw()
ggplot(path_data) +
  geom_point(mapping = aes(x = trial, y = smooth_straightness_to, size = speed_trip_to, color = status), alpha = 0.6) +
  geom_smooth(mapping = aes(x = trial, y = smooth_straightness_to, color = status), method = "lm") +
  scale_size_continuous(range = c(1, 10)) +  # Adjust size scale for visibility
  facet_wrap(~status) +  # Facet by status
  theme_bw() +
  labs(size = "Speed (cm/s)", x = "Trial", y = "Smooth Straightness To", title = "Smooth Straightness vs. Trial with Speed by Status")



## Straightness to ####
path_data$trial_n <- as.numeric(path_data$trial)
min(path_data$smooth_straightness_to)
max(path_data$smooth_straightness_to)

## if max(path_data$smooth_straightness_to) ==1, run the next lines to put it between (0,1)
# N <- nrow(path_data)
# s <- 0.5
# transformed_value <- 1 * (N - 1 + s) / N
# 
# path_data$smooth_straightness_to <- ifelse(path_data$smooth_straightness_to == 1,
#                                            transformed_value,
#                                            path_data$smooth_straightness_to)


## Models ####

#start simple:

get_prior(bf(smooth_straightness_to ~ s(trial_n) + status + season_status + speed_trip_to + season + (1|ID)),
          data = path_data, family = Beta())
model1 <- brm(bf(smooth_straightness_to ~ s(trial_n) + status + season_status + speed_trip_to + season + (1|ID)),
                    data = path_data, family = Beta(),
                    prior = c(prior(normal(0,1), class = "b"),
                              prior(exponential(1), class = "sd")),
                    chains = 4, iter = 4000, cores = 4, seed = 1234,
                    backend = "cmdstanr", threads = threading(2),
                    control = list(adapt_delta = 0.99, max_treedepth = 15))
p1 = pp_check(model1, type = 'stat', stat = 'min') +
  ggtitle('Prior predictive distribution of minimum values') 
p2 = pp_check(model1, type = 'stat', stat = 'max') +
  ggtitle('Prior predictive distribution of maximum values')  +
  #geom_vline(xintercept = maxvalue) +
  theme_classic()
p3 = pp_check(model1, type = 'stat', stat = 'mean') +
  ggtitle('Prior predictive distribution of average values')   +
  #geom_vline(xintercept = meanvalue) +
  theme_classic()
p4 = pp_check(model1, ndraws = 100)
ggpubr::ggarrange(p1,p2,p3, p4)
conditional_effects(model1)

model2 <- brm(formula = smooth_straightness_to ~ s(trial_n) + status * speed_trip_to + season_status + season + (1 | ID),
              data = path_data, family = Beta(),
              prior = c(prior(normal(0, 1), class = "b"), 
                        prior(exponential(1), class = "sd")),
              chains = 4, iter = 4000, cores = 4, seed = 1234)
p1 = pp_check(model2, type = 'stat', stat = 'min') +
  ggtitle('Prior predictive distribution of minimum values') 
p2 = pp_check(model2, type = 'stat', stat = 'max') +
  ggtitle('Prior predictive distribution of maximum values')  +
  #geom_vline(xintercept = maxvalue) +
  theme_classic()
p3 = pp_check(model2, type = 'stat', stat = 'mean') +
  ggtitle('Prior predictive distribution of average values')   +
  #geom_vline(xintercept = meanvalue) +
  theme_classic()
p4 = pp_check(model2, ndraws = 100)
ggpubr::ggarrange(p1,p2,p3, p4)
conditional_effects(model2)

model2b <- brm(formula = smooth_straightness_to ~ s(trial_n) + season * speed_trip_to + status + season_status + (1 | ID),
               data = path_data, family = Beta(),
               prior = c(prior(normal(0, 1), class = "b"), 
                         prior(exponential(1), class = "sd")),
               chains = 4, iter = 4000, cores = 4, seed = 1234)
p1 = pp_check(model2b, type = 'stat', stat = 'min') +
  ggtitle('Prior predictive distribution of minimum values') 
p2 = pp_check(model2b, type = 'stat', stat = 'max') +
  ggtitle('Prior predictive distribution of maximum values')  +
  #geom_vline(xintercept = maxvalue) +
  theme_classic()
p3 = pp_check(model2b, type = 'stat', stat = 'mean') +
  ggtitle('Prior predictive distribution of average values')   +
  #geom_vline(xintercept = meanvalue) +
  theme_classic()
p4 = pp_check(model2b, ndraws = 100)
ggpubr::ggarrange(p1,p2,p3, p4)
model3 <- brm(bf(smooth_straightness_to ~ s(trial_n, by = status) + speed_trip_to + (1|ID)),
                    data = path_data, family = Beta(),
                    prior = c(prior(normal(0,1), class = "b"),
                              prior(exponential(0.2), class = "sds"),
                              prior(exponential(0.7), class = "sd", group = "ID")),
                    chains = 4, iter = 4000, cores = 4, seed = 1234,
                    backend = "cmdstanr", threads = threading(2),
                    control = list(adapt_delta = 0.99, max_treedepth = 15))
p1 = pp_check(model3, type = 'stat', stat = 'min') +
  ggtitle('Prior predictive distribution of minimum values') 
p2 = pp_check(model3, type = 'stat', stat = 'max') +
  ggtitle('Prior predictive distribution of maximum values')  +
  #geom_vline(xintercept = maxvalue) +
  theme_classic()
p3 = pp_check(model3, type = 'stat', stat = 'mean') +
  ggtitle('Prior predictive distribution of average values')   +
  #geom_vline(xintercept = meanvalue) +
  theme_classic()
p4 = pp_check(model3, ndraws = 100)
ggpubr::ggarrange(p1,p2,p3, p4)
conditional_effects(model3)

#the interaction terms (both status * speed_trip_to and season * speed_trip_to) do not explain much of the variability in smooth_straightness_to. 
#The main effects of predictors like season, status, and speed_trip_to are generally weak and not significantly different from zero.
# the smoothing spline captures some non-linear effects of trial number, but the uncertainty is too large to draw strong conclusions. 


get_prior(bf(smooth_straightness_to ~ s(trial_n) + status + speed_trip_to + season + (1 | ID)),
          data = path_data,family = mixture(Beta(), Beta()))

# Define the mixture model with two Beta distributions
fit_mixture <- brm(bf(smooth_straightness_to ~ s(trial_n) + status + speed_trip_to + season + season_status + (1 | ID)),
                   data = path_data,family = mixture(Beta(), Beta()),  # Mixture of two Beta distributions
                   prior = c(prior(normal(0, 1), class = "b", dpar = "mu1"),  
                             prior(exponential(1), class = "sd", dpar = "mu1"),
                             
                             # Priors for the second component (mu2)
                             prior(normal(0, 1), class = "b", dpar = "mu2"),  
                             prior(exponential(1), class = "sd", dpar = "mu2"),
                             
                             # Priors for the mixture weights (theta)
                             prior(dirichlet(1), class = "theta"),
                             
                             # Priors for the precision parameters (phi1 and phi2)
                             prior(gamma(0.01, 0.01), class = "phi1"),
                             prior(gamma(0.01, 0.01), class = "phi2")), # Dirichlet prior for the mixture weights (proportion from each Beta distribution)
                   chains = 4, iter = 4000, cores = 4, seed = 1234,
                   backend = "cmdstanr", threads = threading(2),
                   control = list(adapt_delta = 0.99, max_treedepth = 15))

p1 = pp_check(fit_mixture, type = 'stat', stat = 'min') +
  ggtitle('Prior predictive distribution of minimum values') 
p2 = pp_check(fit_mixture, type = 'stat', stat = 'max') +
  ggtitle('Prior predictive distribution of maximum values')  +
  #geom_vline(xintercept = maxvalue) +
  theme_classic()
p3 = pp_check(fit_mixture, type = 'stat', stat = 'mean') +
  ggtitle('Prior predictive distribution of average values')   +
  #geom_vline(xintercept = meanvalue) +
  theme_classic()
p4 = pp_check(fit_mixture, ndraws = 100)
ggpubr::ggarrange(p1,p2,p3, p4)
conditional_effects(fit_mixture)

fit_mixture2 <- brm(bf(smooth_straightness_to ~ s(trial_n, by = status) + speed_trip_to + season + season_status + (1 | ID)),
                   data = path_data,family = mixture(Beta(), Beta()),  # Mixture of two Beta distributions
                   prior = c(prior(normal(0, 1), class = "b", dpar = "mu1"),  
                             prior(normal(0, 1), class = "b", dpar = "mu2"),  
                             prior(exponential(1), class = "sd", dpar = "mu1"),
                             prior(exponential(1), class = "sd", dpar = "mu2"),
                             prior(dirichlet(1), class = "theta"),
                             prior(gamma(0.01, 0.01), class = "phi1"),
                             prior(gamma(0.01, 0.01), class = "phi2")),
                   chains = 4, iter = 4000, cores = 4, seed = 1234,
                   backend = "cmdstanr", threads = threading(2),
                   control = list(adapt_delta = 0.99, max_treedepth = 15))
p1 = pp_check(fit_mixture2, type = 'stat', stat = 'min') +
  ggtitle('Prior predictive distribution of minimum values') 
p2 = pp_check(fit_mixture2, type = 'stat', stat = 'max') +
  ggtitle('Prior predictive distribution of maximum values')  +
  #geom_vline(xintercept = maxvalue) +
  theme_classic()
p3 = pp_check(fit_mixture2, type = 'stat', stat = 'mean') +
  ggtitle('Prior predictive distribution of average values')   +
  #geom_vline(xintercept = meanvalue) +
  theme_classic()
p4 = pp_check(fit_mixture2, ndraws = 100)
ggpubr::ggarrange(p1,p2,p3, p4)

fit_mixture3 <- brm(bf(smooth_straightness_to ~ s(trial_n, by = season_status, k = 10) + speed_trip_to + (1 | ID)),
                    data = path_data,family = mixture(Beta(), Beta()),  # Mixture of two Beta distributions
                    prior = c(prior(normal(0, 0.5), class = "b", dpar = "mu1"),
                              prior(normal(0, 0.5), class = "b", dpar = "mu2"),
                              prior(exponential(0.5), class = "sds", dpar = "mu1"),  # Allow more flexibility for mu1
                              prior(exponential(0.5), class = "sds", dpar = "mu2"),  # Allow more flexibility for mu2
                              prior(dirichlet(1, 1), class = "theta"),  # Keep mixture components balanced
                              prior(gamma(2, 1), class = "phi1"),
                              prior(gamma(2, 1), class = "phi2")),
                    chains = 4, iter = 4000, cores = 4, seed = 1234,
                    backend = "cmdstanr", threads = threading(2),
                    control = list(adapt_delta = 0.99, max_treedepth = 15))
p1 = pp_check(fit_mixture3, type = 'stat', stat = 'min') +
  ggtitle('Prior predictive distribution of minimum values') 
p2 = pp_check(fit_mixture3, type = 'stat', stat = 'max') +
  ggtitle('Prior predictive distribution of maximum values')  +
  #geom_vline(xintercept = maxvalue) +
  theme_classic()
p3 = pp_check(fit_mixture3, type = 'stat', stat = 'mean') +
  ggtitle('Prior predictive distribution of average values')   +
  #geom_vline(xintercept = meanvalue) +
  theme_classic()
p4 = pp_check(fit_mixture3, ndraws = 100)
ggpubr::ggarrange(p1,p2,p3, p4)
conditional_effects(fit_mixture3)
conditional_smooths(fit_mixture3)

fit_mixture4 <- brm(bf(smooth_straightness_to ~ s(trial_n, by = season) + speed_trip_to + status + season_status + (1 | ID)),
                    data = path_data,family = mixture(Beta(), Beta()),  # Mixture of two Beta distributions
                    prior = c(prior(normal(0, 1), class = "b", dpar = "mu1"),  
                              prior(normal(0, 1), class = "b", dpar = "mu2"),  
                              prior(exponential(1), class = "sd", dpar = "mu1"),
                              prior(exponential(1), class = "sd", dpar = "mu2"),
                              prior(dirichlet(1), class = "theta"),
                              prior(gamma(0.01, 0.01), class = "phi1"),
                              prior(gamma(0.01, 0.01), class = "phi2")),
                    chains = 4, iter = 4000, cores = 4, seed = 1234,
                    backend = "cmdstanr", threads = threading(2),
                    control = list(adapt_delta = 0.99, max_treedepth = 15))
p1 = pp_check(fit_mixture4, type = 'stat', stat = 'min') +
  ggtitle('Prior predictive distribution of minimum values') 
p2 = pp_check(fit_mixture4, type = 'stat', stat = 'max') +
  ggtitle('Prior predictive distribution of maximum values')  +
  #geom_vline(xintercept = maxvalue) +
  theme_classic()
p3 = pp_check(fit_mixture4, type = 'stat', stat = 'mean') +
  ggtitle('Prior predictive distribution of average values')   +
  #geom_vline(xintercept = meanvalue) +
  theme_classic()
p4 = pp_check(fit_mixture4, ndraws = 100)
ggpubr::ggarrange(p1,p2,p3, p4)


path_data <- path_data %>%
  mutate(season_status = factor(season_status, levels = c("summer_wild", "winter_wild", "spring_wild",
                                                          "winter_captive", "spring_captive")))

new_d <- path_data %>%
  distinct(ID, status, trial_n, season, season_status, speed_trip_to)

fitted_all <- fitted(fit_mixture3, newdata = new_d) %>%
  data.frame() %>%
  bind_cols(new_d) %>%
  mutate(season_status = factor(season_status, levels = c("summer_wild", "winter_wild", "spring_wild",
                                                          "winter_captive", "spring_captive")),
         trial = factor(trial_n))

fitted_avg_all <- fitted_all %>%
  group_by(season_status, trial) %>%
  summarize(AvgEstimate = mean(Estimate), Lower = mean(Q2.5), Upper = mean(Q97.5), .groups = 'drop')

ggplot(data = path_data, mapping = aes(x = trial, y = smooth_straightness_to, color = season_status, group = season_status)) +
  geom_jitter(alpha = 0.3) +
  labs(color = "Status", fill = "CI 95%", x = "Trial", y = "Visual Cue - Path Efficiency") +
  scale_color_manual(values = c("summer_wild" = "#D16103", 
                                "winter_wild" = "#0A33A9", 
                                "spring_wild" = "#006400", 
                                "winter_captive" = "#0AA9A9", 
                                "spring_captive" = "#0fa90f"),
                     labels = c("Summer Wild", "Winter Wild", "Spring Wild",
                                "Winter Captive", "Spring Captive")) +
  scale_fill_manual(values = c("summer_wild" = "#D16103", 
                               "winter_wild" = "#0A33A9", 
                               "spring_wild" = "#006400", 
                               "winter_captive" = "#0AA9A9", 
                               "spring_captive" = "#0fa90f"),
                    labels = c("Summer Wild", "Winter Wild", "Spring Wild",
                               "Winter Captive", "Spring Captive")) +
  geom_line(data = fitted_avg_all, aes(x = trial, y = AvgEstimate, group = season_status), linewidth = 1) +
  geom_ribbon(data = fitted_avg_all, aes(x = trial, ymin = Lower, ymax = Upper, fill = season_status, group = season_status), alpha = 0.2, inherit.aes = FALSE) +
  facet_wrap(~season_status, labeller = labeller(season_status = as_labeller(c("summer_wild" = "Summer Wild", 
                                                                               "winter_wild" = "Winter Wild",  
                                                                               "spring_wild" = "Spring Wild",
                                                                               "winter_captive" = "Winter Captive", 
                                                                               "spring_captive" = "Spring Captive")))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 15, vjust = 0),
        axis.title.y = element_text(size = 15, vjust = 2),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(size = 12))

ggplot(data = path_data, mapping = aes(x = trial, y = smooth_straightness_to, color = season_status, group = season_status)) +
  geom_jitter(alpha = 0.3) +
  labs(color = "Status", fill = "CI 95%", x = "Trial", y = "Visual Cue - Path Efficiency") +
  scale_color_manual(values = c("#D16103", "#0AA9A9", "#0A33A9", "#009E73", "#006400"),
                     labels = c("Summer Wild", "Winter Captive", "Winter Wild", "Spring Captive", "Spring Wild")) +
  scale_fill_manual(values = c("#D16103", "#0AA9A9", "#0A33A9", "#009E73", "#006400"),
                    labels = c("Summer Wild", "Winter Captive", "Winter Wild", "Spring Captive", "Spring Wild")) +
  geom_smooth() +
  facet_wrap(~season_status, labeller = labeller(season_status = as_labeller(c("summer_wild" = "Summer Wild",
                                                                               "winter_captive" = "Winter Captive",
                                                                               "winter_wild" = "Winter Wild",
                                                                               "spring_captive" = "Spring Captive",
                                                                               "spring_wild" = "Spring Wild")))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 15, vjust = 0),
        axis.title.y = element_text(size = 15, vjust = 2),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(size = 12))

###Quantify difference in estimates between trial 1, 6, 8 and 10) ####

key_trials <- c(1, 5, 7, 8, max(path_data$trial_n))

key_fitted_avg <- fitted_avg_all %>% 
  filter(trial %in% key_trials) %>% 
  mutate(Odds = exp(AvgEstimate),
         LowerOdds = exp(Lower),
         UpperOdds = exp(Upper))



## Straightness back ####
path_data <- path_data %>%
  group_by(status) %>%
  mutate(speed_trip_back = ifelse(is.na(speed_trip_back),
                                mean(speed_trip_back, na.rm = TRUE),
                                speed_trip_back)) %>%
  ungroup()

fit_back_beta <- brm(bf(smooth_straightness_back ~ s(trial_n, by = status) + status*speed_trip_back + (1|ID)),
                    data = path_data, family = Beta(),
                    prior = c(prior(normal(-1,0.5), class = "b", coef = "statuswinter_captive"),
                              prior(normal(-1,0.5), class = "b", coef = "statuswinter_wild"),
                              prior(normal(0,1), class = "b", coef = "speed_trip_back"),
                              prior(normal(0,1), class = "b", coef = "statuswinter_captive:speed_trip_back"),
                              prior(normal(0,1), class = "b", coef = "statuswinter_wild:speed_trip_back"),
                              prior(normal(0,1), class = "b", coef = "strial_n:statuswinter_wild_1"),
                              prior(exponential(0.2), class = "sds"),
                              prior(exponential(0.65), class = "sd", group = "ID")),
                    chains = 4, iter = 4000, cores = 4, seed = 1234,
                    backend = "cmdstanr", threads = threading(2),
                    control = list(adapt_delta = 0.95, max_treedepth = 15))
new_back <- path_data %>%
  distinct(ID, status, trial, trial_n, season, speed_trip_back)
fitted_back <- fitted(fit_back_beta,
                     newdata = new_back) %>%
  data.frame() %>%
  bind_cols(new_back) %>%
  mutate(ID = str_c("ID[", ID, "]"),
         status = factor(status),
         trial = factor(trial))

fitted_avg_back <- fitted_back %>%
  group_by(status, trial) %>%
  summarize(AvgEstimate = mean(Estimate),
            Lower = mean(Q2.5),
            Upper = mean(Q97.5))

str(fitted_avg_back)
ggplot(data = path_data, mapping = aes(x = trial, y = beta_smooth_straightness_back, color = status, group = status)) +
  geom_jitter(alpha =0.3) +
  labs(color="Status", fill="CI 95%") +
  labs(x = "Trial", y = "Path Integration - Path efficiency") +
  scale_color_manual(values = c("#D16103","#0AA9A9", "#0A33A9"), labels = c("Summer Wild","Winter Captive", "Winter Wild")) +
  scale_fill_manual(values = c("#D16103","#0AA9A9", "#0A33A9"), labels = c("Summer Wild","Winter Captive", "Winter Wild")) +
  #scale_x_discrete(breaks = path_data$trial, labels = path_data$trial) +
  geom_line(data = fitted_avg_back, aes(x = trial, y = AvgEstimate, group = status), linewidth = 1) +
  geom_ribbon(data = fitted_avg_back,
              aes(x = trial, ymin = Lower, ymax = Upper, fill = status, group = status),
              alpha = 0.2, inherit.aes = FALSE) +
  facet_wrap(~status, labeller = labeller(status = as_labeller(c("summer_wild" = "Summer Wild",
                                                                 "winter_captive" = "Winter Captive",
                                                                 "winter_wild" = "Winter Wild")))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_text(size=15, vjust=0),
        axis.title.y = element_text(size=15, vjust=2),
        axis.text.x = element_text(size=12, vjust=0),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(size=12, vjust=0))
