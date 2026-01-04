# load libraries
library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(patchwork)
library(ggtext)
library(posterior)
library(pracma)
library(tidybayes)
library(extrafont)
loadfonts(quiet = TRUE)

# source custom functions
source("R/functions.R")

# data config
min_year = 2010 # minimum year for cross sectional analysis
never_use = c("Northern Ireland") # which entities never used as countries in any analysis

# load data
df <- readRDS("data/wevsdata.rds") %>% filter(!country %in% never_use)

##########################
#### data preparation ####
##########################
df_cross <- df %>% 
  # must have grb measured and age and number of children
  filter(!is.na(grb),
         !is.na(children_nr),
         !is.na(age)) %>% 
  # must be a recent survey
  filter(year >= min_year) %>% 
  # and only one survey per country
  group_by(country) %>% 
  filter(year == max(year)) %>% 
  ungroup %>% 
  # only women in relevant age bracket
  filter(
    sex == "Female",
    age >= 18 & age <= 49
  )%>% 
  # split by grb
  mutate(med_grb = median(grb),
         grb = as.numeric(grb > med_grb))
cutval <- unique(df_cross$med_grb) # 66.66667

# data list for stan
dat <- df_cross
dat$country_num <- as.numeric(factor(dat$country, levels = unique(dat$country)))

datg <- expand_grid(country_num = sort(unique(dat$country_num)),
                    grb = sort(unique(dat$grb))) %>% 
  mutate(ctry_grb_idx = 1:nrow(.))

datd <- dat %>% 
  group_by(country_num, grb, age, country) %>% 
  summarise(y = sum(children_nr),
            m = n(),
            .groups = "drop") %>% 
  left_join(datg, by = join_by(country_num, grb))

asfr_data <- list(
  # scalars
  N = nrow(datd),
  K = length(unique(datd$country_num)),
  J = nrow(datg),
  # by country, grb, age
  nobs = datd$m,
  ys = datd$y,
  ages = datd$age,
  pred_idcs = datd$ctry_grb_idx,
  # predictor values
  pred_grb = datg$grb,
  pred_cty = datg$country_num
)

# descriptive table for cell sizes
tab_asfrsample <- dat %>% 
  mutate(country_year = paste0(country, " (", year, ")"),
         grb = ifelse(grb == 1, 
                      paste0("Higher (> ", round(cutval, 2), ")"),
                      paste0("Lower (\u2264 ", round(cutval, 2), ")"))) %>%
  group_by(country_year, grb) %>% 
  summarise(n = n(), .groups = "drop_last") %>% 
  mutate(p = n/sum(n)) %>% 
  ungroup %>% 
  mutate(
    n = paste0(n, " (",
      format(round(p*100), trim = TRUE), "%)"
    )
  ) %>% 
  select(-p) %>% 
  rename(
    "Country (year)" = country_year,
    "GRB subpopulation" = grb,
    "Sample size" = n
  )

saveRDS(tab_asfrsample, file = "results/TabSM_asfrsample.rds")

############################
#### posterior sampling ####
############################
# run model in cmdstan
mod <- cmdstan_model("Stan/asfrmodel.stan")

fit <- mod$sample(
  data = asfr_data,
  seed = 1,
  init = 0.01,
  chains = 2,
  parallel_chains = 2,
  iter_warmup = 1e3,
  iter_sampling = 1e3,
  adapt_delta = 0.99,
  max_treedepth = 10,
  refresh = 25,
  output_dir = "results/intermediate"
)

saveRDS(fit, "results/intermediate/asfrmodel.rds")

######################
#### diagnostics #####
######################
# define sets of parameters of interest
curvepars <- c("c1", "mu", "sigma1","sigma2")

# print summary of population average and variance parameters
summary(fit$draws(get_variables(fit)[grepl("beta", get_variables(fit))]))
summary(fit$draws(get_variables(fit)[grepl("std", get_variables(fit))]))
summary(fit$draws(get_variables(fit)[grepl("^Omega", get_variables(fit))])) %>% print(., n = 1e2)

# some mcmc checks
fit$diagnostic_summary()
mcmc_trace(fit$draws(), regex_pars = "^Omega")
mcmc_trace(fit$draws(), regex_pars = "std")
mcmc_trace(fit$draws(), pars = paste0("sigma1[", 1:12, "]"))
mcmc_trace(fit$draws(), pars = paste0("c1[", 51:62, "]"))
mcmc_trace(fit$draws(), pars = paste0("b[", 1:12, ",1]"))
mcmc_trace(fit$draws(), pars = paste0("b[1,", 1:8, "]"))
mcmc_trace(fit$draws(), pars = paste0("b[11,", 1:8, "]"))
mcmc_acf(fit$draws(), regex_pars = "std")
draws_df <- fit$draws(curvepars, format = "df")
rhats <- draws_df %>% 
  select(!contains(".")) %>% 
  summarise_all(~rhat(matrix(., ncol = max(draws_df$.chain)))) %>%
  pivot_longer(everything(), values_to = "rhat")
mcmc_rhat_hist(rhats$rhat)
lp_cp <- log_posterior(fit)
np_cp <- nuts_params(fit)
mcmc_nuts_divergence(np_cp, lp_cp)

##############
#### plot ####
##############
# plot config
fontfam <- "Segoe UI"
geqcol <- paletteer::paletteer_d("ggsci::planetexpress_futurama", n = 3)[c(3,1)] 

# constants for rescaling (as used for re-parametrization in stan code)
log_c1_center = log(0.125)
log_mu_center = log(25)
log_sigma_center = log(5)
log_c1_scale = 0.5
log_mu_scale = 0.2  
log_sigma_scale = 0.75

fit_derived <- fit %>%
  spread_draws(beta[i]) %>% 
  ungroup %>% 
  pivot_wider(names_from = i, values_from = beta, names_prefix = "beta") %>% 
  mutate(
    c1_l =     exp(beta1 * log_c1_scale + log_c1_center), # for average country (b = 0), by using beta
    mu_l =     exp(beta2 * log_mu_scale + log_mu_center),
    sigma1_l = exp(beta3 * log_sigma_scale + log_sigma_center),
    sigma2_l = exp(beta4 * log_sigma_scale + log_sigma_center),
    c1_h =     exp(beta5 * log_c1_scale + log_c1_center),
    mu_h =     exp(beta6 * log_mu_scale + log_mu_center),
    sigma1_h = exp(beta7 * log_sigma_scale + log_sigma_center),
    sigma2_h = exp(beta8 * log_sigma_scale + log_sigma_center)
  )

# compute derived quantities
fert_target = 0.5 # age at which avg fertility equals this
min_age = 0 # integration lower...
max_age = 49 # ...and upper limit 

fit_derived <- fit_derived %>%
  mutate(
    # integrate over first and second part of asfr curve
    tfr1_h =  c1_h * ((sqrt(pi) * erf((mu_h - min_age) / sigma1_h) * sigma1_h) / 2),
    tfr1_l =  c1_l * ((sqrt(pi) * erf((mu_l - min_age) / sigma1_l) * sigma1_l) / 2),
    tfr2_h = -c1_h * ((sqrt(pi) * erf((mu_h - max_age) / sigma2_h) * sigma2_h) / 2),
    tfr2_l = -c1_l * ((sqrt(pi) * erf((mu_l - max_age) / sigma2_l) * sigma2_l) / 2),
    # integrate over full asfr curve
    tfr_h = tfr1_h + tfr2_h,
    tfr_l = tfr1_l + tfr2_l,
    # find age at which integral equals fert_target = 0.5
    agef1_h = if_else(tfr1_h >= fert_target, 
                      mu_h - erfinv(erf(mu_h / sigma1_h - min_age / sigma1_h) - (2*fert_target) / (sqrt(pi) * c1_h * sigma1_h)) * sigma1_h,
                      mu_h + erfinv((2*(fert_target - tfr1_h)) / (sqrt(pi) * c1_h * sigma2_h)) * sigma2_h),
    agef1_l = if_else(tfr1_l >= fert_target,
                      mu_l - erfinv(erf(mu_l / sigma1_l - min_age / sigma1_l) - (2*fert_target) / (sqrt(pi) * c1_l * sigma1_l)) * sigma1_l,
                      mu_l + erfinv((2*(fert_target - tfr1_l)) / (sqrt(pi) * c1_l * sigma2_l)) * sigma2_l),
  ) %>%
  # remove helper columns
  select(-contains("tfr1"), -contains("tfr2"))

# plot the asfr curve posterior draws ('average country', i.e., all b at mean)
fit_byage <-  fit_derived %>% 
  expand_grid(age = seq(18,49,0.01)) %>% 
  mutate(asfr_l = c1_l * exp(-((age - mu_l) / ((age <= mu_l) * sigma1_l + (age > mu_l) * sigma2_l))^2),
         asfr_h = c1_h * exp(-((age - mu_h) / ((age <= mu_h) * sigma1_h + (age > mu_h) * sigma2_h))^2))

set.seed(1)
p_asfr <- fit_byage %>% 
  filter(.draw %in% sample(unique(.draw), size = 1e2)) %>% 
  pivot_longer(cols = c(asfr_h, asfr_l), values_to = "pred", names_to = "geq_cat") %>%
  mutate(geq_cat = ifelse(geq_cat == "asfr_h", 
                          paste0("Higher (> ", round(cutval, 2), ")"),
                          paste0("Lower (\u2264 ", round(cutval, 2), ")"))) %>%
  ggplot() +
  geom_line(aes(x = age, y = pred, group = paste0(geq_cat, .draw), col = geq_cat), 
            lwd = 0.25,
            alpha = 0.4) +
  theme_classic(14) +
  scale_color_manual(values = geqcol) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        text = element_text(family = fontfam)) + 
  labs(y = "Age-specific fertility rate\nin the average country",
       x = "Age",
       col = "GRB subpopulation") + 
  coord_cartesian(xlim = c(18,49)) + 
  scale_x_continuous(limits = c(18,49), breaks = c(18, seq(20,45,5), 49)) +
  guides(col = guide_legend(override.aes = list(alpha = 1, lwd = 1)))

# ... with inset tfr
insets_tfr <- fit_derived %>%
  pivot_longer(cols = contains("tfr"),
               values_to = "Total fertility rate", names_to = "geq_cat") %>% 
  ggplot() +
  stat_pointinterval(
    aes(x = `Total fertility rate`, col = geq_cat), .width = c(.9),
    show.legend = FALSE, position = position_dodge(0.5),
    point_interval = "mean_qi"
  ) +
  theme_bw(10) +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid= element_blank(),
        text = element_text(family = fontfam),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank()) +
  scale_color_manual(values = geqcol) +
  labs(y = "")+
  plot_layout(tag_level = 'new')

# ... with inset age
insets_age <- fit_derived %>%
  pivot_longer(cols = contains("agef1"),
               values_to = "Age at average fertility 0.5", names_to = "geq_cat")  %>%
  ggplot() +
  stat_pointinterval(
    aes(x = `Age at average fertility 0.5`, col = geq_cat), .width = c(.9),
    show.legend = FALSE, position = position_dodge(0.5),
    point_interval = "mean_qi"
  ) +
  theme_bw(10) +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid= element_blank(),
        text = element_text(family = fontfam),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank()) +
  scale_color_manual(values = geqcol) +
  labs(y = "")+
  plot_layout(tag_level = 'new')
  
# ... and print
fit_derived %>% 
  mutate(delta_age = fit_derived$agef1_h - fit_derived$agef1_l,
         delta_tfr = fit_derived$tfr_h - fit_derived$tfr_l) %>% 
  summarise(
    mean(delta_age),
    quantile(delta_age, 0.05),
    quantile(delta_age, 0.95),
    mean(delta_tfr),
    quantile(delta_tfr, 0.05),
    quantile(delta_tfr, 0.95)
  ) 
# posterior differences:
# delta_age: 1.52  [1.15; 1.90]
# delta_tfr: -0.173 [-0.256; -0.0899]

# plot posterior probability for rankings
fit_nations <- fit %>%
  spread_draws(b[country_num, i]) %>% 
  ungroup %>% 
  pivot_wider(names_from = i, values_from = c(b), names_prefix = "b_") %>% 
  mutate(
    c1_l =     exp(b_1 * log_c1_scale + log_c1_center),
    mu_l =     exp(b_2 * log_mu_scale + log_mu_center),
    sigma1_l = exp(b_3 * log_sigma_scale + log_sigma_center),
    sigma2_l = exp(b_4 * log_sigma_scale + log_sigma_center),
    c1_h =     exp(b_5 * log_c1_scale + log_c1_center),
    mu_h =     exp(b_6 * log_mu_scale + log_mu_center),
    sigma1_h = exp(b_7 * log_sigma_scale + log_sigma_center),
    sigma2_h = exp(b_8 * log_sigma_scale + log_sigma_center)
  )  %>% 
  # compute derived quantities
  mutate(
    tfr1_h =  c1_h * ((sqrt(pi) * erf((mu_h - min_age) / sigma1_h) * sigma1_h) / 2),
    tfr1_l =  c1_l * ((sqrt(pi) * erf((mu_l - min_age) / sigma1_l) * sigma1_l) / 2),
    tfr2_h = -c1_h * ((sqrt(pi) * erf((mu_h - max_age) / sigma2_h) * sigma2_h) / 2),
    tfr2_l = -c1_l * ((sqrt(pi) * erf((mu_l - max_age) / sigma2_l) * sigma2_l) / 2),
    tfr_h = tfr1_h + tfr2_h,
    tfr_l = tfr1_l + tfr2_l,
    agef1_h = if_else(tfr1_h >= fert_target, 
                      mu_h - erfinv(erf(mu_h / sigma1_h - min_age / sigma1_h) - (2*fert_target) / (sqrt(pi) * c1_h * sigma1_h)) * sigma1_h,
                      mu_h + erfinv((2*(fert_target - tfr1_h)) / (sqrt(pi) * c1_h * sigma2_h)) * sigma2_h),
    agef1_l = if_else(tfr1_l >= fert_target,
                      mu_l - erfinv(erf(mu_l / sigma1_l - min_age / sigma1_l) - (2*fert_target) / (sqrt(pi) * c1_l * sigma1_l)) * sigma1_l,
                      mu_l + erfinv((2*(fert_target - tfr1_l)) / (sqrt(pi) * c1_l * sigma2_l)) * sigma2_l),
  ) %>% 
  select(-contains("tfr1"), -contains("tfr2") )%>%
  # add country name
  left_join(
    datd %>% distinct(country_num, country), by = join_by(country_num)
  )

# plot percentages
p_percent <- fit_nations %>% 
  group_by(.draw) %>% 
  summarise(tfr_l_higher = mean(tfr_l > tfr_h),
            agef1_l_lower = mean(agef1_l < agef1_h, na.rm = TRUE)) %>%
  pivot_longer(cols = contains("_l_"), names_to = "y", values_to = "x") %>% 
  
  mutate(
    y = factor(ifelse(y=="tfr_l_higher",
                "... higher\ntotal fertility rate",
                "... lower age\nat average fertility 0.5")),
    y = fct_rev(y)
  ) %>% 
  ggplot() + 
    stat_pointinterval(
      aes(x = x, y= y), .width = c(.9),
      show.legend = FALSE,
      size = 8,
      height = 0.45,
      point_interval = "mean_qi"
    ) +
  scale_x_continuous(limits = c(0,1), labels = scales::percent_format()) + 
  scale_y_discrete(expand = c(0.12,0.12)) +
  coord_flip() + 
  labs(x = paste0("Proportion of countries where\nsubpopulation with lower GRB (\u2264 ", round(cutval, 2), ") has ..."),
       y = ""
  ) +
  theme_classic(14) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.caption = element_markdown(),
        text = element_text(family = fontfam)) 

png("results/Fig3_asfr.png", width = 5200, height = 2600, res = 500)
p_asfr +
  inset_element(insets_tfr, 0.6, 0.75, 0.99, 0.975) +
  inset_element(insets_age, 0.6, 0.5, 0.99, 0.725) +
  p_percent + plot_layout(widths = c(1, 0.5)) + plot_annotation(tag_levels = "A") & 
  theme(
    legend.position = "inside",
    legend.key.spacing = unit(0.04, "in"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.margin = margin(0.15,0.15,0.15,0.15, unit ="cm"),
    legend.background = element_rect(color = "black"),
    legend.position.inside = c(0.365, 0.175)
  ) 
dev.off()

# ... and print
fit_nations %>% 
  group_by(.draw) %>% 
  summarise(tfr_l_higher = mean(tfr_l > tfr_h),
            agef1_l_lower = mean(agef1_l < agef1_h, na.rm = TRUE)) %>%
  pivot_longer(cols = contains("_l_"), names_to = "y", values_to = "x") %>% 
  group_by(y) %>% 
  summarise(
    mean(x),
    quantile(x, 0.05),
    quantile(x, 0.95)
  )

# posterior probabilities:
# agef1_l_lower: 0.920 [0.861; 0.970]
# tfr_l_higher: 0.823 [0.713; 0.921]

# plot cumulative asfr by country
# ... prepare data frame
fit_check <- fit %>%
  spread_draws(lambda[ctry_grb_idx, age]) %>% 
  ungroup() %>% 
  mutate(age = age-1) %>% # index was shifted by 1
  left_join(
    datd %>% distinct(ctry_grb_idx, country, grb),
    by = join_by(ctry_grb_idx)
  ) %>% 
  left_join(
    dat %>% distinct(country, culzone),
    by = join_by(country)
  ) %>% 
  left_join(
    df_cross %>% 
      group_by(country, age, grb) %>% 
      summarise(children_nr = list(children_nr),
  .groups = "drop"),
    by = join_by(country, age, grb)
  )

# ... plotting
clzs <- sort(unique(fit_check$culzone))
p_compare <- purrr::map(split(clzs, 1:length(clzs)),
                        function(clz){
                        
                          fit_check_filtered <- fit_check %>% 
                            filter(culzone %in% clz) %>%
                            mutate(grb = ifelse(grb == 1, 
                                                paste0("Higher\n(> ", round(cutval, 2), ")"),
                                                paste0("Lower\n(\u2264 ", round(cutval, 2), ")"))) %>%
                            left_join(dat %>% distinct(country, year), by = join_by(country)) %>% 
                            mutate(countrylabel = paste0("**", country, "**<br>(", year, ")")) %>% 
                            filter(age >= 18) # just visualizing from 18, as for average ASFR curve
                            
                            fit_check_obs <- fit_check_filtered %>% 
                              distinct(countrylabel, grb, age, children_nr) %>% 
                              filter(map_dbl(children_nr, ~length(.x))!=0) %>% 
                              mutate(mchild = map_dbl(children_nr, ~mean(.x)),
                                     n =  map_dbl(children_nr, ~length(.x)))
                            
                            fit_check_filtered %>% 
                              group_by(age, grb, countrylabel) %>% 
                              summarise(postmeanlambda = mean(lambda), .groups = "drop") %>% 
                                ggplot() +
                                  geom_line(aes(x = age, y = postmeanlambda, col = grb), lwd = 1) +
                                  theme_bw(12) +
                                  scale_x_continuous(limits = c(18,49), breaks = c(18, seq(30,40,10), 49)) +
                                  scale_y_continuous(limits = c(0,NA)) +
                                  labs(title = clz,
                                      y = "Cumulative ASFR",
                                      size = "Sample size",
                                      x = "Age",
                                      color = "GRB\nsubpopulation") +
                                  theme(text = element_text(family = fontfam),
                                        strip.background = element_blank(),
                                        strip.text = element_markdown()) + 
                                  geom_point(
                                    data = fit_check_obs,
                                    aes(x = age, y = mchild, size = n, color = grb),
                                    shape = 1
                                  ) +
                                  facet_wrap(~countrylabel, ncol = 4) +
                                  scale_color_manual(values = geqcol) +
                              guides(color = guide_legend(override.aes = list(size = 0)))
                        })

for(i in 1:length(p_compare)){
  nctry <- length(unique(p_compare[[i]]$data$countrylabel))
  rows <- ceiling(nctry/4)
  png(paste0("results/FigSM_asfr", i, ".png"), 
      width = 3500, 
      height = 430 + rows * 690,
      res = 365 + 10 * rows)
  print(p_compare[[i]])
  dev.off()
}

###############################################
#### tables with model parameter estimates ####
###############################################
tab_asfrpostpars <- fit_derived %>% 
  summarise(across(.cols = (contains("_l") | contains("_h")) & !contains("tfr") & !contains("agef1"), 
                   ~format(round(mean(.x), 2), nsmall = 2))) %>% 
  rename(
    `$c^0$` = c1_l,
    `$\\mu^0$` = mu_l,
    `$\\sigma_1^0$` = sigma1_l,
    `$\\sigma_2^0$` = sigma2_l,
    `$c^1$` = c1_h,
    `$\\mu^1$` = mu_h,
    `$\\sigma_1^1$` = sigma1_h,
    `$\\sigma_2^1$` = sigma2_h
  )

saveRDS(tab_asfrpostpars, file = "results/TabSM_asfrpostpars.rds")

tab_asfrpostcors <- fit %>%
  spread_draws(Omega[i,j]) %>% 
  group_by(i, j) %>% 
  summarise(post_mean = mean(Omega), .groups = "drop") %>% 
  mutate(i = paste0("$b_{", i, "}$"),
         j = paste0("$b_{", j, "}$")) %>%
  pivot_wider(names_from = j, values_from = post_mean) %>%
  column_to_rownames("i") %>% 
  mutate(across(.cols = everything(), ~format(round(.x, 2), nsmall = 2)))

saveRDS(tab_asfrpostcors, file = "results/TabSM_asfrpostcors.rds")
