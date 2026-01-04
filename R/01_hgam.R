library(tidyverse)
library(mgcv)
library(gratia)
library(gammit)
library(brms)
library(tidybayes)
library(directlabels)
library(ggdist)
library(ggtext)
library(maps)
library(patchwork)
library(ggrepel)
library(paletteer)
library(countrycode)
library(geomtextpath)
library(extrafont)
loadfonts(quiet = TRUE)

# source custom functions
source("R/functions.R")

# data config
min_t = ceiling((2022-1995)/4) # minimum time span between two surveys allowed (1/4 of period)
min_surv = 2 # minimum number of surveyed time points for a country
never_use = c("Northern Ireland") # which entities never used as countries in any analysis

# analysis config
ndraws <- 2e3

# load data
df <- readRDS("data/wevsdata.rds") %>% filter(!country %in% never_use) 
df_census <- readRDS("data/UN_adultpop_bydemogr.rds")
df_census_world <- readRDS("data/UN_adultpop_world.rds")

#########################
#### data processing ####
#########################
dfr <- filter_process_data(df = df, outcome = "grb", min_t = min_t, min_surv = min_surv)
# n                      373855
# countries                  78
# culzones                    9
# years                      26
# surveyed_country_years    270

# % of adult global pop in 2022
adultsinsample <- df_census %>% filter(year == 2022) %>%
  filter(country %in% unique(dfr$country)) %>%
  summarise(population = sum(population)) %>% pull(population)
adultsinworld <- df_census_world %>% filter(year == 2022) %>% pull(population)
(adultsinsample / adultsinworld * 100) %>% round(1) #  86.1%

##################
#### modeling ####
##################
f <- as.formula(
  y_s ~
    s(year_s, country, k = 10, bs = "fs") +
    s(year_s, culzone, k = 10, bs = "fs") +
    s(demogr_grp_country, bs = "re") +
    demogr_grp
)

m <- bam(formula = f, data = dfr, family = "gaussian", nthreads = 15, method = "fREML")
saveRDS(m , file = "results/intermediate/m_hgam.rds")

##########################################
#### sampling and post-stratification ####
##########################################
df_pr <- process_hgam(m = m, df_census = df_census, ndraws = ndraws)

# obtain posterior draws for post-stratified estimates
df_est <- df_pr %>%
  group_by(year, country, culzone, .draw) %>%
  summarise(est = weighted.mean(e_grb, pr_stratum),
            population = sum(population),
            .groups = "drop")

##################
#### plotting ####
##################
# general plot config
fontfam <- "Segoe UI"
hdicols <- paletteer_d("colRoz::v_acanthurus", n = length(unique(dfr$hdicode)), direction = 1)
culzcols <- paletteer_d("ggsci::planetexpress_futurama", n = length(unique(dfr$culzone)))
countrycols <- colorRampPalette(paletteer_d("ggsci::planetexpress_futurama"))(length(unique(dfr$country)))
fillscale <- scale_fill_paletteer_c("scico::vikO", na.value="white", direction =-1, breaks = seq(-40,40,20), limits = c(-75,75))
cfcol <- paletteer_d("ggsci::planetexpress_futurama", n = 5)[c(1,4,5)] 

# histogram
p_hist <- dfr %>% 
  ggplot() + 
  geom_histogram(aes(x = grb), col = "white", fill = "grey50",
                 breaks = seq(0,100,10)) +
  theme_classic(14) +
  theme(text = element_text(family = fontfam)) + 
  labs(x = "GRB",
       y = "Frequency") + 
  scale_y_continuous(labels = ~format(.x, big.mark = ","))

png("results/FigSM_hist.png", width = 2750, height = 2000, res = 475)
p_hist
dev.off()

# compute time trends by group for plot
time_means_by_group <- get_time_means_by_group(df_est = df_est, df = df, ndraws = ndraws)

# ... print estimates for first and last year
time_means_by_group$df_glob %>% 
  filter(year %in% c(min(year), max(year)))
# 1995  52.7 3202302722.  51.8  53.5
# 2022  55.2 4823259018.  53.2  57.1
time_means_by_group$df_avgc %>% 
  filter(year %in% c(min(year), max(year)))
# 1995  52.1  50.9  53.3
# 2022  62.5  61.3  63.6

# plot time trends by group
p_time <- plot_time_means(df_glob = time_means_by_group$df_glob, 
                          df_avgc = time_means_by_group$df_avgc, 
                          df_culz = time_means_by_group$df_culz, 
                          df_hdi  = time_means_by_group$df_hdi)

# plot global map with posterior difference 1995-2022
df_snap <- df_est %>%
  filter(year %in% c(1995,2022)) %>% 
  group_by(country, .draw) %>%
  summarise(
    delta = est[year == 2022] - est[year == 1995],
    .groups = "drop_last"
  ) %>% 
  summarise(m = mean(delta),
            lwr = quantile(delta, 0.05),
            upr = quantile(delta, 0.95),
            .groups = "drop")

world <- as_tibble(map_data("world")) %>% 
  mutate(region = case_when(
    # replace some country names by hand which are stored differently
    subregion == "Macao" ~ "Macao",
    subregion == "Hong Kong" ~ "Hong Kong",
    region    == "Micronesia" ~ "Federated States of Micronesia",
    TRUE ~ region
  )) %>% 
  mutate(country = countrycode(region, "country.name", "country.name"))

df_snap_world <- df_snap %>%
  full_join(world, by = join_by(country)) %>%
  filter(region != "Antarctica")

p_mapdelta <- ggplot(df_snap_world, aes(x = long, y = lat, group = region, subgroup = group)) +
  geom_polygon(aes(fill = m), col = "black", linewidth = 0.4) +
  coord_fixed(1.15) +
  theme_void(14) +
  theme(
    text = element_text(family = fontfam),
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.3),
    legend.key.size = unit(0.3, "in"),
    legend.text = element_text(size = 11, family = fontfam),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  guides(fill = guide_legend(
    override.aes = list(col = "white"),
    ncol = 1,
    reverse = TRUE
  )) +
  labs(fill = "Change in mean GRB\n(1995-2022)") +
  fillscale

png("results/Fig1_change.png", width = 5000, height = 3600, res = 400)
p_time / p_mapdelta + plot_annotation(tag_levels = "A") + 
  plot_layout(heights = c(1,2))
dev.off()

# plot national trends grouped by cultural zone
df_cwise <- df_est %>%
  group_by(year, country) %>%
  summarise(m = mean(est),
            lwr = quantile(est, 0.05),
            upr = quantile(est, 0.95),
            culzone = unique(culzone),
            .groups = "drop")

p_allc <- df_cwise %>%
  ggplot() +
  geom_line(aes(x = year, y = m, col = country, group = country)) +
  geom_ribbon(aes(x = year, ymin = lwr, ymax = upr, fill = country), alpha = 0.15) +
  facet_wrap(~culzone,
             scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = seq(1995,2020,5), limits = c(1995, 2030)) +
  theme_classic(14) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        text = element_text(family = fontfam)) +
  geom_dl(aes(x = year, y = m, label = country, color = country),
          method = list(dl.trans(x = x + 0.025), cex = 0.55, list("last.points"))) +
  labs(x = "Year",
       y = "Mean GRB") + 
  scale_color_manual(values = countrycols) + 
  scale_fill_manual(values = countrycols) + 
  scale_y_continuous(breaks = seq(0,100,20))

png("results/FigSM_countrytrends.png", width = 4750, height = 5250, res = 460)
p_allc
dev.off()

# plot random intercepts for demographic groups
fixed <- tibble(
  var = names(fixef(m)),
  coef_fixed = fixef(m)) %>% 
  mutate(
    var = ifelse(var == "(Intercept)", "demogr_grpFemale: 18 - 34", var),
    coef_fixed = ifelse(var == "demogr_grpFemale: 18 - 34", 0, coef_fixed)
  ) %>% 
  mutate(var = str_remove(var, "demogr_grp"))

rdm <- extract_ranef(m, re = "demogr_grp_country") %>% 
  select(group, coef_rdm = value) %>% 
  separate(group, into = c("var", "country"), "_")

set.seed(1)
fixed_rdm <- left_join(rdm, fixed, by = join_by(var))%>% 
  mutate(var_num = as.numeric(factor(var)),
         var_num_jitt = var_num + runif(nrow(.), -0.175, 0.175),
         label_it = abs(coef_rdm) > 0.15) %>% 
  # get cultural zone labels
  left_join(
    dfr %>% distinct(country, culzone),
    by = join_by(country)
  ) %>% 
  # back-transform
  mutate(across(
    contains("coef_"), function(x){
      x * attr(m$model$y_s, "scaled:scale")
    }
  ))

p_demspread <- ggplot() + 
  geom_vline(aes(xintercept = 0), lty = 2) +
  theme_classic(14) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.caption = element_markdown(),
        text = element_text(family = fontfam)) + 
  labs(y = "Demographic group",
       x = "Average difference in GRB",
       col = "Cultural Zone",
       size = "",
       caption = "*Note:* Reference category is *Female: 18 - 34*.") + 
  geom_point(data = fixed_rdm,
             aes(x = coef_fixed + coef_rdm,
                 y = var_num_jitt,
                 col = culzone, size = "Country-specific effect"),
              alpha = 0.4) + 
  geom_point(data = fixed_rdm %>% distinct(var, .keep_all = TRUE),
             aes(x = coef_fixed, y= var_num, size = "Average effect")) + 
  geom_text_repel(data = fixed_rdm %>% filter(label_it),
                  aes(x = coef_fixed + coef_rdm,
                      y = var_num_jitt,
                      label = country),
                  size = 2.75,
                  box.padding = 0.1,
                  max.overlaps = 12,
                  alpha = 0.5,
                  family = fontfam) + 
  scale_y_continuous(
    labels = unique(fixed_rdm$var),
    breaks = 1:6
  ) + 
  scale_color_manual(values = culzcols) +
  scale_size_manual(values = c(7,2)) +
  guides(size = guide_legend(order = 1, ncol=1),
         col = guide_legend(ncol=2, title.position = "top", title.hjust = 0.5))

png("results/FigSM_demcoefs.png", width = 3750, height = 3750, res = 425)
p_demspread
dev.off()

# plot national trends against survey means
df_des_est <- dfr %>%
  group_by(year, country) %>%
  summarise(e_grb = mean(grb),
            culzone = unique(culzone),
            srv = paste0(paste(paste0(unique(srv), unique(wave)), collapse = ", "), " (", unique(year), ")"),
            .groups = "drop")

clzs <- sort(unique(df_est$culzone))
p_compare <- purrr::map(split(clzs, 1:length(clzs)),
                 function(clz){
                   df_est %>%
                     filter(culzone %in% clz) %>% 
                     ggplot() +
                     stat_lineribbon(aes(x = year, y = est), .width = c(0.5, 0.8, 0.9),
                                     point_interval = "mean_qi", color = "black", 
                                     fill = paletteer_d("ggsci::planetexpress_futurama", n = 5)[4], 
                                     lwd = 1, alpha = 0.25) +
                     geom_point(
                       data = df_des_est %>% filter(culzone %in% clz),
                       aes(x = year, y = e_grb),
                       shape = 1,
                       size = 5
                     ) +
                     geom_label(data = df_des_est %>% filter(culzone %in% clz),
                                aes(x = year, y = e_grb, label = srv), vjust = -1.5, size = 1.25,
                                family = fontfam) +
                     theme_bw(12) +
                     coord_cartesian(ylim = c(0,100)) +
                     labs(title = clz,
                          y = "Mean GRB",
                          x = "Year"
                     ) +
                     theme(text = element_text(family = fontfam),
                           strip.background = element_blank(),
                           strip.text = element_text(face = "bold")) +
                     facet_wrap(~country, ncol = 4) +
                     scale_x_continuous(breaks = seq(1995,2020,5))
                 })

for(i in 1:length(p_compare)){
  nctry <- length(unique(p_compare[[i]]$data$country))
  rows <- ceiling(nctry/4)
  png(paste0("results/FigSM_countries", i, ".png"), 
      width = 3400, 
      height = 430 + rows * 690,
      res = 355 + 12 * rows)
  print(p_compare[[i]])
  dev.off()
}

# plot the sensitivity analysis
# ... create alternative data set with other country-wise inclusion criteria
min_t_alternative = ceiling((2022-1995)/2) # minimum time span between two surveys allowed (1/2 of period)
min_surv_alternative <- 3 # minimum number of surveyed time points for a country 
dfr_s <- filter_process_data(df = df, outcome = "grb", min_t = min_t_alternative, min_surv = min_surv_alternative)
# n                      321863
# countries                  60
# culzones                    9
# years                      26
# surveyed_country_years    229

# ... % of adult global pop in 2022
adultsinsample <- df_census %>% filter(year == 2022) %>%
  filter(country %in% unique(dfr_s$country)) %>%
  summarise(population = sum(population)) %>% pull(population)
adultsinworld <- df_census_world %>% filter(year == 2022) %>% pull(population)
(adultsinsample / adultsinworld * 100) %>% round(1) # 79.7

# ... model fit
m_s <- bam(formula = f, data = dfr_s, family = "gaussian", nthreads = 15, method = "fREML")
saveRDS(m_s , file = "results/intermediate/m_hgam_s.rds")

# ... processing
df_pr_s <- process_hgam(m = m_s, df_census = df_census, ndraws = ndraws)
df_est_s <- df_pr_s %>%
  group_by(year, country, culzone, .draw) %>%
  summarise(est = weighted.mean(e_grb, pr_stratum),
            population = sum(population),
            .groups = "drop")
time_means_by_group_s <- get_time_means_by_group(df_est = df_est_s, df = df, ndraws = ndraws)

# ... plotting
p_time_s <- plot_time_means(df_glob = time_means_by_group_s$df_glob, 
                            df_avgc = time_means_by_group_s$df_avgc, 
                            df_culz = time_means_by_group_s$df_culz, 
                            df_hdi  = time_means_by_group_s$df_hdi)

png("results/FigSM_changesens.png", width = 5000, height = 1500, res = 360)
p_time_s
dev.off()

############################################
#### association with population growth ####
############################################
# posterior draws for global post-stratified estimates 1995 and 2022, current populations
df_est_curr <- df_est %>% filter(year %in% c(1995,2022)) %>% 
  group_by(year, .draw) %>%
  summarise(est_curr = weighted.mean(x = est, w = population),
            pop = sum(population),
            .groups = "drop")

# posterior draws for global post-stratified estimate 2022, constant populations from the past
df_est_cons <- df_pr %>%
  filter(year %in% c(1995,2022)) %>% 
  group_by(demogr_grp_country) %>%
  mutate(pr_stratum_2022 = unique(pr_stratum[2022]),
         population_1995 = unique(population[1995])) %>%
  group_by(year, country, culzone, .draw) %>%
  summarise(est = weighted.mean(e_grb, pr_stratum_2022),
            population = sum(population_1995),
            .groups = "drop") %>% 
  group_by(year, .draw) %>%
  summarise(est_cons = weighted.mean(x = est, w = population),
            pop_cons = sum(population),
            .groups = "drop") %>% 
  filter(year == 2022)

# plot posterior means with current vs constant population weights
df_cons <- left_join(df_est_curr, df_est_cons, by = join_by(year, .draw)) %>%
  pivot_longer(cols = contains("est_"), names_to = "type", values_to = "est") %>% 
  filter(!(year == 1995 & type == "est_cons")) %>% 
  mutate(type = case_when(
    type == "est_curr" ~ "Observed",
    type == "est_cons" ~ "Counterfactual:\nconstant national population shares"
  ))

p_cons <- df_cons %>% 
  ggplot() + 
  stat_pointinterval(
    aes(x = factor(year), y = est, col = type, fill = type), 
    alpha = 0.5,
    position = position_dodge(0.25),
    size = 8,
    .width = c(.9),
    point_interval = "mean_qi"
  ) + 
  coord_flip(ylim = c(49,63)) +
  theme_minimal(14) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        plot.caption = element_text(color = cfcol[1], size = 12),
        text = element_text(family = fontfam)) +
  labs(y = "Global mean GRB",
       x = "Year") + 
  scale_color_manual(values = cfcol) + 
  scale_fill_manual(values = cfcol) + 
  geom_text(
    data = df_cons %>%
      group_by(year, type) %>%
      summarise(m = mean(est), .groups = "drop"),
    aes(
      x = factor(year),
      y = m,
      label = ifelse(year == 1995, NA_character_, type),
      col = type
    ),
    nudge_x = c(NA, -0.35, 0.3),
    size = 3.5,
    show.legend = FALSE,
    family = fontfam,
    check_overlap = TRUE
  )

# ... print estimate
df_est_cons %>% 
  group_by(year) %>% 
  summarise(m = mean(est_cons),
            lwr = quantile(est_cons, 0.05),
            upr = quantile(est_cons, 0.95))
# 2022  57.0  55.0  59.1

# plot correlation of population size and GRB at two different time points
df_cor <- df_est %>% 
  filter(year %in% c(1995,2022)) %>% 
  group_by(year, .draw) %>% 
  summarise(r = cor(est, population, method = "spearman"),
            .groups = "drop")

p_cor <- df_cor  %>% 
  ggplot() + 
  stat_pointinterval(
    aes(x = factor(year), y= r), 
    alpha = 0.5, size = 8, .width = c(.9),
    point_interval = "mean_qi", key_glyph = draw_key_rect
  ) +
  coord_flip() +
  theme_minimal(14) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        text = element_text(family = fontfam)) +
  labs(y = "Spearman correlation of country's\nadult population size and mean GRB",
       x = "")

# ... print estimates
df_cor  %>% 
  group_by(year) %>% 
  summarise(m = mean(r),
            lwr = quantile(r, 0.05),
            upr = quantile(r, 0.95))
# 1995 -0.0261 -0.0933  0.0416
# 2022 -0.264  -0.315  -0.212 

# plot a bayesian regression estimate with measurement error
# ... prepare dataframe for modeling
df_mod <- df_est %>%
  mutate(log_population = log(population)
  ) %>% select(-population) %>% 
  pivot_wider(names_from = year, values_from = c(est, log_population)) %>% 
  group_by(country) %>% 
  summarise(
    mgrb = mean(est_1995),
    sgrb = sd(est_1995),
    log_population_1995 = as.numeric(unique(log_population_1995)),
    log_population_2022 = as.numeric(unique(log_population_2022)),
    .groups = "drop")

# ... define brms formula
formula1 <- bf(log_population_2022 ~ offset(log_population_1995) + mi(mgrb)) + gaussian(identity)
formula2 <- bf(mgrb | mi(sdy = sgrb) ~ 0) + gaussian(identity)

# ... run sampling
fit_pop <- brm(
  formula   = formula1 + formula2 + set_rescor(FALSE),
  data      = df_mod,
  control   = list(adapt_delta = 0.99),
  warmup    = 1e3,
  iter      = 2e3,
  chains    = 2,
  cores     = 2,
  seed      = 1,
  backend   = "cmdstan",
  save_pars = save_pars(latent = FALSE)
)

# ... create predictions
pred_fct <- fit_pop %>%
  spread_draws(b_logpopulation2022_Intercept,
               bsp_logpopulation2022_mimgrb) %>% 
  expand_grid(x_grb = seq(0,100,1)) %>% 
  mutate(y_fct = exp(b_logpopulation2022_Intercept + x_grb * bsp_logpopulation2022_mimgrb))

# ... plotting
p_fct <- pred_fct %>% 
  group_by(x_grb) %>%
  summarise(y = mean(y_fct),
            lwr = quantile(y_fct, 0.05),
            upr = quantile(y_fct, 0.95)) %>% 
  ggplot() + 
  geom_texthline(yintercept = 1, lty = 2, label = "No change", size = 4.5, family = fontfam) +
  geom_line(aes(x = x_grb, y = y), lwd = 1) +
  geom_ribbon(aes(x = x_grb, ymin = lwr, ymax = upr), alpha = 0.2) +
  scale_y_continuous(breaks = seq(0.5,2.5,0.5)) +
  scale_x_continuous(breaks = seq(0,100,20)) +
  labs(x = "Country's mean GRB in 1995",
       y = "\nExpected factor change in\ncountry's adult population size\n1995 to 2022") +
  theme_classic(14) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        text = element_text(family = fontfam))

# compose all plots
png("results/Fig2_pop.png", width = 4400, height = 3400, res = 510)
p_fct / (p_cons + p_cor) + plot_layout(heights = c(1.1,1)) + plot_annotation(tag_levels = "A")
dev.off()
