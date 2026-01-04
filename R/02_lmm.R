library(lme4)
library(tidyverse)
library(broom.mixed)
library(kableExtra)
library(ggtext)
library(patchwork)
library(marginaleffects)
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
    age >= 40 & age <= 49
  )

# set factor reference levels
df_cross <- df_cross %>% 
  mutate(across(.cols = c(education, income_lvl, townsize), factor),
         education = relevel(education, ref = "Primary or less"),
         income_lvl = relevel(income_lvl, ref = "Low"),
         townsize = relevel(townsize, ref = "<5,000")
  )

#######################
#### model fitting ####
#######################
m_cross_null <- model_fert_hierarchical(
  df_cross,
  outcome = "children_nr",
  covariates = NULL,
  with_age = FALSE
)

m_cross_demogr <- model_fert_hierarchical(
  df_cross,
  outcome = "children_nr",
  covariates = c("education", "townsize", "income_lvl"),
  with_age = TRUE
)

m_cross_demogr_relig <- model_fert_hierarchical(
  df_cross,
  outcome = "children_nr",
  covariates = c("education", "townsize", "income_lvl", "religiosity"),
  with_age = TRUE
)

restab_for_export <- build_cross_regrtable(
  list(m_cross_null$m, m_cross_demogr$m, m_cross_demogr_relig$m)
)

cat(restab_for_export, file = "results/Tab1_nrchildren.html")

##############
#### plot ####
##############
# plot config
fontfam <- "Segoe UI"
geqcol <- rev(paletteer::paletteer_d("ggsci::planetexpress_futurama", n = 3)[c(3,1)] )

# predictions
prds <- predictions(m,
            newdata = datagrid(grb_indiv = seq(-1, 1, 0.1), 
                               grb_mean = c(-1, 1)),
            conf_level = 0.9,
            re.form = NA) %>% as_tibble()

c <- attributes(m_cross_null$Xs$grb)$`scaled:center`
s <- attributes(m_cross_null$Xs$grb)$`scaled:scale`
m <- m_cross_null$m

prds_plot <- prds %>% 
  mutate(across(contains("grb_"), ~(.x*s)+c)) %>% # back-transform
  mutate(grb_mean = factor(round(grb_mean,2)))

# plotting
p_prds <- prds_plot %>% 
  ggplot() + 
  geom_line(aes(x = grb_indiv, y=estimate, group= grb_mean, col= grb_mean), lwd  = 1) + 
  geom_ribbon(aes(x = grb_indiv, ymin=conf.low, ymax = conf.high, group= grb_mean, fill= grb_mean),
              alpha= 0.1, show.legend = FALSE) + 
  theme_classic(14) +
  theme(text = element_text(family = fontfam),
        legend.position = "right",
        plot.caption = element_markdown()
        ) +
  scale_y_continuous(limits = c(1,3.5)) +
  scale_color_manual(values = geqcol) +
  scale_fill_manual(values = geqcol) +
  labs(y = "Average number of children",
       x = "GRB",
       col = "Country\nmean GRB",
       caption = paste0("*Note:* The GRB values of ", 
                        paste(sort(unique(prds_plot$grb_mean)), collapse = " and "), 
                        " correspond to mean \u00B1 1 SD<br>of the sample used for model estimation."))

png("results/FigSM_nrchildren.png", width = 2400, height = 1800, res = 360)
p_prds
dev.off()
