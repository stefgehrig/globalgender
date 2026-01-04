library(tidyverse)
library(patchwork)
library(lavaan)
library(psych)
library(extrafont)
loadfonts(quiet = TRUE)

# source custom functions
source("R/functions.R")

# data config
min_t = ceiling((2022-1995)/4) # minimum time span between two surveys allowed (1/4 of period)
min_surv = 2 # minimum number of surveyed time points for a country
never_use = c("Northern Ireland") # which entities never used as countries in any analysis

# plot config
fontfam <- "Segoe UI"

# load data
df <- readRDS("data/wevsdata.rds") %>% filter(!country %in% never_use) 

#########################
#### data processing ####
#########################
dfr <- filter_process_data(df = df, outcome = "grb", min_t = min_t, min_surv = min_surv)

#########################
#### factor loadings ####
#########################
# fit model
f <- 'GRB =~ gender_jobs + gender_lead + gender_univ'
m <- cfa(model = f, data = dfr, group = "culzone")
fitMeasures(m, c("cfi", "tli", "rmsea", "srmr")) # perfect due to zero degrees of freedom (saturated)

# get parameter estimates and prepare labels
std_params <- standardizedSolution(m)
loadings <- std_params %>%
  filter(op == "=~", lhs == "GRB") %>%
  mutate(rhs = case_when(
    rhs == "gender_jobs" ~ "Item: Jobs",
    rhs == "gender_lead" ~ "Item: Politics",
    rhs == "gender_univ" ~ "Item: Education"
  )) %>% 
  select(group, rhs, est = est.std)

group_labels <- lavInspect(m, "group.label")
loadings$group_name <- group_labels[loadings$group]

# plotting
p1 <- ggplot(loadings, aes(x = group_name, y = est)) +
  geom_point() +
  facet_wrap(~ rhs, ncol = 3) +
  labs(x = "Cultural Zone", y = "Standardized factor loading") +
  theme_bw(12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_continuous(limits = c(-0,1)) + 
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        text = element_text(family = fontfam)) 

##########################
#### cronbach's alpha ####
##########################
# global
alpha(dfr[,c("gender_jobs", "gender_lead", "gender_univ")])$total["raw_alpha"] %>% round(2) # 0.66

# by zone
alphas <- map_dfr(unique(dfr$culzone), function(x){
  tibble(
    culzone = x,
    alpha = as.numeric(alpha(
      dfr[dfr$culzone==x, c("gender_jobs", "gender_lead", "gender_univ")]
    )$total["raw_alpha"])
  )
})

# plotting
p2 <- ggplot(alphas, aes(x = culzone, y = alpha)) +
  geom_point() +
  labs(x = "Cultural Zone", y = "Cronbach's \u03B1") +
  theme_bw(12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_continuous(limits = c(-0,1)) + 
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        text = element_text(family = fontfam)) 

# compose combined plot
p <- p1 + p2 + plot_layout(width = c(3,1)) + plot_annotation(tag_levels = "A")

png("results/FigSM_measure.png", width = 3600, height = 1800, res = 375)
p
dev.off()
