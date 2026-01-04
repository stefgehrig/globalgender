# back-transformation of scaled variable
unscale <- function(z, center = attr(z, "scaled:center"), scale = attr(z, "scaled:scale")) {
  z <- z * scale
  z <- z + center
  return(z)
}

# add columns wih age groups and demographic groups
add_age_groups <- function(data){
  data %>% 
    mutate(age_grp = case_when(
      age <= 34  ~ "18 - 34",
      age <= 59  ~ "35 - 59",
      age >= 60  ~ "60 +",
    )) %>% 
    mutate(
      demogr_grp = paste0(sex, ": ", age_grp)
    )
}

# process and filter the individual level data for fitting the hgam
filter_process_data <- function(df, outcome, min_t, min_surv){
  
  # create outcome
  df$y <- df %>% pull(!!as.symbol(outcome))
  # filter
  df <- df %>% filter(!is.na(sex), !is.na(y), !is.na(age), age >= 18) %>% 
    filter(!country %in% never_use)
  nrow(df) # 416,021
  # create age groups and demographic groups
  df <- add_age_groups(df)
  # arrange
  df <- df %>% 
    arrange(country, age_grp, sex, year)
  
  # filter for reduced data: only use if minimum time span is satisfied
  dfr <- df %>% 
    group_by(country) %>% 
    mutate(timespan = max(year)-min(year),
           numbersurvs = length(unique(year))) %>% 
    ungroup %>% 
    filter(timespan >= min_t,
           numbersurvs >= min_surv) 
  
  # centering and scaling
  dfr$y_s <- scale(dfr$y)
  dfr$year_s <- scale(dfr$year)
  
  # create factors
  dfr$culzone <- factor(dfr$culzone)
  dfr$country <- factor(dfr$country)
  dfr$demogr_grp <- factor(dfr$demogr_grp)
  dfr <- dfr %>% mutate(demogr_grp_country = factor(paste0(
    demogr_grp, "_", country
  )))
  
  dfr %>%
    summarise(
      n = nrow(.),
      countries = length(unique(country)),
      culzones = length(unique(culzone)),
      years = length(unique(year)),
      surveyed_country_years = length(unique(paste0(country, year)))
    ) %>% t() %>% print()
  
  return(dfr)
  
}

# create posterior distributions from hgam and census data for post-stratification
process_hgam <- function(m, df_census, ndraws){
  
  # create grid for predictions for all population groups
  newdat <- expand_grid(
    year_s = sort(unique(m$model$year_s)),
    country = levels(m$model$country),
    demogr_grp = levels(m$model$demogr_grp)
  ) %>%
    left_join(
      m$model %>% distinct(country, culzone),
      by = join_by(country)
    ) %>%
    mutate(demogr_grp_country = factor(paste0(
      demogr_grp, "_", country
    )))
  
  # get a posterior expectation
  df_pr <- gratia::fitted_samples(m,
                                  data = newdat,
                                  seed = 1,
                                  method = "gaussian",
                                  n = ndraws,
                                  mvn_method = "mgcv",
                                  unconditional = TRUE,
                                  freq = FALSE)
  
  # merge it to the data, and backtransform outcome and year
  df_pr <- df_pr %>%
    left_join(
      newdat %>%
        mutate(.row = 1:nrow(.)),
      by = join_by(.row)
    ) %>%
    mutate(
      year  = round(unscale(year_s,  center = attr(m$model$year_s, "scaled:center"), scale = attr(m$model$year_s, "scaled:scale"))),
      e_grb =       unscale(.fitted, center = attr(m$model$y_s,    "scaled:center"), scale = attr(m$model$y_s,    "scaled:scale"))
    ) %>%
    select(-year_s, -.row, -.fitted)
  
  # merge with post-stratification weights (which add to 1 per country-year)
  df_pr <- df_pr %>%
    left_join(
      df_census %>%
        transmute(country,
                  year,
                  demogr_grp = paste0(sex, ": ", age_grp),
                  sex, age_grp,
                  population, pr_stratum),
      by = join_by(country, year, demogr_grp)
    )
  
  return(df_pr)
  
}

# compute posterior means over time for different groupings
get_time_means_by_group <- function(df_est, df, ndraws){
  
  # global pop average
  df_glob <- df_est %>%
    group_by(year, .draw) %>%
    summarise(est = weighted.mean(est, population),
              population = sum(population) / ndraws, # ndraws counts repeated entries of countries here
              .groups = "drop")  %>%
    group_by(year) %>%
    summarise(m = mean(est),
              population = sum(population),
              lwr = quantile(est, 0.05),
              upr = quantile(est, 0.95))
  
  # average country
  df_avgc <- df_est %>% 
    group_by(year, .draw) %>% 
    summarise(est = mean(est), .groups = "drop_last") %>% 
    summarise(m = mean(est),
              lwr = quantile(est, 0.05),
              upr = quantile(est, 0.95),
              .groups = "drop")
  
  # cultural zones
  df_culz <- df_est %>%
    group_by(year, culzone, .draw) %>%
    summarise(est = weighted.mean(est, population),
              population = sum(population) / ndraws,
              .groups = "drop")  %>%
    group_by(year, culzone) %>%
    summarise(m = mean(est),
              population = sum(population),
              lwr = quantile(est, 0.05),
              upr = quantile(est, 0.95),
              .groups = "drop")
  
  # hdi
  df_hdi <- df_est %>%
    left_join(
      df %>% 
        select(country, hdicode) %>% 
        distinct(),
      by = join_by(country)
    ) %>%
    group_by(year, hdicode, .draw) %>%
    summarise(est = weighted.mean(est, population),
              population = sum(population) / ndraws,
              .groups = "drop")  %>%
    group_by(year, hdicode) %>%
    summarise(m = mean(est),
              population = sum(population),
              lwr = quantile(est, 0.05),
              upr = quantile(est, 0.95),
              .groups = "drop") %>% 
    mutate(hdicode = factor(hdicode, ordered = TRUE, levels = c("Low", "Medium", "High", "Very High")))
  
  res <- list(
    "df_glob" = df_glob,
    "df_avgc" = df_avgc,
    "df_culz" = df_culz,
    "df_hdi" = df_hdi
  )
  
  return(res)
  
}

# plot posteriors of time trends for different groupings
plot_time_means <- function(df_glob, df_avgc, df_culz, df_hdi){
  
  # global
  p_glob <- ggplot() +
    geom_line(data = df_glob, aes(x = year, y = m, linewidth = population), 
              col = "grey50", lty = 1, lineend = "round") +
    geom_ribbon(data = df_glob, aes(x = year, ymin = lwr, ymax = upr),
                fill = "grey50", alpha = 0.15) +
    geom_line(data = df_avgc, aes(x = year, y = m), 
              col = "black", lwd = 1, lty = 2, lineend = "round") +
    geom_ribbon(data = df_avgc, aes(x = year, ymin = lwr, ymax = upr),
                fill = "black", alpha = 0.15) +
    geom_dl(data = df_glob, aes(x = year, y = m), label = "Population avg.",
            method = list(dl.trans(x = x + 0.3), cex = 0.9, "last.qp")) +
    geom_dl(data = df_avgc, aes(x = year, y = m), label = "Avg. country",
            method = list(dl.trans(x = x + 0.3), cex = 0.9, "last.qp")) +
    theme_classic(14) +
    theme(plot.subtitle = element_text(face = "bold"),
          text = element_text(family = fontfam),
          legend.position = "inside",
          legend.key.spacing.y = unit(0, "in"),
          legend.text = element_text(size = 9),
          legend.margin = margin(0.15,0.15,0.15,0.15, unit ="cm"),
          legend.background = element_rect(color = "black"),
          legend.title = element_text(margin = margin(b = 1), size = 11),
          legend.position.inside = c(0.5,0.85)) +
    scale_x_continuous(breaks = seq(1995,2020,5), limits = c(1995, 2032)) +
    labs(x = "Year",
         y = "Mean GRB",
         subtitle= "Globally",
         linewidth = "Adult population") +
    scale_y_continuous(limits = c(0,100))+ 
    coord_cartesian(ylim = c(30,90)) +
    scale_linewidth_continuous(range = c(0.5,6),
                               limits = c(1*1e7, 5*1e9),
                               breaks = c(1*1e9, 3*1e9, 5*1e9),
                               labels = c("1 Billion", "3 Billion", "5 Billion")) +
    guides(linewidth = guide_legend(
      nrow = 1
    ))
  
  # cultural zones
  p_culz <- df_culz %>% 
    ggplot() + 
    geom_line(aes(x = year, y = m, col = culzone, group = culzone, linewidth= population),
              , lineend = "round") +
    geom_ribbon(aes(x = year, ymin = lwr, ymax = upr, fill = culzone), alpha = 0.15) +
    theme_classic(14) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none",
          plot.subtitle = element_text(face = "bold"),
          axis.title.y = element_blank(),
          text = element_text(family = fontfam)) +
    labs(x = "Year",
         subtitle = "by Cultural Zone",
         y = "Mean GRB") +
    scale_x_continuous(breaks = seq(1995,2020,5), limits = c(1995, 2032)) +
    scale_y_continuous(limits = c(0,100)) + 
    scale_color_manual(values = culzcols) +
    scale_fill_manual(values  = culzcols) +
    geom_dl(aes(x = year, y = m, label = culzone, col = culzone),
            method = list(dl.trans(x = x + 0.1), cex = 0.75, list("last.qp")))+ 
    coord_cartesian(ylim = c(30,90)) +
    scale_linewidth_continuous(range = c(0.5,6),
                               limits = c(1*1e7, 5*1e9))
  
  # hdi
  p_hdi <- df_hdi %>% 
    ggplot() + 
    geom_line(aes(x = year, y = m, col = hdicode, group = hdicode, linewidth= population),
              lineend = "round") +
    geom_ribbon(aes(x = year, ymin = lwr, ymax = upr, fill = hdicode), alpha = 0.15) +
    theme_classic(14) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none",
          plot.subtitle = element_text(face = "bold"),
          axis.title.y = element_blank(),
          text = element_text(family = fontfam)) +
    labs(x = "Year",
         subtitle = "by Human Development Index (2009)",
         y = "Mean GRB") +
    scale_x_continuous(breaks = seq(1995,2020,5), limits = c(1995, 2032)) +
    scale_y_continuous(limits = c(0,100)) + 
    scale_color_manual(values = hdicols) +
    scale_fill_manual(values  = hdicols) +
    geom_dl(aes(x = year, y = m, label = hdicode, col = hdicode),
            method = list(dl.trans(x = x + 0.15), cex = 1, list("last.qp")))+ 
    coord_cartesian(ylim = c(30,90)) +
    scale_linewidth_continuous(range = c(0.5,6),
                               limits = c(1*1e7, 5*1e9))
  
  p <- (p_glob + p_culz + p_hdi)
  return(p)
  
}

# fit hierarchical model for fertility
model_fert_hierarchical <- function(data, outcome, geq_var, covariates, with_age){
  
  # prepare model equation and data
  form          <- reformulate(c(geq_var, covariates), outcome)
  data_complete <- data %>% drop_na(all_of(c(outcome, geq_var, covariates)))
  
  # create model matrix with country-centered predictors included
  X <- Xs <- bind_cols(model.matrix(form, data =  data_complete), data_complete %>% select(country)) %>% 
    select(-`(Intercept)`) %>% 
    mutate(across(where(is.numeric), function(x) if(all(x %in% c(0,1))) x else scale(x))) %>% # scale all continuous pred except age
    # children_nr as outcome variable is NOT scaled
    mutate(across(.cols = everything(), mean, .names = "{.col}_mean"), .by = "country") %>% 
    select(-country)
  
  colnames(X)[!grepl("_mean", colnames(X))] <- paste0(colnames(X)[!grepl("_mean", colnames(X))], "_indiv")
  # create deviations: subtract mean from individual values
  X[, grepl("_indiv", colnames(X))] <- X[, grepl("_indiv", colnames(X))] - X[, grepl("_mean", colnames(X))]
  
  # avoid matrix columns
  X <- apply(X, 2, as.numeric)
  data_complete <- bind_cols(data_complete, X)
  
  # fit model
  form_full <- reformulate(c(paste0("`", colnames(X), "`"), if(isTRUE(with_age)) "age" else NULL, "(1 + geq_indiv |country)"), outcome)
  
  m <- lmer(form_full, data = data_complete, REML = TRUE,
            control = lmerControl(optimizer = "bobyqa"))
  res <- list(m=m, Xs=Xs)
  return(res)
}

# function specialized to buil the specific regression table for cross
# sectional fertility models
build_cross_regrtable <- function(mods){
  
  stopifnot(length(mods) == 3)
  
  # collect model results
  names(mods) <- 1:3
  restab <- map_dfr(mods, .id = "model", function(x){
    bind_cols(
      bind_rows(
        tidy(lmerTest::as_lmerModLmerTest(x), effects = "fixed"),
        insight::get_variance(x) %>% 
          as_tibble() %>% 
          select(var.intercept, var.slope, var.residual) %>% 
          pivot_longer(cols = contains("var"), names_to = "term", values_to = "estimate") %>% 
          mutate(estimate = sqrt(estimate)) %>% 
          mutate(term = paste0(term, "_sd"))
      ),
      bind_cols(glance(x),
                ncountries = summary(x)$ngrps)
    ) %>% 
      mutate(cond_r2 = suppressWarnings(performance::r2_nakagawa(x)$R2_conditional)) %>% 
      select(-effect)
  })
  
  # build a nice regression table
  restab_for_export <- restab %>% 
    mutate(
      predgroup = 
        case_when(
          grepl("education", term) ~ "Education",
          grepl("income_lvl", term) ~ "Income level",
          grepl("townsize", term) ~ "Town size",
          grepl("religiosity", term) ~ "Religiosity",
          term == "age" ~ "Age",
          grepl("_sd", term) ~ "Variance components",
          TRUE ~ " "
        ),
      modlevel = 
        case_when(
          grepl("_mean", term) ~ "between",
          grepl("_indiv", term) ~ "within",
          term == "age" ~ "pooled",
          TRUE ~ " "
        ),
      term = str_remove_all(term, "`"),
      term = str_remove_all(term, "education"),
      term = str_remove_all(term, "income_lvl"),
      term = str_remove_all(term, "townsize"),
      term = str_remove_all(term, "_mean"),
      term = str_remove_all(term, "_indiv"),
      term = case_when(
        term == "(Intercept)" ~ "Intercept",
        term == "var.intercept_sd" ~ "Country intercept SD",
        term == "var.slope_sd" ~ "Country slope SD",
        term == "var.residual_sd" ~ "Residual SD",
        term == "age" ~ "Age",
        term == "religiosity" ~ "Religiosity",
        grepl("geq", term) ~ "GRB",
        TRUE ~ term
      ),
      estimate_se = paste0(format(round(estimate, 3),nsmall = 3,trim = FALSE), 
                           ifelse(!is.na(std.error),
                                  paste0(" (", format(round(std.error, 3),nsmall = 3,trim = TRUE), ")"),
                                  "")
      ),
      p.value = ifelse(is.na(p.value), " ", 
                       ifelse(p.value < 1e-3, "<0.001", 
                              format(round(p.value, 3), nsmall = 3, trim = TRUE)))
    ) %>% 
    select(
      model, term, estimate_se, p.value, predgroup, modlevel
    ) %>% 
    pivot_wider(names_from = model, values_from = c(estimate_se, p.value)) %>% 
    bind_rows(
      tibble(term = c("N obs.", "N countries", "Cond. R2"),
             predgroup = c("Model information"),
             modlevel = c(" "),
             estimate_se_1 = c(restab %>% filter(model=="1") %>% pull(nobs) %>% unique, 
                               restab %>% filter(model=="1") %>% pull(ncountries) %>% unique,
                               restab %>% filter(model=="1") %>% pull(cond_r2) %>% unique),
             estimate_se_2 = c(restab %>% filter(model=="2") %>% pull(nobs) %>% unique, 
                               restab %>% filter(model=="2") %>% pull(ncountries) %>% unique,
                               restab %>% filter(model=="2") %>% pull(cond_r2) %>% unique),
             estimate_se_3 = c(restab %>% filter(model=="3") %>% pull(nobs) %>% unique, 
                               restab %>% filter(model=="3") %>% pull(ncountries) %>% unique,
                               restab %>% filter(model=="3") %>% pull(cond_r2) %>% unique)
      ) %>% mutate(across(where(is.numeric), 
                          ~ifelse(.x > 1,
                                  format(round(.x, 0), big.mark = ",", 
                                         nsmall = 0,
                                         scientific=FALSE),
                                  format(round(.x, 3), nsmall = 3))))
    ) %>% 
    mutate(across(.cols = everything(), ~ifelse(is.na(.x), " ", .x))) %>% 
    # ordering for display
    mutate(predgroup = 
             factor(predgroup, ordered = TRUE, levels = 
                      c(" ",
                        "Age",
                        "Education",
                        "Income level",
                        "Town size",
                        "Religiosity",
                        "Variance components",
                        "Model information")),
           modlevel = 
             factor(modlevel, ordered = TRUE, levels = 
                      c(" ", "within", "between", "pooled")),
           term = 
             factor(term, ordered = TRUE, levels = 
                      c("Intercept",
                        "GRB",
                        "Age",
                        "Lower secondary",
                        "Upper secondary",
                        "Tertiary",
                        "Medium",
                        "High",
                        "5,000-20,000",
                        "20,000-100,000",
                        "100,000-500,000",
                        ">500,000",
                        "Religiosity",
                        "Country intercept SD",
                        "Country slope SD",
                        "Residual SD",
                        "N obs.", 
                        "N countries",
                        "Cond. R2"))) %>% 
    arrange(predgroup, modlevel, term) %>% 
    relocate(predgroup, term, modlevel, contains("_1"), contains("_2"), contains("_3"))
  
  restab_for_export  %>% 
    select(-predgroup) %>% 
    kable(
      align = "l",
      col.names = c("Term", "Effect level", rep(c("Estimate (Std. error)", "p-value"), 3)),
      caption = "<p align='left';>Results from linear mixed-effects models for number of children among women aged 40 - 49 years with country random intercepts. 
    GRB and religiosity have been scaled to have mean zero and SD of one. 
    Reference categories are 'Primary or less' (Education), 'Low' (Income level) and '<5,000' (Town size).
    The effect of GRB is modeled both within and between countries by mean-centering the variable
    within countries and including the country mean as additional predictor. Where appropriate, this was also
    done for other predictors (it was not appropriate for Age, since the data was already restricted to the same age range between all countries).</p>"
    ) %>% 
    kable_styling() %>% 
    add_header_above(c(" " = 2, "Model 1 (only GRB)" = 2, "Model 2 (with demographics)" = 2, "Model 3 (with religiosity)" = 2)) %>% 
    group_rows(index = table(restab_for_export$predgroup)) %>% 
    row_spec(which(grepl("GRB", restab_for_export$term)), background = "Gainsboro")
  
  
}
