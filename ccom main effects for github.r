## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)

## ----data_prep-----------------------------------------------------------
## setup
library(lme4)
library(lmerTest)
library(glmmTMB)
library(mitml)
library(tidyverse)
library(broom)
library(broom.mixed)
library(cowplot)
library(kableExtra)

ccom <- read_csv("~/Dropbox/PhD/Research/Couples Communication/Food task/Data assembly/Data assembly final/clean/all2.csv")
ccom <- ccom %>%
  mutate(couple = as.factor(couple),
         food = factor(food, levels = c("Radish", "Cookie")),
         pfood = factor(pfood, levels = c("Radish", "Cookie")),
         Total.pinsBin = as.numeric(ccom$Total.pins > 0),
         Total.pins_trunc = ifelse(Total.pins > 0,
                                   Total.pins,
                                   as.numeric(NA)))

couple_ratings <- read_csv("~/Dropbox/PhD/Research/Couples Communication/Food task/Data assembly/Data assembly final/clean/couple_ratings.csv")
couple_ratings$couple <- as.factor(couple_ratings$couple)

indiv_data <- ccom %>%
  select(ID,
         couple,
         food,
         pfood,
         DEBQemotr,
         # beh ratings
         pos.factor,
         p.pos.factor,
         neg.factor,
         p.neg.factor,
         a.Unconstructive.,
         a.p.Unconstructive.,
         a.Positive.mood.,
         a.p.Positive.mood.,
         # self report
         PFQ,
         CTI) %>%
  mutate_at(vars(-ID, -couple, -food, -pfood),
            funs(scale(.))) %>%
  gather(variable, value, -ID, -couple, -food, -pfood, -DEBQemotr)

couple_data <- couple_ratings %>% select(ID,
                                         couple,
                                         ndep,
                                         couple.debq.emotr,
                                         Mutual.avoidance.,
                                         p.Mutual.avoidance.,
                                         Negative.Reciprocity.,
                                         p.Negative.Reciprocity.,
                                         Positive.Reciprocity.,
                                         p.Positive.Reciprocity.,
                                         Vulnerability.empathy.support.,
                                         p.Vulnerability.empathy.support.) %>%
   mutate_at(vars(-ID, -couple, -ndep),
            funs(scale(-1 / .^(1/2)))) %>%
  gather(variable, value, -ID, -couple, -ndep, -couple.debq.emotr)

# read in indiv imps---------------------------
# for affect change analyses

blimps_indiv <- read_csv("/Users/Ben/Dropbox/PhD/Research/Couples Communication/Paper draft/tidied analyses/imputation/individual/ccom-prq-out-debq2.csv",
                   na = "999", col_names = c("imp",
                                             "couple",
                                             "ID",
                                             "PFQ",
                                             "Total.pins",
                                             "pos.factor",
                                             "AG",
                                             "DEBQrest",
                                             "DEBQexte",
                                             "prePosA",
                                             "preNegA",
                                             "postPosA",
                                             "postNegA",
                                             "satisfied",
                                             "content",
                                             "happy",
                                             "committed",
                                             "dedicated",
                                             "devoted",
                                             "intimate",
                                             "close",
                                             "connected",
                                             "trust",
                                             "count.on",
                                             "dependable",
                                             "passionate",
                                             "lustful",
                                             "sexual",
                                             "love",
                                             "adore",
                                             "cherish",
                                             "IMSr",
                                             "food",
                                             "pfood",
                                             "DEBQemot",
                                             "prq.mean",
                                             "prq.mean_p",
                                             "ims_p",
                                             "foodxpfood",
                                             "prqxfood",
                                             "prqxpfood",
                                             "imsxfood",
                                             "imsxpfood",
                                             "prq_pxfood",
                                             "prq_pxpfood",
                                             "ims_pxfood",
                                             "ims_pxpfood"))
blimps_indiv <- blimps_indiv %>% filter(couple != 0)

# split into separate files
blimplist_indiv <- split(blimps_indiv, blimps_indiv$imp)

# storing prq item names
prq.names <- c("satisfied", "content", "happy",
               "committed", "dedicated", "devoted",
               "intimate", "close", "connected",
               "trust", "count.on", "dependable",
               "passionate", "lustful", "sexual",
               "love", "adore", "cherish")

blimplist_indiv <- map(blimplist_indiv,
                 # sum prq items
                 # calculate affect change
                 ~mutate(rowwise(.),
                         prq = mean(c(!!! syms(prq.names))),
                         posAchange = postPosA - prePosA,
                         negAchange = postNegA - preNegA) %>%
                   ungroup() %>%
                   rename(ims = IMSr) %>%
                   # center
                   mutate_at(vars(prq, ims, posAchange, negAchange),
                      funs(c = scale(.))) %>%
                   group_by(couple) %>%
                   mutate(p.prq_c = rev(prq_c),
                          p.ims_c = rev(ims_c)) %>%
                   ungroup()
                 )
# https://stackoverflow.com/a/48846441/
# explains this weird evaluation to get vector of column names

blimplist_indiv <- as.mitml.list(blimplist_indiv)

print_kable <- function(table, ...) {
  kable(table, ...) %>% 
    kable_styling(full_width = FALSE)
}

## ----main_effects--------------------------------------------------------
# individual models, except pins and affect-change---------------
indiv_me_models <- indiv_data %>%
  nest(-variable) %>%
  mutate(initial_fit = map(data,
                           ~ lmer(formula = value ~ food*pfood + DEBQemotr + (1 | couple),
                                  data = .x)),
         initial_tidied = map(initial_fit, tidy))

indiv_me_initial_sig <- indiv_me_models %>%
  unnest(initial_tidied) %>%
  select(variable, term, p.value) %>%
  filter(term %in% c("foodCookie",
                     "pfoodCookie",
                     "foodCookie:pfoodCookie")) %>%
  mutate(sig = ifelse(p.value < .05, 1, 0)) %>%
  mutate(term = str_sub(term, end = -7),
         term = paste(term, "ini_fit_sig", sep = "_")) %>%
  select(-p.value) %>%
  spread(term, sig)

indiv_me_models <- full_join(indiv_me_models, indiv_me_initial_sig)

# update models
# for models with only one significant predictor, drop the other
indiv_me_models <- indiv_me_models %>%
  mutate(fit = case_when(
    # if interaction term is significant, keep original
    "foodCookie:pfoodCookie_ini_fit_sig" == 1 ~ initial_fit,
    # if only one main effect is significant, drop the other
    food_ini_fit_sig == 1 & pfood_ini_fit_sig == 0 ~
      map(data,
          ~ lmer(formula = value ~ food + DEBQemotr + (1 | couple),
                 data = .x)),
    food_ini_fit_sig == 0 & pfood_ini_fit_sig == 1 ~
      map(data,
          ~ lmer(formula = value ~ pfood + DEBQemotr + (1 | couple),
                 data = .x)),
    # else keep initial model
    TRUE ~ initial_fit),
    tidied = map(fit, tidy))

indiv_me_fit <- indiv_me_models %>%
  unnest(tidied) %>%
  filter(effect == "fixed" & ! term %in% c("(Intercept)", "DEBQemotr")) %>%
  select(-ends_with("ini_fit_sig"), -group, -effect)



## ----pin_models----------------------------------------------------------
# pin models---------------------

pin_count_me_model <- glmmTMB(formula = Total.pins_trunc ~ food +
                           pfood + DEBQemotr + (1 | couple),
                         data = ccom,
                         family = "truncated_poisson")
pin_count_me_fit <- pin_count_me_model %>%
  tidy() %>%
  mutate(variable = "pin_count")

pin_bin_me_model <- glmer(Total.pinsBin ~ food + pfood + DEBQemotr + (1 | couple),
                          data=ccom, family= binomial)

pin_bin_me_fit <- pin_bin_me_model %>%
  tidy() %>%
  mutate(variable = "pin_bin")

pin_me_fit <- bind_rows(pin_count_me_fit, pin_bin_me_fit) %>%
  filter(effect == "fixed" & ! term %in% c("(Intercept)", "DEBQemotr")) %>%
  select(-component, -group, -effect)

## ----affect_change-------------------------------------------------------
# tidying function
tidy_mitml <- function(x, ...) {
  x <- x$estimates
  x_tidied <- as_tibble(x, rownames = "term")
  x_tidied <- x_tidied %>% rename(estimate = Estimate,
                                  std.error = Std.Error,
                                  statistic = t.value,
                                  p.value = "P(>|t|)",
                                  riv = RIV,
                                  fmi = FMI)
  x_tidied
  }

# positive affect change-----------------
posa_change_model <- with(blimplist_indiv, lmer(formula = posAchange_c ~ DEBQemot + food*pfood + (1 | couple)))
model.df <- summary(posa_change_model[[1]]) %>%
  .$coefficients %>%
  .[ , "df"] %>%
  as.vector()
posa_change_fit <- testEstimates(posa_change_model, df = model.df) %>%
  tidy_mitml()
posa_change_fit

# interaction non-sig, dropping
posa_change_model <- with(blimplist_indiv, lmer(formula = posAchange_c ~ DEBQemot + food + pfood + (1 | couple)))
model.df <- summary(posa_change_model[[1]]) %>%
  .$coefficients %>%
  .[ , "df"] %>%
  as.vector()
posa_change_fit <- testEstimates(posa_change_model, df = model.df) %>%
  tidy_mitml() %>%
  select(-riv, -fmi) %>%
  mutate(variable = "pos_a_change") %>%
  filter(! term %in% c("(Intercept)", "DEBQemot"))
posa_change_fit

# negative affect change-----------------
nega_change_model <- with(blimplist_indiv, lmer(formula = negAchange_c ~ DEBQemot + food*pfood + (1 | couple)))
model.df <- summary(nega_change_model[[1]]) %>%
  .$coefficients %>%
  .[ , "df"] %>%
  as.vector()
nega_change_fit <- testEstimates(nega_change_model, df = model.df) %>%
  tidy_mitml()
nega_change_fit

# interaction non-sig, dropping
nega_change_model <- with(blimplist_indiv, lmer(formula = negAchange_c ~ DEBQemot + food + pfood + (1 | couple)))
model.df <- summary(nega_change_model[[1]]) %>%
  .$coefficients %>%
  .[ , "df"] %>%
  as.vector()
nega_change_fit <- testEstimates(nega_change_model, df = model.df) %>%
  tidy_mitml() %>%
  select(-riv, -fmi) %>%
  mutate(variable = "neg_a_change") %>%
  filter(! term %in% c("(Intercept)", "DEBQemot"))
nega_change_fit


## ----couple_ratings------------------------------------------------------
# couple ratings--------------------------
couple_me_models <- couple_data %>%
  nest(-variable) %>%
  mutate(fit = map(data,
                   ~ lm(formula = value ~ ndep + couple.debq.emotr,
                                  data = .x)),
         tidied = map(fit, tidy))

couple_me_fit <- couple_me_models %>%
  unnest(tidied) %>%
  filter(! term %in% c("(Intercept)", "DEBQemotr"))

## ----combining_results---------------------------------------------------
# binding pin models with other individual models
indiv_me_fit <- bind_rows(indiv_me_fit, pin_me_fit,
                          posa_change_fit, nega_change_fit)

# binding couple ratings with other main effects results
me_fit <- bind_rows(indiv_me_fit, couple_me_fit)

# keeping the beta and pvalue columns only
me_fit <- me_fit %>%
  select(variable, term, estimate, p.value)

# extracting results for manuscript----------------

# behavior ratings---------------------------------
beh_ratings_fit <- me_fit %>%
  filter(variable %in% c("pos.factor",
                         "p.pos.factor",
                         "neg.factor",
                         "p.neg.factor",
                         "a.Unconstructive.",
                         "a.p.Unconstructive.",
                         "a.Positive.mood.",
                         "a.p.Positive.mood.",
                         "Mutual.avoidance.",
                         "p.Mutual.avoidance.",
                         "Negative.Reciprocity.",
                         "p.Negative.Reciprocity.",
                         "Positive.Reciprocity.",
                         "p.Positive.Reciprocity.",
                         "Vulnerability.empathy.support.",
                         "p.Vulnerability.empathy.support."))
filter(beh_ratings_fit, p.value < .05)

# next lowest p value
beh_ratings_fit %>%
  arrange(p.value) %>%
  slice(1:2)

# affect change------------------------------------
me_fit %>%
  filter(variable %in% c("pos_a_change", "neg_a_change"))

# positive feelings for partner
me_fit %>%
  filter(variable == "PFQ")

# aggressive impulses
me_fit %>%
  filter(variable == "CTI")

# pins
me_fit %>%
  filter(variable %in% c("pin_count", "pin_bin"))

print_kable(me_fit)

