---
title: "ccom main effects"
author: "Ben Shulman"
output:
  html_document:
    keep_md: true
    toc: true
---



# prep


```r
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
```

# fitting individual models

## individual variables less pin and affect-change

I fit the pin models separately because pin-insertion is a count variable.


```r
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
```

## pin variables

Fitting these models separately to handle the pin variable's count distribution with a hurdle model.


```r
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
```

## affect change variables

Fitting these models separately because I used multiple imputation to handle the missingness on participants' affect ratings.


```r
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
```

```
## # A tibble: 5 x 8
##   term        estimate std.error statistic    df p.value    riv    fmi
##   <chr>          <dbl>     <dbl>     <dbl> <dbl>   <dbl>  <dbl>  <dbl>
## 1 (Intercept)  1.04      0.827       1.26  106.    0.210 0.173  0.149 
## 2 DEBQemot     0.00167   0.00689     0.242  72.9   0.809 0.377  0.280 
## 3 food        -0.229     0.530      -0.432 122.    0.667 0.0970 0.0892
## 4 pfood       -0.248     0.553      -0.449  97.7   0.654 0.213  0.179 
## 5 food:pfood  -0.163     0.337      -0.484 116.    0.629 0.125  0.112
```

```r
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
```

```
## # A tibble: 2 x 7
##   term  estimate std.error statistic    df p.value variable    
##   <chr>    <dbl>     <dbl>     <dbl> <dbl>   <dbl> <chr>       
## 1 food    -0.475     0.181     -2.62  89.5 0.0102  pos_a_change
## 2 pfood   -0.493     0.180     -2.75  84.9 0.00737 pos_a_change
```

```r
# negative affect change-----------------
nega_change_model <- with(blimplist_indiv, lmer(formula = negAchange_c ~ DEBQemot + food*pfood + (1 | couple)))
model.df <- summary(nega_change_model[[1]]) %>%
  .$coefficients %>%
  .[ , "df"] %>%
  as.vector()
nega_change_fit <- testEstimates(nega_change_model, df = model.df) %>%
  tidy_mitml()
nega_change_fit
```

```
## # A tibble: 5 x 8
##   term        estimate std.error statistic    df p.value    riv    fmi
##   <chr>          <dbl>     <dbl>     <dbl> <dbl>   <dbl>  <dbl>  <dbl>
## 1 (Intercept)  1.08      0.869       1.25   64.1  0.217  0.0699 0.0658
## 2 DEBQemot     0.00648   0.00807     0.803  48.5  0.426  0.682  0.415 
## 3 food        -0.769     0.575      -1.34   70.4  0.185  0.0746 0.0698
## 4 pfood       -1.01      0.576      -1.75   67.5  0.0852 0.0948 0.0873
## 5 food:pfood   0.593     0.364       1.63   63.1  0.108  0.0776 0.0725
```

```r
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
```

```
## # A tibble: 2 x 7
##   term  estimate std.error statistic    df p.value variable    
##   <chr>    <dbl>     <dbl>     <dbl> <dbl>   <dbl> <chr>       
## 1 food     0.125     0.194     0.647  86.6   0.519 neg_a_change
## 2 pfood   -0.118     0.176    -0.667 122.    0.506 neg_a_change
```

# couple ratings

I model the couple couple-level ratings separately because these don't require multi-level models; they can be fit using OLS.


```r
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
```

# extracting results for manuscript


```r
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
```

```
## # A tibble: 1 x 4
##   variable                term  estimate p.value
##   <chr>                   <chr>    <dbl>   <dbl>
## 1 p.Positive.Reciprocity. ndep    -0.441 0.00839
```

```r
# next lowest p value
beh_ratings_fit %>%
  arrange(p.value) %>%
  slice(1:2)
```

```
## # A tibble: 2 x 4
##   variable                term       estimate p.value
##   <chr>                   <chr>         <dbl>   <dbl>
## 1 p.Positive.Reciprocity. ndep         -0.441 0.00839
## 2 p.pos.factor            foodCookie   -0.422 0.123
```

```r
# affect change------------------------------------
me_fit %>%
  filter(variable %in% c("pos_a_change", "neg_a_change"))
```

```
## # A tibble: 4 x 4
##   variable     term  estimate p.value
##   <chr>        <chr>    <dbl>   <dbl>
## 1 pos_a_change food    -0.475 0.0102 
## 2 pos_a_change pfood   -0.493 0.00737
## 3 neg_a_change food     0.125 0.519  
## 4 neg_a_change pfood   -0.118 0.506
```

```r
# positive feelings for partner
me_fit %>%
  filter(variable == "PFQ")
```

```
## # A tibble: 3 x 4
##   variable term                   estimate p.value
##   <chr>    <chr>                     <dbl>   <dbl>
## 1 PFQ      foodCookie              -0.283   0.310 
## 2 PFQ      pfoodCookie             -0.480   0.0851
## 3 PFQ      foodCookie:pfoodCookie   0.0977  0.823
```

```r
# aggressive impulses
me_fit %>%
  filter(variable == "CTI")
```

```
## # A tibble: 3 x 4
##   variable term                   estimate p.value
##   <chr>    <chr>                     <dbl>   <dbl>
## 1 CTI      foodCookie              -0.303    0.197
## 2 CTI      pfoodCookie              0.0456   0.844
## 3 CTI      foodCookie:pfoodCookie   0.352    0.286
```

```r
# pins
me_fit %>%
  filter(variable %in% c("pin_count", "pin_bin"))
```

```
## # A tibble: 4 x 4
##   variable  term        estimate p.value
##   <chr>     <chr>          <dbl>   <dbl>
## 1 pin_count foodCookie   -0.422   0.115 
## 2 pin_count pfoodCookie   0.225   0.397 
## 3 pin_bin   foodCookie    0.657   0.0902
## 4 pin_bin   pfoodCookie  -0.0658  0.862
```

```r
print_kable(me_fit)
```

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> variable </th>
   <th style="text-align:left;"> term </th>
   <th style="text-align:right;"> estimate </th>
   <th style="text-align:right;"> p.value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> pos.factor </td>
   <td style="text-align:left;"> foodCookie </td>
   <td style="text-align:right;"> 0.0691277 </td>
   <td style="text-align:right;"> 0.8055033 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pos.factor </td>
   <td style="text-align:left;"> pfoodCookie </td>
   <td style="text-align:right;"> -0.1703655 </td>
   <td style="text-align:right;"> 0.5440298 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pos.factor </td>
   <td style="text-align:left;"> foodCookie:pfoodCookie </td>
   <td style="text-align:right;"> 0.0352467 </td>
   <td style="text-align:right;"> 0.9390749 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> p.pos.factor </td>
   <td style="text-align:left;"> foodCookie </td>
   <td style="text-align:right;"> -0.4221666 </td>
   <td style="text-align:right;"> 0.1229511 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> p.pos.factor </td>
   <td style="text-align:left;"> pfoodCookie </td>
   <td style="text-align:right;"> -0.3862752 </td>
   <td style="text-align:right;"> 0.1560564 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> p.pos.factor </td>
   <td style="text-align:left;"> foodCookie:pfoodCookie </td>
   <td style="text-align:right;"> 0.3516026 </td>
   <td style="text-align:right;"> 0.4178435 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> neg.factor </td>
   <td style="text-align:left;"> foodCookie </td>
   <td style="text-align:right;"> -0.1636154 </td>
   <td style="text-align:right;"> 0.5491147 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> neg.factor </td>
   <td style="text-align:left;"> pfoodCookie </td>
   <td style="text-align:right;"> 0.0733768 </td>
   <td style="text-align:right;"> 0.7878043 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> neg.factor </td>
   <td style="text-align:left;"> foodCookie:pfoodCookie </td>
   <td style="text-align:right;"> 0.4182144 </td>
   <td style="text-align:right;"> 0.3402856 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> p.neg.factor </td>
   <td style="text-align:left;"> foodCookie </td>
   <td style="text-align:right;"> 0.2489451 </td>
   <td style="text-align:right;"> 0.4281399 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> p.neg.factor </td>
   <td style="text-align:left;"> pfoodCookie </td>
   <td style="text-align:right;"> 0.3852624 </td>
   <td style="text-align:right;"> 0.2206021 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> p.neg.factor </td>
   <td style="text-align:left;"> foodCookie:pfoodCookie </td>
   <td style="text-align:right;"> -0.2147111 </td>
   <td style="text-align:right;"> 0.6799660 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> a.Unconstructive. </td>
   <td style="text-align:left;"> foodCookie </td>
   <td style="text-align:right;"> -0.2460237 </td>
   <td style="text-align:right;"> 0.3891559 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> a.Unconstructive. </td>
   <td style="text-align:left;"> pfoodCookie </td>
   <td style="text-align:right;"> -0.1647115 </td>
   <td style="text-align:right;"> 0.5634684 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> a.Unconstructive. </td>
   <td style="text-align:left;"> foodCookie:pfoodCookie </td>
   <td style="text-align:right;"> 0.4411174 </td>
   <td style="text-align:right;"> 0.3203363 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> a.p.Unconstructive. </td>
   <td style="text-align:left;"> foodCookie </td>
   <td style="text-align:right;"> 0.0793759 </td>
   <td style="text-align:right;"> 0.7762909 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> a.p.Unconstructive. </td>
   <td style="text-align:left;"> pfoodCookie </td>
   <td style="text-align:right;"> 0.0736962 </td>
   <td style="text-align:right;"> 0.7909161 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> a.p.Unconstructive. </td>
   <td style="text-align:left;"> foodCookie:pfoodCookie </td>
   <td style="text-align:right;"> -0.0086233 </td>
   <td style="text-align:right;"> 0.9838699 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> a.Positive.mood. </td>
   <td style="text-align:left;"> foodCookie </td>
   <td style="text-align:right;"> -0.1687008 </td>
   <td style="text-align:right;"> 0.5465348 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> a.Positive.mood. </td>
   <td style="text-align:left;"> pfoodCookie </td>
   <td style="text-align:right;"> -0.0040843 </td>
   <td style="text-align:right;"> 0.9883299 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> a.Positive.mood. </td>
   <td style="text-align:left;"> foodCookie:pfoodCookie </td>
   <td style="text-align:right;"> 0.0670859 </td>
   <td style="text-align:right;"> 0.8847742 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> a.p.Positive.mood. </td>
   <td style="text-align:left;"> foodCookie </td>
   <td style="text-align:right;"> -0.0943199 </td>
   <td style="text-align:right;"> 0.7339943 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> a.p.Positive.mood. </td>
   <td style="text-align:left;"> pfoodCookie </td>
   <td style="text-align:right;"> -0.0775624 </td>
   <td style="text-align:right;"> 0.7793572 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> a.p.Positive.mood. </td>
   <td style="text-align:left;"> foodCookie:pfoodCookie </td>
   <td style="text-align:right;"> -0.1484663 </td>
   <td style="text-align:right;"> 0.7462254 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PFQ </td>
   <td style="text-align:left;"> foodCookie </td>
   <td style="text-align:right;"> -0.2825061 </td>
   <td style="text-align:right;"> 0.3098041 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PFQ </td>
   <td style="text-align:left;"> pfoodCookie </td>
   <td style="text-align:right;"> -0.4800941 </td>
   <td style="text-align:right;"> 0.0851267 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PFQ </td>
   <td style="text-align:left;"> foodCookie:pfoodCookie </td>
   <td style="text-align:right;"> 0.0977403 </td>
   <td style="text-align:right;"> 0.8231246 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CTI </td>
   <td style="text-align:left;"> foodCookie </td>
   <td style="text-align:right;"> -0.3030566 </td>
   <td style="text-align:right;"> 0.1970405 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CTI </td>
   <td style="text-align:left;"> pfoodCookie </td>
   <td style="text-align:right;"> 0.0455959 </td>
   <td style="text-align:right;"> 0.8444962 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CTI </td>
   <td style="text-align:left;"> foodCookie:pfoodCookie </td>
   <td style="text-align:right;"> 0.3516320 </td>
   <td style="text-align:right;"> 0.2861697 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pin_count </td>
   <td style="text-align:left;"> foodCookie </td>
   <td style="text-align:right;"> -0.4224814 </td>
   <td style="text-align:right;"> 0.1146920 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pin_count </td>
   <td style="text-align:left;"> pfoodCookie </td>
   <td style="text-align:right;"> 0.2248435 </td>
   <td style="text-align:right;"> 0.3973964 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pin_bin </td>
   <td style="text-align:left;"> foodCookie </td>
   <td style="text-align:right;"> 0.6571991 </td>
   <td style="text-align:right;"> 0.0901995 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pin_bin </td>
   <td style="text-align:left;"> pfoodCookie </td>
   <td style="text-align:right;"> -0.0658154 </td>
   <td style="text-align:right;"> 0.8623353 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pos_a_change </td>
   <td style="text-align:left;"> food </td>
   <td style="text-align:right;"> -0.4749164 </td>
   <td style="text-align:right;"> 0.0102191 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pos_a_change </td>
   <td style="text-align:left;"> pfood </td>
   <td style="text-align:right;"> -0.4931710 </td>
   <td style="text-align:right;"> 0.0073683 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> neg_a_change </td>
   <td style="text-align:left;"> food </td>
   <td style="text-align:right;"> 0.1254490 </td>
   <td style="text-align:right;"> 0.5193018 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> neg_a_change </td>
   <td style="text-align:left;"> pfood </td>
   <td style="text-align:right;"> -0.1175714 </td>
   <td style="text-align:right;"> 0.5060179 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Mutual.avoidance. </td>
   <td style="text-align:left;"> ndep </td>
   <td style="text-align:right;"> -0.0580805 </td>
   <td style="text-align:right;"> 0.7360520 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Mutual.avoidance. </td>
   <td style="text-align:left;"> couple.debq.emotr </td>
   <td style="text-align:right;"> 0.1104476 </td>
   <td style="text-align:right;"> 0.4093607 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> p.Mutual.avoidance. </td>
   <td style="text-align:left;"> ndep </td>
   <td style="text-align:right;"> -0.0815538 </td>
   <td style="text-align:right;"> 0.6327243 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> p.Mutual.avoidance. </td>
   <td style="text-align:left;"> couple.debq.emotr </td>
   <td style="text-align:right;"> 0.0296559 </td>
   <td style="text-align:right;"> 0.8206299 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Negative.Reciprocity. </td>
   <td style="text-align:left;"> ndep </td>
   <td style="text-align:right;"> 0.1498707 </td>
   <td style="text-align:right;"> 0.3737823 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Negative.Reciprocity. </td>
   <td style="text-align:left;"> couple.debq.emotr </td>
   <td style="text-align:right;"> -0.0307436 </td>
   <td style="text-align:right;"> 0.8132561 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> p.Negative.Reciprocity. </td>
   <td style="text-align:left;"> ndep </td>
   <td style="text-align:right;"> 0.0379143 </td>
   <td style="text-align:right;"> 0.7036227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> p.Negative.Reciprocity. </td>
   <td style="text-align:left;"> couple.debq.emotr </td>
   <td style="text-align:right;"> -0.0390528 </td>
   <td style="text-align:right;"> 0.6096356 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Positive.Reciprocity. </td>
   <td style="text-align:left;"> ndep </td>
   <td style="text-align:right;"> -0.1763492 </td>
   <td style="text-align:right;"> 0.3100357 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Positive.Reciprocity. </td>
   <td style="text-align:left;"> couple.debq.emotr </td>
   <td style="text-align:right;"> 0.0987189 </td>
   <td style="text-align:right;"> 0.4623681 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> p.Positive.Reciprocity. </td>
   <td style="text-align:left;"> ndep </td>
   <td style="text-align:right;"> -0.4410113 </td>
   <td style="text-align:right;"> 0.0083900 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> p.Positive.Reciprocity. </td>
   <td style="text-align:left;"> couple.debq.emotr </td>
   <td style="text-align:right;"> 0.1189843 </td>
   <td style="text-align:right;"> 0.3405977 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Vulnerability.empathy.support. </td>
   <td style="text-align:left;"> ndep </td>
   <td style="text-align:right;"> -0.0125050 </td>
   <td style="text-align:right;"> 0.9400215 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Vulnerability.empathy.support. </td>
   <td style="text-align:left;"> couple.debq.emotr </td>
   <td style="text-align:right;"> 0.0735823 </td>
   <td style="text-align:right;"> 0.5683483 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> p.Vulnerability.empathy.support. </td>
   <td style="text-align:left;"> ndep </td>
   <td style="text-align:right;"> -0.0857688 </td>
   <td style="text-align:right;"> 0.6405513 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> p.Vulnerability.empathy.support. </td>
   <td style="text-align:left;"> couple.debq.emotr </td>
   <td style="text-align:right;"> -0.1451520 </td>
   <td style="text-align:right;"> 0.3047429 </td>
  </tr>
</tbody>
</table>
