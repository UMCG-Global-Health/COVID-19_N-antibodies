library(tidyverse)
library(nnet)
library(splines)
library(ggplot2)

#Example code for visualization of splines for ct values with 3 knots

df <- df %>%
  filter(swab_pos %in% c('Positive swab during', 'Positive swab before and during')) %>% 
  group_by(pseudo_id) %>%
  mutate(
    closest_swab = days_first_pos_all[which.min(abs(days_first_pos_all - (n_ab_pos_date)))],  
    seropositivity_binary = ifelse(n_ab_class %in% c('increasing', 'de- and increasing') & 
                                     swab_pos %in% c('Positive swab during', 'Positive swab before and during')& 
                                     days_first_pos_all == closest_swab, 'both', 'missed by n'))


#exclusion based on exclusion criteria
df <- df %>%  
  filter(flag_zero_before_mult == 0 | is.na(flag_zero_before_mult)) %>%       # Exclude participants with s-ab infection before study period
  filter(flag_zero_mult == 0) %>%                                             # Exclude participants with s-ab infection during study period
  filter(dummy_within_60 == 1 | seropositivity_binary == 'missed by n') %>%   # Exclude those with more than 60 days between infection dates
  filter(!is.na(vaccination)) %>%   
  filter(!is.na(ct_value)) %>% 
  filter(!is.na(symptoms)) %>% 
  filter(type != 'Alpha epoch')

# Create spline terms for age and ct
splineterms_age <- data.frame(ns(df$age_at_visit, knots=quantile(df$age_at_visit, probs = c(0.28, 0.5, 0.72)), 
                                 Boundary.knots = quantile(df$age_at_visit,probs=c(0.05,0.95))))
splineterms_ct <- data.frame(ns(df$ct_value, knots= c(15, 29, 34)))
df <- cbind(df, splineterms_age, splineterms_ct)
names(df)[30:33] <- c("spline1_age","spline2_age", 'spline3_age', 'spline4_age')
names(df)[34:37] <- c("spline1_ct","spline2_ct", 'spline3_ct', 'spline4_ct')

# Logistic regression
model <- multinom(seropositivity_binary ~  sex+ spline1_age + spline2_age + spline3_age + spline4_age+
                    spline1_ct + spline2_ct +spline3_ct + spline4_ct+
                    ethnicity_wo + ever_hcw + ever_lthc + vaccination +
                    symptoms + type,
                  data = df)

X <- model.matrix(model)
X_df <- data.frame(X)

# Fixing covariates to get smooth curve
X_df$sexMale <- 1
X_df$spline1_age <- 0.28
X_df$spline2_age <- -0.09
X_df$spline3_age <- 0.26
X_df$spline4_age <- 0.003
X_df$ethnicity_woWhite <- 1
X_df$ever_hcwYes <- 0
X_df$ever_lthcYes <- 0
X_df$vaccination1.vaccination <- 1
X_df$vaccination2.vaccinations..3.months.and..6.months <- 0
X_df$vaccination2.vaccinations.6..months <- 0
X_df$vaccination3.4.vaccinations <- 0
X_df$vaccinationnot.vaccinated <- 0
X_df$symptomsYes.amn.ageu <- 1
X_df$symptomsYes.cgh.fev.not.amn.ageu <- 0
X_df$symptomsYes.other <- 0
X_df$typeBA.1.epoch <- 0
X <- as.matrix(X_df)

# Model coefficients, with expicit zero row added for reference category & transposed
betahat <- t(rbind(0, coef(model))) 

# Predictions on the link scale (with explicit zero for reference level included here, sometimes this level is dropped)
preds_link <- X %*% betahat 
colnames(preds_link) <- model$lev
write.pred <- cbind(df, preds_link)
write.pred <- write.pred[order(write.pred$ct_value),]
write.pred.long <- tidyr::gather(write.pred, outcome, measurement, both:`missed by n`)
ggplot(write.pred.long, aes(ct_value, measurement)) + geom_line() + facet_wrap(~outcome)
aic <- round(model$AIC)