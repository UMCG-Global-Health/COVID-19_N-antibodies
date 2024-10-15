library(tidyverse)
library(sjPlot)

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

# Create dummies for heterogeneity effect of age and ct values
df <- df %>% mutate(dummy_age_1 = ifelse(age_at_visit > 3, 1, 0),
                    age_at_visit1 = (age_at_visit -3)*dummy_age_1, 
                    dummy_age_2 = ifelse(age_at_visit > 6, 1, 0),
                    age_at_visit2 = (age_at_visit -6)*dummy_age_2, 
                    dummy_ct_1 = ifelse(ct_value >20, 1, 0), 
                    ct_value1 = (ct_value-20) *dummy_ct_1,
                    dummy_ct_2 = ifelse(ct_value >30, 1,0),
                    ct_value2 = (ct_value -30)*dummy_ct_2)


# Transform age and ct values according to knots
df <- df %>%  mutate(age1 = age_at_visit *(!dummy_age_1) + 3*dummy_age_1, 
                     age2 = (age_at_visit-3) *(dummy_age_1 & !dummy_age_2) + 3*dummy_age_2,
                     age3 = (age_at_visit - 6)* dummy_age_2 ,
                     ct1 = ct_value *(!dummy_ct_1) + 20*dummy_ct_1, 
                     ct2 = (ct_value-20) *(dummy_ct_1 & !dummy_ct_2) + 10*dummy_ct_2,
                     ct3 = (ct_value - 30)* dummy_ct_2)

# Logistic regression
model_final <- glm(seropositivity_binary ~ age1 + age2 + age3+
                              sex+ ethnicity_wo + ever_hcw + ever_lthc + vaccination+
                              ct1 +  ct2 + ct3 +  symptoms + type, 
                            data = df, family = 'binomial')


tab_model(model_final, file = 'model.doc')

# Model for heterogeneity value of age and ct value
model_final <- glm(seropositivity_binary ~ age_at_visit + age_at_visit1 + age_at_visit2+
                       ct_value +  ct_value1 + ct_value2  + symptoms +
                       sex+ ethnicity_wo + ever_hcw + ever_lthc + vaccination + type , 
                     data = df, family = 'binomial')

summary(model_final)

#######univariate analyses

model1 <- glm(seropositivity_binary ~ age1 + age2 + age3, 
              data = df, family = 'binomial')

tab_model(model1, file = 'model_age.doc')

model2 <- glm(seropositivity_binary ~ sex, 
              data = df, family = 'binomial')

tab_model(model2, file = 'model_sex.doc')

model3 <- glm(seropositivity_binary ~ ethnicity_wo, 
              data = df, 
              family = 'binomial')

tab_model(model3, file = 'model_eth.doc')


model4 <- glm(seropositivity_binary ~ ever_hcw, 
              data = df, 
              family = 'binomial')

tab_model(model4, file = 'model_hcw.doc')


model5 <- glm(seropositivity_binary ~ ever_lthc, 
              data = df, 
              family = 'binomial')

tab_model(model5, file = 'model_lthc.doc')


model6 <- glm(seropositivity_binary ~ vaccination, 
              data = df, 
              family = 'binomial')

tab_model(model6, file = 'model_vac.doc')

model7 <- glm(seropositivity_binary ~ symptoms, 
              data = df, 
              family = 'binomial')

tab_model(model7, file = 'model_symp.doc')


model8 <- glm(seropositivity_binary ~ ct1 + ct2 + ct3, 
              data = df, 
              family = 'binomial')

tab_model(model8, file = 'model_ct.doc')

model9 <- glm(seropositivity_binary ~ type, 
              data = df, 
              family = 'binomial')

tab_model(model9, file = 'model_type.doc')
