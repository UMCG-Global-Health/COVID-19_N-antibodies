# Create N-antibody positivity based on fixed threshold of 30 ng/ml
library(tidyverse)

# Create variables for before, during and never
df1 <- df %>%  
  group_by(pseudo_id) %>% 
  mutate(thermo_fisher_before = ifelse(first(assay_mabs_n) >= 30, TRUE, FALSE), 
         thermo_fisher_during = ifelse(any(assay_mabs_n >= 30), TRUE, FALSE), 
         thermo_fisher_never = ifelse(all(assay_mabs_n < 30), TRUE, FALSE)) 

# Create variables for before and during
df2 <- df %>%
  group_by(pseudo_id) %>%
  mutate(
    first_mabs = first(assay_mabs_n),
    below_30 = assay_mabs_n < 30,
    index = row_number()
  ) %>%
  arrange(pseudo_id, index) %>%
  mutate(
    any_below_30 = cumsum(below_30) > 0,
    mabs_after_below_30 = ifelse(any_below_30, assay_mabs_n, NA)
  ) %>%
  summarise(
    thermo_fisher_before_during = first(first_mabs >= 30 & any(below_30) & any(mabs_after_below_30 >= 30)))

# Classify threshold based groups
df <- df1 %>% 
  left_join(df2) %>% 
  mutate(thermo_fisher_positivity = case_when(thermo_fisher_before_during == TRUE  ~ 'infected before + during',
                                            (thermo_fisher_before == TRUE)  ~ 'before', 
                                            thermo_fisher_before == FALSE & thermo_fisher_during == TRUE ~ 'infected', 
                                            thermo_fisher_never == TRUE ~ 'never infected' ))
