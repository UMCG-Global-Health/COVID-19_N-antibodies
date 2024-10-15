library(tidyverse)

# Difference between two measurements per participant
df <- df %>% 
  group_by(pseudo_id) %>% 
  mutate(antibody.difference = assay_mabs_n - lag(assay_mabs_n))

# Find max difference and days of the maximum difference between two measurements and the day before
max <- df %>%  
  group_by(pseudo_id) %>% 
  summarise(max = max(antibody.difference, na.rm=TRUE), 
            max_day = days[which.max(antibody.difference)], 
            day_before_max= days[which.max(antibody.difference)-1])

# Create midpoint and (hypothetical) N-antibody infection date
df <- df %>% 
  left_join(max, by = 'pseudo_id') %>% 
  mutate(midpoint = (max_day + day_before_max)/2, 
         n_ab_pos_date = ifelse(n_ab_class %in% c('increasing', 'de- and increasing') & 
                                  swab_pos %in% c('Positive swab during', 'Positive swab before and during'), midpoint - 14, NA))
