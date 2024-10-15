library(dplyr)

selection <- df %>%  
  group_by(pseudo_id) %>% 
  summarize(N= n()) %>% 
  filter(N>3)

# Exclude all participants with less than 4 measurements
df <- df %>%
  filter(pseudo_id %in% selection$pseudo_id) 

df_ab_less_ten <- df %>%   
  group_by(pseudo_id) %>%  
  filter(all(assay_mabs_n ==10)) %>% 
  select(assay_mabs_n)

# Exclude all participants with all antibodies <=10
df <- df %>% 
  filter(!(pseudo_id %in% df_ab_less_ten$pseudo_id))

df_ab_higher_twohundred <- df %>%    
  group_by(pseudo_id) %>%  
  filter(all(assay_mabs_n ==200)) %>% 
  select(assay_mabs_n)

# Exclude all participants with all antibodies =>200
df <- df %>% 
  filter(!(pseudo_id %in% df_ab_higher_twohundred$pseudo_id))
