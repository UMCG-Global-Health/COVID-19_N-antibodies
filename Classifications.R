library(tidyverse)

## Classification of clusters. Cluster numbers refer to original cluster numbers 
## and are therefore not the same as the ordered cluster numbers in the figures of the supplementary files.

df <- df %>%  
  mutate(group_id = case_when(as_clusters_id == 0 ~ 'flat', 
                              as_clusters_id == 1 ~ 'flat' , 
                              as_clusters_id == 2 ~ 'increasing', 
                              as_clusters_id == 3 ~ 'decreasing', 
                              as_clusters_id == 4 ~'increasing', 
                              as_clusters_id == 5 ~ 'de- and increasing',
                              as_clusters_id == 6 ~ 'decreasing', 
                              as_clusters_id == 7 ~ 'increasing', 
                              as_clusters_id == 8 ~ 'flat', 
                              as_clusters_id == 9 ~ 'decreasing', 
                              as_clusters_id == 10 ~ 'decreasing', 
                              as_clusters_id == 11 ~ 'decreasing', 
                              as_clusters_id == 12 ~ 'increasing'
                              ), 
         group_log = case_when(as_clusters_log == 0 ~ 'decreasing', 
                               as_clusters_log == 1 ~ 'de- and increasing', 
                               as_clusters_log == 2 ~ 'increasing', 
                               as_clusters_log == 3 ~ 'flat', 
                               as_clusters_log == 4 ~ 'increasing', 
                               as_clusters_log == 5 ~ 'flat', 
                               as_clusters_log == 6 ~ 'increasing', 
                               as_clusters_log == 7 ~ 'flat', 
                               as_clusters_log == 8 ~ 'flat', 
                               as_clusters_log == 9 ~ 'flat', 
                               as_clusters_log == 10 ~ 'decreasing', 
                               as_clusters_log == 11 ~ 'increasing', 
                               as_clusters_log == 12 ~ 'decreasing'))

# Set final classification 
df <- df %>% 
  mutate(group_final =ifelse(group_id == group_log, group_id, NA),
         group =ifelse(group_id == group_log, 0, 1))

# Classify trajectories with all N-antibodies <=10 
df1 <- df %>%  
  group_by(pseudo_id) %>% 
  filter(all(assay_mabs_n == 10)) %>% 
  mutate(group_final = 'All ≤10', 
         group = 2)

# Classify trajectories with all N-antibodies >=200 
df2 <- df %>%  
  group_by(pseudo_id) %>% 
  filter(all(assay_mabs_n == 200)) %>% 
  mutate(group_final = 'All ≥200', 
         group = 3)

df3 <- df %>%  
  group_by(pseudo_id) %>% 
  filter(!all(assay_mabs_n == 200) & !all(assay_mabs_n == 10)) 

df <- rbind(df1, df2, df3)  
