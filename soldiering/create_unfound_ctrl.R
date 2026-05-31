library(tidyverse)

blatt.full <- read_csv("consequences_analysis_full.csv")

blatt.full$found

blatt.full %>% 
  dplyr::filter(found == 0) %>% 
  dplyr::summarise(
  n = n(),
  num_treat = sum(abd == 1),
  num_ctrl = sum(abd == 0)
)

blatt.full.unfound.ctrl <- blatt.full %>% 
  dplyr::filter(found == 0 & abd == 0)

write_csv(blatt.full.unfound.ctrl, "unfound_ctrl.csv")
