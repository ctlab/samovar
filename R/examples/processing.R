library(samovaR)
library(tidyverse)

# download and prepare data
samovar_raw <- GMrepo_type2data(number_to_process = 2000)

samovar_data <- samovar_raw %>%
  samovar_preprocess()

# Similar to:
## samovar_raw %>%
##   teatree_trim() %>%
##   tealeaves_pack() %>%
##   teabag_brew() %>%
##   concotion_pour()

# generate
new_data <- samovar_data %>%
  samovar_boil(N = 100)
