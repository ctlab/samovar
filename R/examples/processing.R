# download and prepare data
samovar <- GMrepo_type2data(number_to_process = 2000) %>%
  teatree_trim(treshhold_species = 3, treshhold_samples = 3, treshhold_amount = 10^(-3)) %>%
  tealeaves_pack(normalization_function = function(x) log10(x+1)) %>%
  teabag_brew(min_cluster_size = 4, max_cluster_size = 6) %>%
  concotion_pour()

# generate
new_data <- samovar %>%
  samovar_boil(n = 100)


