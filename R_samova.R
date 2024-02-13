source("scripts/source.R")
teatree <- GMrepo_type2data(test = T)[[1]] 
teabag <- teatree %>% res_normalize(treshhold_amount = 10^(-3))
samovar <- build_samovar(teabag, k_means = 200)
tea <- boil(samovar)