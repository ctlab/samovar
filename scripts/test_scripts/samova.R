#get the environment
source("scripts/source.R")

#get and load runs from GMrepo

#You can use 2 steps download: first, download runIDs from GMrepo, and second - download samples based on GMrepo runs, using GMrepo_type2run(...) %>% GMrepo_run2data(...), or do it on one stage using GMrepo_type2data(...). For test, use option test = True
#runs <- GMrepo_type2run(mesh_ids = "D006262", number_to_process = 10000) #idk why but this results to runs with no species list
runs <- read.table("scripts/test/health_runIDs.txt") %>% unlist # - or simple use downloaded runs list
runs <- GMrepo_type2run(test = T) %>% unlist %>% GMrepo_run2data(runs, number_to_process = 250)
res_list <- GMrepo_type2data(test = T)
rm(runs)
#here are some differences with normal data in metadata but fix later

#on this stage you can filter results using metadata 
data <- res_list[[1]] 

#plot and visualize data scaling
data %>% res_trim %>% plot_data_with_treshhold
data %>% res_trim %>% plot_data_n2amount(normalisation_function = log10)
data %>% res_trim %>% plot_data_n2amount(normalisation_function = function(x) log10(1/(1-log10(x)))*(1-2*log10(x)))

#apply
data <- data %>% res_trim (treshhold_amount = 10^(-3),)
teabag <- data %>% 
  #res_normalize()
  res_normalize(normalisation_function = function(x) log10(1/(1-log10(x)))*(1-2*log10(x)))

teabag$data %>% plot_data_with_treshhold

#analyze data
teabag$data %>% res_cluster_dendro
teabag$data %>% res_cluster_PCA2D

#make structure to analyze
samovar <- build_samovar(teabag, k_means = 100)

#create new data!
new_reads <- boil(data, data_scaled, samovar, prepared = T)



