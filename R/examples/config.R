library(samovar)
library(yaml)

# Example config
tf <- tempfile()

write_yaml(
  list(
    treshhold_amount = 10^(-5),
    plot_log = F,
    min_cluster_size = 5,
    N = 5,
    N_reads = 100
  ),
  tf
)

config_samovar <- unpack_config(tf)

## build samovar directly from the table
data <- matrix(
  rlnorm(625),
  dimnames = list(
    rownames = paste0("sp_", 1:25),
    colnames = paste0("sample_", 1:25)
  ),
  nrow = 25)

samovar_data <- table2samovar(data)

# Run with config
config_samovar$samovar_data <- samovar_data
samovar <- do.call(samovar_preprocess, config_samovar)
samovar_new <- samovar_boil(samovar, N = config_samovar$N)

new_data <- samovar_new$data * config_samovar$N_reads

heatmap(as.matrix(new_data))
