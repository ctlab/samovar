library(samovar)
library(yaml)

# Example config
tf <- tempfile()

write_yaml(
  list(
    metadata_filter = F,
    treshhold_amount = 10^(-5),
    plot_log = F,
    min_cluster_size = 5
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


