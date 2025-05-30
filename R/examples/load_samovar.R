library(samovaR)

## build samovar directly from the table
data <- matrix(
  runif(25),
  dimnames = list(
    rownames = paste0("sp_", 1:5),
    colnames = paste0("sample_", 1:5)
  ),
  nrow = 5)

samovar_data <- table2samovar(data)
print(samovar_data)

# or read abundance table
tf <- tempfile()
write.csv(data, tf)

samovar_data <- read_samovar(tf)
print(samovar_data)
