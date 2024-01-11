res_normalize = function(
    data,
    normalisation_function = log10
) {
  data <- normalisation_function (data)
  min_data <- min(data, na.rm = T)
  data <- data - min_data
  data[is.na(data)] <- 0
  return(data)
}
