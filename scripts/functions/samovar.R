samovar <- function(data, reps = 1, ...) {
  data <- data %>% res_trim (...)
  data_scaled <- data %>% res_normalize(...)
  samovar <- build_samovar(...)
  boil(data, data_scaled, samovar, reps = reps, ...)
}