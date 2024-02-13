res_normalize = function(
    data,
    normalisation_function = function(x) log10(x + 1),
    treshhold_amount = 10^(-5), 
    treshhold_presence = 2,
    drop_unclassified = T
) {
  #trim
  data <- data %>% res_trim(treshhold_amount, treshhold_presence, drop_unclassified)
  
  #normalize
  data <- normalisation_function (data)
  min_data <- min(data, na.rm = T)
  max_data <- max(data, na.rm = T)
  data <- data - min_data
  data <- data / (max_data - min_data)
  data[is.na(data)] <- 0
  
  cat("---  Data rescaled to [0,1]  ---\n")
  
  #get inverse
  inverse = function (f, lower, upper) {
    function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1] %>% 
      unlist %>% as.numeric
  }
  reverse_function = inverse(function (x) normalisation_function(x), 
                             lower = 0, upper = 1)
  
  #reverse_function = function(x) reverse_function(x * max_data + min_data)
  

  return(list("data" = data, 
              "parameters" = list (
                "trim" = list(
                  "treshhold_amount" = treshhold_amount, 
                  "treshhold_presence" = treshhold_presence, 
                  "drop_unclassified" = drop_unclassified
                ),
                "normalize" = list(
                  "normalisation_function" = normalisation_function, 
                  "reverse_function" = reverse_function, 
                  "min_data" = min_data, 
                  "max_data" = max_data)
                )
              )
         )
}
