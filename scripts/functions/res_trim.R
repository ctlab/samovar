res_trim <- function(data, 
                     treshhold_amount = 10^(-5), 
                     treshhold_presence = 2,
                     drop_unclassified = T) {
  
  #make tables for further analysis
  row.names(data) <- data[,2]
  data <- data[,-c(1:2)]
  
  if (drop_unclassified) {
    data <- data %>% 
      subset(rownames(data) %>% str_detect("nclassified", negate = T)) %>%
      subset(rownames(data) %>% str_detect("nknown", negate = T))
  }
  
  #make % from abundance
  sum_data <- data %>% apply(2, sum, na.rm = T)
  data <- data / sum_data
  
  #remove values below treshhold_amount
  data[data < treshhold_amount] <- NA
  
  #remove all taxa with abundance below treshhold_presence
  sumNA <- !is.na(data)
  keep_rows <- apply(sumNA, 1, sum) > treshhold_presence
  keep_cols <- apply(sumNA, 2, sum) > 0
  
  data <- data[keep_rows,keep_cols]
  
  return(data)
}
