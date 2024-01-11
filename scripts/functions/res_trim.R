res_trim <- function(data, treshhold_amount = 10^(-5), treshhold_presence = 2) {
  
  #make tables for further analysis
  row.names(data) <- data[,2]
  data <- data[,-c(1:2)]
  
  #make % from abundance
  sum_data <- data %>% apply(2, sum, na.rm = T)
  data <- data / sum_data
  
  #remove values below treshhold_amount
  data[data < treshhold_amount] <- NA
  
  #remove all taxa with abundance below treshhold_presence
  sumNA <- is.na(data)
  sumNA <- apply(sumNA, 1, sum)
  sumNA <- which(sumNA > treshhold_presence)
  
  data <- data[sumNA,]
  
  return(data)
}