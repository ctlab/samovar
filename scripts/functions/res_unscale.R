res_unscale <- function(z, data, data_scaled) {
  resF <- data.frame(x = as.numeric(unlist(data_scaled)),
                     y = as.numeric(unlist(data)))
  resF <- resF[order(resF$x),]
  
  if(z <= 0) {
    return(0)
  } else {
    wz <- which(resF$x > z)[1]
    res_sc <- resF[c(wz-1, wz),]
    
    rs <- res_sc$y[1] -
      ((res_sc$x[2] - z) / (res_sc$x[2] - res_sc$x[1])) *
      (res_sc$y[1] - res_sc$y[2])
    
    return(rs)
  }
}