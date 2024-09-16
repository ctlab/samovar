sample_scale <- function (res_sp, trsh = 10^(-5)) {
  #remove all taxa with abundance < 2
  sumNA <- !is.na(res_sp[,-c(1:2)])
  sumNA <- apply(sumNA, 1, sum)
  sumNA <- which(sumNA > 2)

  res_sp <- res_sp[sumNA,]

  #make tables for further analysis
  species <- res_sp[,2]
  taxa_id <- res_sp[,1]

  res_sp <- res_sp[,-c(1:2)]

  #make % from abundance
  res_sp <- res_sp / 100

  #scale
  res_sp[res_sp < trsh] <- NA
  sum0 <- !is.na(res_sp)
  sum0 <- apply(sum0, 1, sum)
  res_sp <- res_sp[sum0 > 2,]
  species <- species[which(sum0 > 2)]

  res_stat <- res_stat[which(sum0 > 2),]

  #scaling function
  res_scale <- function(x) {
    x <- log10(1/(1-log10(x)))*(1-2*log10(x))
  }

  res_sp_scale <- res_scale(res_sp)
  rownames(res_sp_scale) <- species
}

#to be implemented
inverse <- function(inv0, lower = -100, upper = 100) {
  inv1 <- function (f, lower = -100, upper = 100) {
    function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1] }
  inv2 <- inv1(function (x) log10(x), 0, 100)
  return(inv2)
}
