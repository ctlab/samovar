plot_data_with_treshhold <- function(
    data,
    treshhold_amount = 10^(-4),
    split_n = 0.3,
    normalisation_function = function(x) x
){
  resF <- data.frame(amount = unlist(data))
  resF$amount <- resF$amount %>% normalisation_function
  resF$distr <- (data > treshhold_amount) %>% 
    apply(1, sum, na.rm = T) %>% 
    rep(each = ncol(data)) / ncol(data)
  
  resF$amount[resF$amount < treshhold_amount] <- NA
  resF = resF[!is.na(resF$amount),]
  resF = resF[order(resF$amount), ]
  resF$N = nrow(resF):1
  
  ggplot() +
    geom_jitter(mapping = aes(y = resF$N[resF$distr<=split_n], resF$amount[resF$distr<=split_n]), 
                color = "lightyellow",
                alpha = 0.5, size = 0.3, 
                width = max(resF$amount)/25, height = max(resF$N)/25) +
    geom_jitter(mapping = aes(y = resF$N[resF$distr>split_n], resF$amount[resF$distr>split_n], 
                              color = resF$distr[resF$distr>split_n]),
                alpha = 0.5, size = 0.3, 
                width = max(resF$amount)/25, height = max(resF$N)/25) +
    
    geom_smooth(aes(y = resF$N, resF$amount), method = "glm", 
                method.args = list(family = "Gamma"),
                color = "red",
                linewidth = 0.5,
                linetype = 2,
                formula = 'y~x') +
    scale_color_gradient2("consistency", low = "lightyellow", mid = "lightgreen", 
                          high = "darkblue",midpoint = 0.4) +
    xlab("means") + ylab("index") + ggtitle("Means distribution") +
    theme_minimal()
}

#####
plot_data_n2amount <- function (
    data,
    normalisation_function = function(x) x,
    split_n
) {
  resF <- unlist(data)
  resF <- resF[order(resF, decreasing = F)]
  resF <- resF[resF > 0]
  resF <- resF[!is.na(resF)]
  N <- 1:length(resF)
  
  ggplot(mapping = aes(normalisation_function(resF), N)) +
    geom_point(size = 0.5) +
    geom_smooth(method = "glm", formula = 'y~x') +
    xlab("amount") + ylab("N") + 
    ylim(0, length(resF)) +
    theme_minimal()
}