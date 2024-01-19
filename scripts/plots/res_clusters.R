res_cluster_dendro <- function(
    data_scaled,
    species = T,
    k_means = 1
){
  if (!species) data_scaled <- data_scaled %>% t
  
  ggcl <- data_scaled %>% 
    dist(method = "euclidian") %>% 
    hclust (method = "ward.D2") %>%
    ggdendrogram(leaf_labels = F, segments = T) +
    theme(axis.text = element_text(size = 0.1))
    
  ggplotly(ggcl)
}

res_cluster_PCA2D <- function(
    data_scaled,
    species = T,
    k_means = 1
){
  if (!species) data_scaled <- data_scaled %>% t
  
  samples = colnames(data_scaled)
  clust = data_scaled %>% t %>%
    dist(method = "euclidian") %>% 
    hclust (method = "ward.D2") 

  clusters <- cutree(clust, k = k_means)
  PCA <- data_scaled %>% t %>% prcomp()
  
  gg <- ggplot(mapping = aes (PCA[["x"]][,1], PCA[["x"]][,2], text = samples)) +
    geom_point(aes(colour = clusters), show.legend = F) +
    theme_minimal() +
    xlab(paste0("1 comp, sdev = ", round(PCA[["sdev"]][1], digits = 2))) + 
    ylab(paste0("2 comp, sdev = ", round(PCA[["sdev"]][2], digits = 2)))
  
  ggplotly(gg, tooltip = "text")
  
}
