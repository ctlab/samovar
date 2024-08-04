cname <- function(x) {
  x %>% 
    str_remove_all("\\..*") %>% 
    str_replace_all(".*Phi.*", "Phix") %>% 
    str_replace_all("Esch.*", "Ecoli") %>% 
    str_replace_all("Enterobact.*", "Ecoli") %>% 
    str_replace_all("Homo.*", "Hsap") %>% 
    str_replace_all("Sach.*", "Scer") %>% 
    str_replace_all(".*uncl.*", "unclassified") %>% 
    str_replace_all("root.*", "root") %>% 
    str_replace_all("other.*", "other")
}

# K2 ----
read_k2 <- function(path, sep = "\t") {
  
  res <- tibble()
  for (i in dir(path, pattern = "out", full.names = T)) {
    
    cat(basename(i),"\t")
    
    if (file.size(i) > 0) {
      tmp <- read.table(i, sep = sep, fill = T, header = F) %>% 
        mutate(sample = i %>%
                 str_remove_all(".*\\/"))
      res <- rbind(res, tmp)
    }
  }
  return(res)
}

# K2 to abund

read_kraken <- function(path) {
  res <- read_k2(path) %>% 
    mutate(V3 = cname(V3)) %>%  
    count(V3, sample) %>% 
    pivot_wider(names_from = V3, id_cols = sample, values_from = n) %>% 
    sapply(function(x) {
      x[is.na(x)] = 0 
      return(x)} ) %>% 
    as.data.frame() %>% 
    column_to_rownames("sample") %>% t %>%  
    apply(c(1,2), as.numeric)
}

read_true <- function(path) {
  res <- read_k2(path) %>% 
    mutate(V2 = cname(V2)) %>%  
    count(V2, sample) %>% 
    pivot_wider(names_from = V2, id_cols = sample, values_from = n) %>% 
    sapply(function(x) {
      x[is.na(x)] = 0 
      return(x)} ) %>% 
    as.data.frame() %>% 
    column_to_rownames("sample") %>% t %>% 
    apply(c(1,2), as.numeric)
}

# kaiju ----

read_kaiju <- function(path) {
  res <- read_k2(path) %>% 
    mutate(V4 = ifelse(V4 == "", "unclassified", cname(V4))) %>%
    count(V4, sample) %>% 
    pivot_wider(names_from = V4, id_cols = sample, values_from = n) %>% 
    sapply(function(x) {
      x[is.na(x)] = 0 
      return(x)} ) %>% 
    as.data.frame() %>%
    column_to_rownames("sample") %>% t %>% 
    apply(c(1,2), as.numeric)
}

# mp4

read_metaphlan <- function (path) {
  
}

samovar2iss <- function (path, 
                         type, res_true = F, 
                         default_path = "benchmarking") {
  
  mesh_id <- type
  
  gsave <- function(gg, postf, height = 1200) {
    ggsave(paste0(default_path,"/",mesh_id, "_", postf, ".png"), gg,
           width = 2000, height = height, units = "px")
  }
  
  
  
  if(type == "k2") res <- read_kraken(path)
  if(type == "kaiju") res<- read_kaiju(path)
  if(type == "mp4") res <- read_metaphlan(path)
  
  if (type == "true") {
    res <- read_true(path) 
    return(res) }
  
  # viz ----
  ch <- c("Ecoli", "Scer", "Phix", "Hsap", "other", "unclassified", "root")
  
  # network
  res_count <- res %>%
    apply(1, mean) %>% 
    as.data.frame() %>% 
    rownames_to_column("sp") %>% 
    mutate("sp_show" = ifelse(sp %in% ch, sp, NA))
  
  #UMAP
  gg <- umap::umap(res)
  
  gg2 <- gg$layout %>% 
    as.data.frame() %>%
    cbind(res_count) %>% 
    ggplot() +
    geom_point(aes(V1, V2, size = ., color = .)) +
    ggrepel::geom_text_repel(aes(V1, V2, label = sp_show)) +
    scale_size_continuous(guide = NULL) +
    theme_minimal() +
    scale_color_gradient(low = "lightblue", high = "darkblue")
  
  gsave(gg2, "umap")
  
  #stats
  
  if (!isFALSE(res_true)) {
    
    #find most related groups
    res_true <- res_true[,colnames(res)]
    
    # DISTANCES_INIT
    distmat <- tibble()
    cormat <- tibble()
    for (sp in unique(rownames(res_true)) ) {
      distances <- res %>% 
        apply(1, function(x) dist(rbind(x, res_true[sp,]))) %>% 
        as.data.frame() %>% 
        rownames_to_column("classified") %>% 
        mutate("true" = sp)
      distmat <- rbind(distmat, distances)
      
      
      distances <- res %>% 
        apply(1, function(x) cor(x, res_true[sp,])) %>% 
        as.data.frame() %>% 
        rownames_to_column("classified") %>% 
        mutate("true" = sp)
      cormat <- rbind(cormat, distances)
    }
    
    gg <- htmp(distmat, 1,3,2, scale = minmaxscale) +
      ggtitle("before correction") +
      scale_color_gradient(name = "dist", low = "lightblue", high = "darkblue")
    gsave(gg, "dist_before")
    
    gg <- htmp(cormat, 1,3,2, scale = F) +
      ggtitle("before correction") +
      scale_fill_gradient2(name = "cor")
    gsave(gg, "cor_before")
    
    #calculate F1
    f1_init <- f1_score(res, res_true)
    
    cat("\nF1-score initial:", f1_init)
    
    
    # calculate R2 ----
    stat_table <- tibble(stat = "f1-score", type = "before", val = f1_init) 
    
    for (sp in rownames(res_true)) {
      if(sp %in% rownames(res)) {
        sp1 <- res_true[sp,]
        sp2 <- res[sp,]
        
        reg <- glm(sp1~sp2) 
        R2 <-  with(summary(reg), 1 - deviance/null.deviance)
        stat_table <- rbind(stat_table,
                            data.frame(stat = paste0(sp, " R^2"),type= "before", val =  R2))
      } else {
        stat_table <- rbind(stat_table,
                            data.frame(stat = paste0(sp, " R^2"),type= "before", val =  0))
      }
    }
    
    
    #re-assign taxa
    newres <- res
    
    for (sp in rownames(newres)) {
      if (!(sp %in% rownames(res_true))) {
        tmp_res <- newres
        
        tmp_cormat <- cormat[cormat$classified == sp,]
        new_sp <- tmp_cormat$true[which.max(tmp_cormat$.)[1]]
        
        if (!(new_sp %in% rownames(tmp_res))) {
          tmp_res <- rbind(tmp_res, 0)
          rownames(tmp_res)[nrow(tmp_res)] <- new_sp
        }
        
        
        tmp_res[new_sp,] <- tmp_res[sp,]
        tmp_res <- tmp_res[-which(rownames(tmp_res) == sp),]
        
        if(f1_score(tmp_res, res_true) > f1_init) newres <- tmp_res
      }
    }
    
    res <- newres
    # DISTANCES_AFTER
    distmat <- tibble()
    cormat <- tibble()
    for (sp in unique(rownames(res_true)) ) {
      distances <- res %>% 
        apply(1, function(x) dist(rbind(x, res_true[sp,]))) %>% 
        as.data.frame() %>% 
        rownames_to_column("classified") %>% 
        mutate("true" = sp)
      distmat <- rbind(distmat, distances)
      
      
      distances <- res %>% 
        apply(1, function(x) cor(x, res_true[sp,])) %>% 
        as.data.frame() %>% 
        rownames_to_column("classified") %>% 
        mutate("true" = sp)
      cormat <- rbind(cormat, distances)
    }
    
    gg <- htmp(distmat, 1,3,2, scale = minmaxscale) +
      scale_color_gradient(name = "dist", low = "lightblue", high = "darkblue") +
      ggtitle("after correction")
    gsave(gg, "dist_after")
    
    gg <- htmp(cormat, 1,3,2, scale = F) +
      ggtitle("after correction") +
      scale_fill_gradient2(name = "cor")
    gsave(gg, "cor_after")
    
    # f1 score after
    
    f1_after <- f1_score(newres, res_true)
    
    cat("\nF1-score after taxa correction:", f1_after)
    
    # calculate R2 ----
    stat_table <- rbind(stat_table, 
                        tibble(stat = "f1-score", type = "after", val = f1_after))
    
    for (sp in rownames(res_true)) {
      if(sp %in% rownames(res)) {
        sp1 <- res_true[sp,]
        sp2 <- res[sp,]
        
        reg <- glm(sp1~sp2) 
        R2 <-  with(summary(reg), 1 - deviance/null.deviance)
        stat_table <- rbind(stat_table,
                            data.frame(stat = paste0(sp, " R^2"),type= "after", val =  R2))
      } else {
        stat_table <- rbind(stat_table,
                            data.frame(stat = paste0(sp, " R^2"),type= "after", val =  0))
      }
    }
    
    #viz stats
    gg <- stat_table %>% 
      mutate(stat = fct_inorder(stat)) %>%
      ggplot(aes(x = type)) +
      geom_col(aes(y = val, fill = val)) + 
      geom_text(aes(y = val+0.2, label = round(val, 2))) +     
      facet_wrap(~stat) +
      theme_minimal() +
      scale_fill_gradient("", low = "lightblue", high = "darkblue") +
      theme(strip.text.x.top = ggtext::element_markdown()) +
      xlab("") + ylab("")
    
    gsave(gg, "stat")
  }
  
  # generate ----
  
  res %>% 
    apply(2, function(x) x / sum(x)) %>% 
    table2iss(mesh_id = mesh_id, 
              output = mesh_id,
              default_path = default_path,
              number_of_clusters = 1
              #normalization_function = function(x) log10(x+1)
    )
  
  return(res)
}