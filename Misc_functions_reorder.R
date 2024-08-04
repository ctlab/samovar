reord <- function(df, cols, rows, vals, scale_function = scale) {
  tmp_col = colnames(df)
  tmp_row = df[,tmp_col[cols]] %>% unique
    
  df_wide <- df %>% 
    pivot_wider(id_cols = cols, names_from = rows, values_from = vals) %>% 
    column_to_rownames(tmp_col[cols]) 
    
  if(!isFALSE(scale_function)) {
    df_wide <- df_wide %>% 
      apply(2, scale_function) %>% 
      "rownames<-"(tmp_row)
  }
  
  dc1 <- df_wide %>% 
    dist %>% hclust
  
  dc2 <- df_wide %>% 
    t %>% dist %>% hclust
  
  df_wide <- df_wide[dc1$order,dc2$order] %>%
    as.data.frame() %>% 
    rownames_to_column(tmp_col[cols]) %>% 
    pivot_longer(-1, values_to = tmp_col[vals], names_to = tmp_col[rows]) %>% 
    as.data.frame()
  
  df_wide[,tmp_col[rows]] <- fct_inorder(df_wide[,tmp_col[rows]])
  df_wide[,tmp_col[cols]] <- fct_inorder(df_wide[,tmp_col[cols]])
  
  df_wide
}

minmaxscale <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

trunc <- function (x) {
  x %>% 
    str_remove_all(" \\(.*") %>% 
    str_extract("(^[a-zA-Z]* [a-zA-Z]*)|(^[a-zA-Z]*)")
}

htmp <- function(df, ...) {
  gg <- df  %>% 
    reord(...) %>% 
    mutate(classified = trunc(classified) %>% fct_inorder) %>% 
    ggplot(aes(classified, true)) +
    geom_tile(aes(fill = .)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, face = "italic", size = 5)) 
}

dist_to_true <- function(res, res_true, dist_f) {
  distmat <- tibble()
  for (sp in unique(rownames(res_true)) ) {
    distances <- res %>% 
      apply(1, dist_f) %>% 
      as.data.frame() %>% 
      rownames_to_column("classified") %>% 
      mutate("true" = sp)
    distmat <- rbind(distmat, distances)
  }
  return(distmat)
}


gsave <- function(gg, postf, height = 1200) {
  ggsave(paste0(default_path,"/",mesh_id, "_", postf, ".png"), gg,
         width = 2000, height = height, units = "px")
}

merge_with_true <- function(res, res_true) {
  tmp <- res %>%
    as.data.frame() %>% 
    rownames_to_column("sp") %>% 
    pivot_longer(-1) %>% 
    mutate(newal = ifelse(sp %in% ch, sp, "other"))
  
  tmp2 <- res_true %>% 
    as.data.frame() %>% 
    rownames_to_column("sp") %>% 
    pivot_longer(-1, values_to = "true")
  
  tmp <- full_join(tmp, tmp2, by = c("sp", "name"))
  tmp$value[is.na(tmp$value)] <- 0
  tmp$true[is.na(tmp$true)] <- 0
  
  return(tmp)
}

f1_score <- function(res, res_true) {
  tmp <- merge_with_true(res, res_true) %>% 
    mutate(neg = abs(value - true)) %>% 
    mutate(pos = true - neg) 
  tmp$pos[tmp$pos < 0] <- 0
  
  tp <- sum(tmp$pos)
  fpfn <- sum(tmp$neg)
  
  f1 <- tp / (tp + 1/2 * fpfn)
  f1 %>% round(2)
}


