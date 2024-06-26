# Setup ----
default_path <- "~/meshid" #path to save all output
sample_amount <- 500 #number of samples
minimal_abundance <- 10^(-4) #minimal abundance per species to determine as not the noise
mesh_id = "D006262" #only implemented
number_of_clusters = 20 #number of clusters to split the data
inner_model <- "gaussian" #model of glm connections within clusters
inter_model <- "gaussian" #model of glm connections between clusters

initial <- "Bifidobacterium bifidum" #initialized sp
initial_level <- 0.002 # initialized level
generated_amount <- 10 #number of generated samples

# Add libraries ----
cat("Initiallized\n")
liblist <- c("tidyverse",
             "corrplot",
             "viridis",
             "plotly",
             "ggtree",
             "ape",
             "httr",
             "jsonlite",
             "xml2",
             "tsne",
             "cluster",
             "Matrix",
             "progress",
             "ggnewscale",
             "htmlwidgets")

for (i in liblist) {
  if (!require(i, quietly = T, character.only = T))
    install.packages(i)
  library(i, character.only = T)
  cat("--- ", i, " loaded  ---\n")
}

rm(i, liblist)

# Configure misc functions ----
trsh <- minimal_abundance

progress_function <- function (iters) {
  if (length(iters) != 1) iters = length(iters)
  
  progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                   total = iters,
                   complete = "=",
                   incomplete = "-",
                   current = ">",
                   clear = FALSE,
                   width = 100)
}

gsave <- function(gg, postf, height = 1200) {
  ggsave(paste0(default_path,"/",mesh_id, "_", postf, ".png"), gg,
         width = 2000, height = height, units = "px")
}

psave <- function(pp, postf){
  htmlwidgets::saveWidget(as_widget(pp), paste0(default_path,"/",mesh_id, "_", postf, ".html"))
}

which.max.coord <- function(df, pred) {
  df <- as.data.frame(df)
  df[,-pred] <- 0
  df[pred,] <- 0
  
  df_row <- apply(df, 1, max)
  df_row <- which.max(df_row)[1]
  
  df_col <- apply(df, 2, max)
  df_col <- which.max(df_col)[1]
  
  return(c(df_col, df_row))
}

prob_cl <- function(a, b){
  and <- ((a > 0) & (b > 0)) %>% sum
  or  <- ((a > 0) | (b > 0)) %>% sum
  return(and/or)
}

composition <- function(data_generated) {
  data_generated[is.na(data_generated)] <- 0
  data_generated[data_generated < trsh] <- 0
  
  split_samples <- hclust(dist(data_generated %>% t))
  
  data_generated <- data_generated[,split_samples[["order"]]]
  
  gg <- data_generated %>% 
    mutate("mean" = apply(.,1,mean)) %>% 
    dplyr::arrange(.$mean) %>%
    rownames_to_column("sp") %>% 
    pivot_longer(cols = -1) %>% 
    subset(value > 0) %>% 
    mutate(name = fct_inorder(name)) %>% 
    mutate(sp = fct_inorder(.$sp)) %>% 
    ggplot(aes(name, value, fill = sp, text = sp)) +
    geom_col(position = "stack", show.legend = F, color = "white", linewidth = .5, alpha = .7) +
    xlab("sample") +
    scale_fill_viridis_d() +
    theme_void()
  
  pp <- ggplotly(gg) %>% style(showlegend = FALSE)
  
}

# Read data ----
cat("Read data\n")
## read data from file or get it from web
run_info <- readLines("https://raw.githubusercontent.com/dsmutin/samovar/main/scripts/test/health_runIDs.txt")[1:sample_amount] 

## set config to ignore certificate
httr::set_config(httr::config(ssl_verifypeer = FALSE))
options(RCurlOptions = list(ssl_verifypeer = FALSE))
options(rsconnect.check.certificate = FALSE)

pb <- progress_function(length(run_info))
res_sp <- tibble(taxa = character()) 

## Read abundance data frame from run list ----
for (RUN in run_info) {
  pb$tick()
  
  # get run info from GMrepo
  query <- httr::POST("https://gmrepo.humangut.info/api/getFullTaxonomicProfileByRunID",
                      body = list("run_id" = RUN), encode = "json"
  ) %>%
    httr::content() %>%
    xml2::xml_text() %>%
    jsonlite::fromJSON()
  
  if (!is.null(query[["species"]])) {
    df_sp <- query[["species"]][, -2] %>%
      summarise(N = sum(relative_abundance), .by = scientific_name)
    colnames(df_sp) <- c("taxa", RUN)
    res_sp <- full_join(res_sp, df_sp, by = "taxa")
  }
}
rm(df_sp)
rm(query)

res_sp <- res_sp %>% column_to_rownames("taxa")

# Vizualize composition ----
cat("Composition vizualization")
pp <- res_sp %>% 
  composition

psave(pp, "composition_initial")


# Data processing ----
cat("Prepare\n")
#remove all taxa with abundance < 2
sumNA <- !is.na(res_sp[,-c(1:2)])
sumNA <- apply(sumNA, 1, sum)
sumNA <- which(sumNA > 2)

res_sp <- res_sp[sumNA,]

#make % from abundance
res_sp <- res_sp / 100
res_sp[is.na(res_sp)] <- 0
#make tables for further analysis
species <- res_sp %>% rownames

# we further need tables with median and sd numbers for each group

res_stat <- data.frame("sp" = rownames(res_sp),
                       "mn" = apply(res_sp, 1, mean, na.rm = T),
                       "sd" = apply(res_sp, 1, sd, na.rm = T))

trsh = minimal_abundance
resF <- unlist(res_sp) %>% as.numeric()
distr <- (res_sp > trsh) %>% apply(1, sum) %>% rep(each = length(res_sp[1,])) / length(res_sp[1,])
distr <- distr[resF > trsh]
resF <- resF[resF > trsh]

ord = resF %>% order(decreasing = T)

resF <- resF[ord]
distr <- distr[ord]
N <- 1:length(resF)

split_n <- 0.25
gg <- ggplot() +
  geom_jitter(mapping = aes(y = N[distr<=split_n], resF[distr<=split_n]), 
              color = "lightyellow",
              alpha = 0.5, size = 0.3, 
              width = max(resF)/25, height = max(N)/25) +
  geom_jitter(mapping = aes(y = N[distr>split_n], resF[distr>split_n], 
                            color = distr[distr>split_n]),
              alpha = 0.5, size = 0.3, 
              width = max(resF)/25, height = max(N)/25) +
  geom_smooth(aes(y = N, resF), method = "glm", 
              method.args = list(family = "Gamma"),
              color = "red",
              linewidth = 0.5,
              linetype = 2,
              formula = 'y~x') +
  scale_color_gradient2("consistency", low = "lightyellow", mid = "lightgreen", 
                        high = "darkblue",midpoint = 0.4) +
  xlab("mean abundance value") + ylab("index") + ggtitle("Means distribution") +
  theme_minimal()

gsave(gg, "distribution")

# model re-normalization
sum0 <- apply(res_sp > trsh, 1, sum)
res_sp <- res_sp[sum0 > 2,]
species <- species[which(sum0 > 2)]

res_stat <- res_stat[which(sum0 > 2),]

res_scale <- function(x) {
  x <- log10(1/(1-log10(x+1)))*(1-2*log10(x+1))
  return(x)
}

res_sp_scale <- res_scale(res_sp)

resF <- unlist(res_sp) %>% as.numeric()
resF <- resF[resF > trsh]
ord = resF %>% order(decreasing = T)
resF <- resF[ord]
resF2 <- resF %>% res_scale

gg <- data.frame(res = resF, type = "before", N = 1:length(resF)) %>% 
  rbind(data.frame(res = resF2, type = "after", N = 1:length(resF))) %>% 
  ggplot(mapping = aes(sample = res)) +
  geom_qq()+
  geom_qq_line() +
  facet_wrap(~type, scales = "free") +
  #geom_smooth(method = "glm", formula = 'y~x') +
  theme_minimal() +
  ggtitle("QQ-plot")

gsave(gg, "QQ_plot")

# Get scaled stats ----
min_res_sp_scale <- min(res_sp_scale, na.rm = T)
res_sp_scale <- res_sp_scale - min_res_sp_scale
res_sp_scale[is.na(res_sp_scale)] <- 0
res_sp_scale <- as.matrix(res_sp_scale)

res_stat_scale <- data.frame("sp" = species,
                             "mn" = apply(res_sp, 1, mean, na.rm = T),
                             "sd" = apply(res_sp, 1, sd, na.rm = T))

# Generate inverse function ----
res_unscale <- function(z) {
  resF <- data.frame(x = as.numeric(unlist(res_sp_scale)),
                     y = as.numeric(unlist(res_sp)))
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

# Clusterization ----
## HClust ----
cat("Clustering\n")
clust <- res_sp_scale %>% 
  dist(method = "euclidian") %>% 
  hclust (method = "ward.D2")

KM = number_of_clusters
kmeans <- cutree(clust, k = KM)

## Dendrogramm ----
d <- data.frame(label=clust$labels, 
                cluster = kmeans, 
                amount = res_stat$mn, 
                presence = res_sp %>% apply(1, function(x) sum(x > trsh)))
trs <- full_join(ape::as.phylo(clust), d, by = "label")
trs <- list("Cluster" = trs, "Mean abundance" = trs, "Presence" = trs)
class(trs) <- 'treedataList'

gg <- ggtree(trs, layout = "circ") +
  facet_wrap(~(.id %>% fct_rev)) +
  geom_tippoint(data=td_filter(.id == "Cluster"), mapping = aes(color = cluster)) +
  scale_color_viridis_c("",direction = -1) +
  ggnewscale::new_scale_colour()  + 
  geom_tippoint(data=td_filter(.id == "Mean abundance"), mapping = aes(color = amount)) +
  scale_colour_gradient(name = "", low='lightyellow', high='red') +
  ggnewscale::new_scale_colour()  + 
  geom_tippoint(data=td_filter(.id == "Presence"), mapping = aes(color = presence)) +
  scale_colour_gradient(name = "", low='lightyellow', high='red') +
  theme(legend.position = "bottom")
#pp <- ggplotly(gg) 

gsave(gg, "cluster")

## PCA, etc ----
fit <- kmeans(res_sp_scale, KM)

PCA <- prcomp(res_sp_scale %>% apply(1, scale) %>% t )

gg <- ggplot(mapping = aes (PCA[["x"]][,1], PCA[["x"]][,2], text = species)) +
  geom_point(aes(colour = kmeans), show.legend = F) +
  theme_minimal() +
  xlab(paste0("1 comp, sdev = ", round(PCA[["sdev"]][1], digits = 2))) + 
  ylab(paste0("2 comp, sdev = ", round(PCA[["sdev"]][2], digits = 2))) +
  scale_color_viridis_c(direction = -1)

pp <- ggplotly(gg, tooltip = "text")
psave(pp, "PCA")

cat("Perform tSNE\n")

tSNE2D <- res_sp %>% apply(2, scale) %>% 
  tsne(k = 2, initial_dims = KM) %>% data.frame

gg <- ggplot(tSNE2D, text = species) + 
  geom_point(aes(x=X1, y=X2, color = kmeans), 
             show.legend = T, alpha = 0.5) +
  xlab("") + ylab("") +
  scale_color_continuous("", type = "viridis") +
  theme_minimal()

pp <- ggplotly(gg, tooltip = "text")
psave(pp, "tSNE_2D")

tSNE3D <- res_sp %>% apply(2, scale) %>% 
  tsne(k = 3) %>% data.frame

pp <- plot_ly(data = tSNE3D,
        x =  ~X1, y = ~X2, z = ~X3,
        color = kmeans) %>%
  add_markers(size = 8) %>%
  layout( 
    xaxis = list(
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'), 
    yaxis = list(
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'))

psave(pp, "tSNE_3D")

# Build samovar ----
## Inner-clusters ----
cat("Number of clusters:", KM, "\n")
cat("Samovar initialized")
# get R square tables and GLMs for different clusters
Rsq_cl <- list()
Rsq_pr <- list()


# subsetting tables and make 
for (k in 1:KM) {
  clust_num <- (kmeans == k)
  cat("\nCluster", k, "\n")
  if (sum(clust_num) > 1) {
    df_k <- res_sp_scale[clust_num,]
    ldfk <- length(df_k[,1])
    
    df_r2 <- matrix(ncol = ldfk, nrow = ldfk)
    df_pr <- matrix(ncol = ldfk, nrow = ldfk)
    
    pb <- progress_function(ldfk)
    for(i in 1:ldfk){
      pb$tick()
      j <- i
      for(j in i:ldfk) {
        if(i == j) {
          df_r2[i,j] <- 0
        } else {
          
          dfkj <- df_k[j,]
          dfki <- df_k[i,]
          dfNA <- !((dfki>0) & (dfkj>0))
          dfBOTH <- ((dfki>0) | (dfkj>0)) %>% sum
          dfkj <- dfkj[!dfNA]
          dfki <- dfki[!dfNA]
          
          #get the prob = W(a&b)/W(a|b)
          probij <- (sum(!dfNA) / dfBOTH)
          if(sum(!dfNA) > 1) {
            glmij <- glm(dfki ~ dfkj, family = inner_model)
            df_r2[i,j] <- with(summary(glmij), 1 - deviance/null.deviance)
          } else {
            df_r2[i,j] <- 0
          }
          df_pr[i,j] <- probij
        }
        j <- j + 1
      }
    }
    df_r2 <- Matrix::forceSymmetric(df_r2)
    df_pr <- Matrix::forceSymmetric(df_pr)
  } else {
    df_r2 <- 0
    df_pr <- 0
  }
  
  Rsq_cl[[k]] <- df_r2
  Rsq_pr[[k]] <- df_pr
  
  #cat("%", k, "%")
}

cat("Generate inter-cluster connections\n")
## Inter-clusters ----
# get R square tables and GLMs for between clusters
Rsq_il <- matrix(ncol = KM, nrow = KM)
Glm_il <- list()

pb <- progress_function(KM)
for (i in 1:KM){
  pb$tick()
  df_k1 <- res_sp_scale[kmeans == i,]
  if(!is.null(ncol(df_k1))) df_k1 <- apply(df_k1, 2, mean) %>% unlist %>% as.numeric
  j <- i
  
  while (j <= KM) {
    
    if(i == j) {
      Rsq_il[i,j] <- 0
      j <- j + 1
    } else {
      df_k2 <- res_sp_scale[kmeans == j,]
      if(!is.null(ncol(df_k2))) df_k2 <- apply(df_k2, 2, mean) %>% unlist %>% as.numeric
      
      glmij <- glm(df_k1 ~ df_k2, family = inter_model)
      Rsq_il[i,j] <- with(summary(glmij), 1 - deviance/null.deviance)
      Glm_il[[as.character(i)]][[as.character(j)]] <- glmij
      j <- j + 1
    }
  }
}

Rsq_il <- Matrix::forceSymmetric(Rsq_il) %>% as.matrix

# Generation ----
cat("\nSamovar built, generation started...")
data_generated <- data.frame(sp = res_sp %>% rownames)

for (i in 1:generated_amount) {
  seed_rand <- runif(1, 0, 10000) %>% round
  cat("\nSample", i, "of", generated_amount, "- seed -",seed_rand, "- prediction\n")
  set.seed(seed_rand)
  
  init <- initial
  init_level <- initial_level
  
  res <- data.frame()
  pb <- progress_function(KM)
  cl_todo <- 1:KM
  cl_done <- c()
  
  #check is there some cluster else to do and start with init cluster
  cl <- kmeans[which(init == species)]
  
  #rescale init
  init_level <- res_scale(init_level) - min_res_sp_scale
  
  while (length(cl_todo)>0) {
    pb$tick()
    #make predictions for relevant cluster
    #for now: easiest prediction just based on one group. 
    #in future, find the best path must be based on ML
    
    #first: make cluster tables for R2 and means 
    sp_cl <- species[cl == kmeans]
    df_k <- res_sp_scale[cl == kmeans,]
    df_r_cl <- Rsq_cl[[cl]] %>% as.matrix
    df_r_pr <- Rsq_pr[[cl]] %>% as.matrix
    df_r_pr[is.na(df_r_pr)] <- 0
    
    res_cl <- data.frame("sp" = sp_cl,
                         "res_scale" = NA)
    
    #prediction of initial level 
    if(!is.na(init)){
      #for represented cluster
      pred <- data.frame(from = which(init == sp_cl), to = which(init == sp_cl))
      res_cl[which(sp_cl == init),2] <- init_level
      
      #second: start predictions based on biggest R2
      while(sum(is.na(res_cl[,2]))) {
        
        #prediction graph
        pr <- which.max.coord(df_r_pr, pred$to)
        pred <- rbind(pred, pr)
        
        if(sum(pr == c(1,1)) == 2) {
          res_cl[is.na(res_cl)] <- 0
        }
        
        #prediction of presence
        if(runif(1) <= df_r_pr[pr[1], pr[2]]) {
          
          df <- df_r_cl
          df[-pred$from,] <- 0
          best <- which.max(df[,pr[2]])
          
          #fit GLM
          dfi <- df_k[best,]
          dfj <- df_k[pr[2],]
          
          dfij <- ((dfi > 0) & (dfj > 0)) %>% which
          
          if(length(dfij) < 2){
            #edit in future; for now if there are less than 2 connected observations for 2 groups level of 2 group <- 0
            res_cl[pr[2],2] <- 0
            
          } else {
            dfij <- data.frame("i" = as.numeric(dfi[dfij]), 
                               "j" = as.numeric(dfj[dfij]))
            glmij <-glm(j~i, data = dfij, family = inner_model)
            res_cl[pr[2],2] <- predict.glm(glmij, newdata = data.frame(i = res_cl[pr[1],2]))
          }
          
        } else {
          #the best way is bayesian correction if prob <- 0, so we should re-celculate:
          # prob <- w(!a&b)/w(!a|b)
          
          #but for now, lets just remove that probabilities
          df_r_pr[pr[2],] <- 0
          df_r_pr[,pr[2]] <- 0
          
          res_cl[pr[2],2] <- 0
        }
        
      }
      
    } else {
      #for not represented cluster
      res_cl$res_scale <- 0
    }
    
    #add results
    res <- rbind(res, res_cl)
    cl_todo <- cl_todo[-cl]
    cl_done <- c(cl_done, cl)
    
    #select new cluster for generation and predict initial levels
    #it is also better to use bayesian correction
    if(length(cl_todo) > 0) {
      cl_new <- which.max.coord(Rsq_il, cl_done)
      
      if((cl_new[1] == 1) & (cl_new[2] == 1)) break
      
      df_k1 <- df_k
      if(!is.null(ncol(df_k1))) df_k1 <- apply(df_k1, 2, mean)
      
      df_k <- res_sp_scale[cl_new[2] == kmeans,]
      df_k2 <- df_k
      if(!is.null(ncol(df_k2))) df_k2 <- apply(df_k2, 2, mean)
      
      df_i12 <- (df_k1 > 0) & (df_k2 > 0)
      df_o12 <- (df_k1 > 0) | (df_k2 > 0)
      
      #random generation of cluster presence based on Bayesian model
      if(runif(1) < (sum(df_i12)/sum(df_o12))) {
        
        #mean of cluster level
        dfk12 <- data.frame("i" = as.numeric(df_k2[df_i12]),
                            "j" = as.numeric(df_k1[df_i12]))
        glmij <- glm(i~j, data = dfk12, family = inter_model)
        cl_level <- predict.glm(glmij, newdata = data.frame(j = mean(res_cl$res_scale)), type = "response")
        
        #select init
        sp_cl <- species[cl_new[2] == kmeans]
        
        if(length(sp_cl) == 1) {
          res <- rbind(res, c(sp_cl, cl_level))
        } else {
          probs <- apply(df_k, 1, prob_cl, df_k1)
          probs <- probs / sum(probs)
          init <- sample(sp_cl, size = 1, prob = probs)
          df_ki <- df_k[which(init == sp_cl),] %>% as.numeric
          df_k2i <- data.frame("i" = as.numeric(df_ki),
                               "j" = as.numeric(df_k2))
          glmij <- glm(i~j, data = dfk12, family = inter_model)
          init_level <- predict.glm(glmij, newdata = data.frame(j = mean(res_cl$res_scale)), type = "response")
        }
        
      } else {
        sp_cl <- species[cl_new[2] == kmeans]
        init <- NA
      }
      
      cl <- cl_new[2]
    }
  }
  
  #rescale
  res$res_scale <- as.numeric(res$res_scale)
  
  cat("rescaling...")
  res$res_scale <- sapply(res$res_scale, res_unscale)

  res <- res[!duplicated(res[,1]),]
  res$res_scale[res$res_scale < 0] <- 0
  if(sum(res$res_scale) > 1 ) res$res_scale = res$res_scale / sum (res$res_scale)
  res$res_scale[res$sp == "Unknown"] <- 1-sum(res$res_scale[res$res_scale != "Unknown"])
  
  colnames(res)[2] <- i
  data_generated <- full_join(data_generated, res, by = "sp")
}

# Viz results ----

data_generated <- data_generated %>% 
  column_to_rownames("sp") %>% 
  mutate_all(as.numeric)

data_generated[data_generated < trsh] <- 0

pp <- data_generated %>% 
  composition

psave(pp, "composition_generated")

# PCA
merged_df <- cbind(res_sp, data_generated)
PCA <- prcomp(merged_df %>% apply(1, scale))

cols <- c(
  rep("initial", ncol(res_sp)),
  rep("generated", ncol(data_generated))
)

gg <- ggplot(mapping = aes (PCA[["x"]][,1], PCA[["x"]][,2], text = colnames(merged_df))) +
  geom_point(aes(colour = cols), show.legend = F) +
  theme_minimal() +
  xlab(paste0("1 comp, sdev = ", round(PCA[["sdev"]][1], digits = 2))) + 
  ylab(paste0("2 comp, sdev = ", round(PCA[["sdev"]][2], digits = 2))) +
  scale_color_viridis_d(direction = -1)

pp <- ggplotly(gg, tooltip = "text")
psave(pp, "generatedPCA")

# New species tSNE
keep_df <- apply(data_generated, 1, sum)>0
tSNE2D <- data_generated %>% 
  subset(apply(., 1, sum)>0) %>% 
  apply(2, scale) %>% 
  tsne(k = 2, initial_dims = KM) %>% 
  data.frame

gg <- ggplot(tSNE2D, text = species[keep_df]) + 
  geom_point(aes(x=X1, y=X2, color = kmeans[keep_df]), 
             show.legend = T, alpha = 0.5) +
  xlab("") + ylab("") +
  scale_color_continuous("", type = "viridis") +
  theme_minimal()

pp <- ggplotly(gg, tooltip = "text")
psave(pp, "tSNE_2D_generated")

# Save ----
cat("Write data...")

write.csv(data_generated, file = paste0(default_path,"/",mesh_id, "_generated.csv"))

cat("All done")