library(samovaR)

data <- read_annotation_dir("data/test_annotations/")
gglist <- viz_annotation(data)

samovar_list <- annotation2samovar(data)
