library(samovaR)

data <- read_annotation_dir(system.file("extdata", "test_annotations", package = "samovaR"))
gglist <- viz_annotation(data)

samovar_list <- annotation2samovar(data)
