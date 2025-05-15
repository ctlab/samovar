# samova.R v.0.5 <a href=""><img src="img/logos/logo_stable.png" align="right" width="150" ></a> 
### Artificial metagenome generation and automatic benchmarking

There is a fundamental problem in modern ***metagenomics***: there are huge differences between methodological approaches that strongly influence the results, while remaining outside the attention of researchers. 

The use of golden practice and open code, while allowing data to be analyzed reproducibly, locks scientists into a single, far from perfect approach, with its own bias.

Therefore, we propose an approach that utilizes de novo generation of the artificial metagenomes - `SamovaR`.

## Components

- **R** package `samova.R` for the artificial abundance table generation
- **Python** + **bash** pipeline for the automated benchmarking

## Usage
Basic usage described in <a href="./vignettes">**vignettes**</a> and <a href="https://github.com/ctlab/samovar/wiki">**wiki**</a>

You can also try the generator with <a href="https://dsmutin.shinyapps.io/samovaR/">**web** shiny app</a>

## Installation

To get, run in R:

```r
devtools::install_github("https://github.com/dsmutin/samovar/")
```

## R algorithm summary
<img src="data/img/additional/algo.png">


## Usage

<a href="https://html-preview.github.io/?url=https://github.com/ctlab/samovar/samovar.html">See</a> or <a href="vignettes/samovar-basic.Rmd">source</a>

``` r
library(samovaR)

# download data
teatree <- GMrepo_type2data(number_to_process = 2000)

# filter
tealeaves <- teatree %>%
  teatree_trim(treshhold_species = 3, treshhold_samples = 3, treshhold_amount = 10^(-3))

# normalizing
teabag <- tealeaves %>%
  tealeaves_pack()

# clustering
concotion <- teabag %>%
  teabag_brew(min_cluster_size = 4, max_cluster_size = 6)

# building samovar
samovar <- concotion %>%
  concotion_pour()

# generating new data
new_data <- samovar %>%
  samovar_boil(n = 100)
```

## Documentation
This package can generate new abundance tables based on known abundance tables. It can be done in a few steps: download or import, filtration, normalization, clusterization, generating prediction graphs and generation itself.

Full documentation and tests are avialable in R functions help or in manual and vignettes folder.

## Algorithm details
<img src="data/img/additional/details.png">

## References
- Chechenina А., Vaulin N., Ivanov A., Ulyantsev V. Development of in-silico models of metagenomic communities with given properties and a pipeline for their generation. Bioinformatics Institute 2022/23 URL: https://elibrary.ru/item.asp?id=60029330
- Dai, D. et al. "GMrepo v2: a curated human gut microbiome database with special focus on disease markers and cross-dataset comparison". Nucleic Acids Res (2022). Volume 50, Issue D1, Pages D777–D784.


## Dependencies
### R package
- **main**
  - tidyverse
  - cluster
  - scclust
  - Matrix
  - methods
  - here
- **vizualization** 
  - ggplot
  - plotly
  - tsne
  - ggnewscale
- **API**
  - httr
  - jsonlite
  - xml2


### Automated benchmarking
- **R packages**
- InSilicoSeq
