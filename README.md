# samova.R <a href=""><img src="https://github.com/dsmutin/samovar/blob/main/additional/samovar_new_logo.png" align="right" width="150" ></a> 
### R package for generating model metagenomes with specified properties

*Metagenomics* serves as a fundamental approach in the analysis of biological communities. 

The field continuously witnesses the emergence of numerous novel tools, which in turn necessitates the validation of these tools to address the crucial challenge at hand. In light of this, we have devised an artificial data generation tool `SAMOVAR`, aimed at enhancing the development of algorithms and expediting scientific discoveries. This addon implement creation of abundance files 


## Installation

To get the tool clone the git repository:

``` r
devtools::install_github("https://github.com/dsmutin/samovar/tree/beta")
```

## Usage

<a href="https://html-preview.github.io/?url=https://github.com/ctlab/samovar/blob/beta/samovar.html">See</a> or <a href="vignettes/samovar-basic.Rmd">source</a>

``` r
library(samovaR)

# download data
teatree <- GMrepo_type2data(number_to_process = 2000)

# filter
tealeaves <- teatree %>%
  teatree_trim(treshhold_species = 3, treshhold_samples = 3, treshhold_amount = 10^(-3))

# normalizing
## if you build teatree by your own, rescaling stage when building via teatree$rescale() or assigning teatree$min_value and teatree$max_value is required
## good approximation to normal distribution is required for glm generating methods
teabag <- tealeaves %>%
  tealeaves_pack(normalization_function = function(x) log10(x+1))

# clustering
concotion <- teabag %>%
  teabag_brew(min_cluster_size = 4, max_cluster_size = 6)
# remember: if you want to refilter, it is better to re-do welding stage to avoid crashes in future!

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

## References
- Chechenina А., Vaulin N., Ivanov A., Ulyantsev V. Development of in-silico models of metagenomic communities with given properties and a pipeline for their generation. Bioinformatics Institute 2022/23 URL: https://elibrary.ru/item.asp?id=60029330
- Dai, D. et al. "GMrepo v2: a curated human gut microbiome database with special focus on disease markers and cross-dataset comparison". Nucleic Acids Res (2022). Volume 50, Issue D1, Pages D777–D784.


## Required R packages
- Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686
- Taiyun Wei and Viliam Simko (2021). R package 'corrplot': Visualization of a Correlation Matrix (Version 0.92). Available from https://github.com/taiyun/corrplot
- Simon Garnier, Noam Ross, Robert Rudis, Antônio P. Camargo, Marco Sciaini, and Cédric Scherer (2023). viridis(Lite) - Colorblind-Friendly Color Maps for R. viridis package version 0.6.4.
- C. Sievert. Interactive Web-Based Data Visualization with R, plotly, and shiny. Chapman and Hall/CRC, Florida, 2020.
- de Vries A, Ripley BD (2022). _ggdendro: Create Dendrograms and Tree Diagrams Using 'ggplot2'_. R package version 0.1.23, <https://CRAN.R-project.org/package=ggdendro>.
- Wickham H (2023). _httr: Tools for Working with URLs and HTTP_. R package version 1.4.6, <https://CRAN.R-project.org/package=httr>.
- Ooms J (2014). “The jsonlite Package: A Practical and Consistent Mapping Between JSON Data and R Objects.” _arXiv:1403.2805 [stat.CO]_. <https://arxiv.org/abs/1403.2805>.
- Wickham H, Hester J, Ooms J (2023). _xml2: Parse XML_. R package version 1.3.4, <https://CRAN.R-project.org/package=xml2>.
- Donaldson J (2022). _tsne: T-Distributed Stochastic Neighbor Embedding for R (t-SNE)_. R package version 0.1-3.1, <https://CRAN.R-pro
