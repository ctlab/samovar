# samova.R: benchmarking <a href=""><img src="additional/logo_benchmarking.png" align="right" width="150" ></a> 
### Vignette for automatic benchmarking of taxonomy annotation tools

Diffrent metagenomics tools produce different results. If most of them can produce some metrics for evaluation of their results, different WGS annotators can not. Most evaluations are based on indirected metrics: unclassified an "noizy" taxa abundance, while overall classification quality should be based on accuracy of taxa assignment. 

Follow general steps or see vignette from example

## Installation

For the benchmarking, you need to install: samovaR for artificial abundance tables generation, InSilicoSeq to generate raw .fastq metagenomics data based on fasterq-dump obtained genomes:

```bash
wget https://github.com/ctlab/samovar/edit/benchmarking
```

In your R session, source all functions:
```R
source(PATH/TO/SAMOVAR/scripts/source.R)
```

Obtain InSilicoSeq and fasterq-dump as follows:

```bash
conda install bioconda::insilicoseq
conda install fasterq-dump
```

Or use possible choices from product websites

## Usage

### Generate or select data

For usage in annotator selection, choose data from one populations group, otherwise artificial populations may be chimeric

For benchmarking, generate data with some artificial species relations

### Annotation

Select different annotation algorithms (i.e. different programms and/or different databases) and make generalized bash functions to use tham on the whole directory and obtain one file for each sample with at least 2 fields: 
- sequence name
- assigned taxa (ncbi taxID prefered)

Run all annotators on selected data

### Import

Import results to R. If you combine them by sequence name field, you than can produce some raw cross-validation tests to see, which taxa are classified differently via different algorythms.

After importing and comparison, summarize data by assigned taxa and sample name fields, and than pivot to produce abundance matrix with `samples x species` dimensions

### Generate new data

Iteratively for each species generate new compositions:

![](https://github.com/ctlab/samovar/blob/main/additional/validation_plots/D006262_composition_generated.png)

And use it for iss generation

### Perform cross-validation



### Evaluate raw scores

### Merge annotations

### References
- Chechenina А., Vaulin N., Ivanov A., Ulyantsev V. Development of in-silico models of metagenomic communities with given properties and a pipeline for their generation. Bioinformatics Institute 2022/23 URL: https://elibrary.ru/item.asp?id=60029330
- Dai, D. et al. "GMrepo v2: a curated human gut microbiome database with special focus on disease markers and cross-dataset comparison". Nucleic Acids Res (2022). Volume 50, Issue D1, Pages D777–D784.


### Required R packages
- Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686
- Taiyun Wei and Viliam Simko (2021). R package 'corrplot': Visualization of a Correlation Matrix (Version 0.92). Available from https://github.com/taiyun/corrplot
- Simon Garnier, Noam Ross, Robert Rudis, Antônio P. Camargo, Marco Sciaini, and Cédric Scherer (2023). viridis(Lite) - Colorblind-Friendly Color Maps for R. viridis package version 0.6.4.
- C. Sievert. Interactive Web-Based Data Visualization with R, plotly, and shiny. Chapman and Hall/CRC, Florida, 2020.
- de Vries A, Ripley BD (2022). _ggdendro: Create Dendrograms and Tree Diagrams Using 'ggplot2'_. R package version 0.1.23, <https://CRAN.R-project.org/package=ggdendro>.
- Wickham H (2023). _httr: Tools for Working with URLs and HTTP_. R package version 1.4.6, <https://CRAN.R-project.org/package=httr>.
- Ooms J (2014). “The jsonlite Package: A Practical and Consistent Mapping Between JSON Data and R Objects.” _arXiv:1403.2805 [stat.CO]_. <https://arxiv.org/abs/1403.2805>.
- Wickham H, Hester J, Ooms J (2023). _xml2: Parse XML_. R package version 1.3.4, <https://CRAN.R-project.org/package=xml2>.
- Donaldson J (2022). _tsne: T-Distributed Stochastic Neighbor Embedding for R (t-SNE)_. R package version 0.1-3.1, <https://CRAN.R-project.org/package=tsne>.
- Maechler, M., Rousseeuw, P., Struyf, A., Hubert, M., Hornik, K.(2022).  cluster: Cluster Analysis Basics and Extensions. R package version 2.1.4.
- Bates D, Maechler M, Jagan M (2023). _Matrix: Sparse and Dense Matrix Classes and Methods_. R package version 1.6-1.1, <https://CRAN.R-project.org/package=Matrix>.
- Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, Xie Y, Allen J, McPherson J, Dipert A, Borges B (2023). _shiny: Web Application Framework for R_. R package version 1.7.5.1, <https://CRAN.R-project.org/package=shiny>.
