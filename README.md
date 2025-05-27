# SamovaR <a href=""><img src="data/img/logos/logo_stable.png" align="right" width="150" ></a> 
### Artificial metagenome generation and automatic benchmarking


[![R package](https://github.com/ctlab/samovar/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ctlab/samovar/actions/workflows/R-CMD-check.yaml)
[![python package](https://github.com/ctlab/samovar/actions/workflows/python-package.yml/badge.svg)](https://github.com/ctlab/samovar/actions/workflows/python-package.yaml)

There is a fundamental problem in modern ***metagenomics***: there are huge differences between methodological approaches that strongly influence the results, while remaining outside the attention of researchers. 

The use of golden practice and open code, while allowing data to be analyzed reproducibly, locks scientists into a single, far from perfect approach, with its own bias.

Therefore, we propose an approach that utilizes de novo generation of the artificial metagenomes - `SamovaR`.

## Installation

### Quick Installation

<b><font color="red">Warning:</font></b> beta

Use installation script:

```bash
git clone https://github.com/ctlab/samovar
cd samovar
chmod +x install.sh
./install.sh
```

*Attention: the script automatically detects custom R library paths from `.Renviron` (R_LIBS) or `.Rprofile` (libPaths())*

### Manual Installation

Install R package:

```r
# If you have a custom library path, specify it:
.libPaths('/path/to/your/R/library')
devtools::install_github("https://github.com/ctlab/samovar/")
```

Install python package:

```bash
git clone https://github.com/ctlab/samovar
cd samovar
pip install -e .
```

## Examples

To check installation or dive into what SamovaR do, run

```bash
cd samovar
bash workflow/pipeline.sh
```

## Components

- **R** package `samova.R` for the artificial abundance table generation
- **Python** + **bash** pipeline for the automated benchmarking

## Project Structure

```mermaid
%%{init: {'theme': 'base', 'themeVariables': { 'fontSize': '16px', 'fontFamily': 'arial', 'primaryColor': '#fff', 'primaryTextColor': '#000', 'primaryBorderColor': '#000', 'lineColor': '#000', 'secondaryColor': '#fff', 'tertiaryColor': '#fff'}}}%%
graph LR
    A[SamovaR] --> G1[Abundance table generation]
    G1 --> B[R Package]
    A --> G2[Automated re-profiling]
    G2 --> C[snakemake + Python Pipeline]
    G1 --> G[Shiny App]

    B --> B1[R/]
    B --> B2[man/]
    B --> B3[vignettes/]

    C --> C1[workflow/]
    C --> C2[src/]

    G --> H[shiny/]
```

## Usage
### Cross-validation and re-profiling

To run worklow, specify your config files and run:
```bash
cd samovar
bash workflow/pipeline.sh
```

```mermaid
%%{init: {'theme': 'base', 'themeVariables': { 'fontSize': '16px', 'fontFamily': 'arial', 'primaryColor': '#fff', 'primaryTextColor': '#000', 'primaryBorderColor': '#000', 'lineColor': '#000', 'secondaryColor': '#fff', 'tertiaryColor': '#fff'}}}%%
graph TD
    subgraph Input
        subgraph Metagenomes
            A1[FastQ files]
            A2([InSilicoSeq config])
        end
        A3([Annotation config])
        A4([Generation config])
    end

    subgraph Processing
        Metagenomes --> C[Initial annotation]
        A3 --> C
        C --> E
        A4 --> E[Metagenome generation]
        E --> F[Re-annotation]
        A4 --> F
    end

    subgraph Results
        C --> C1[Cross-validation]
        F --> G1[Annotators scores]
        F --> ML
        subgraph Re-profiling
            C --> R
            ML --> R[Corrected results]
        end
    end

    style Input fill:#90ee9020,stroke:#333,stroke-width:2px
    style Metagenomes fill:#b2ee9020,stroke:#333,stroke-width:2px
    style Processing fill:#ee90bf20,stroke:#333,stroke-width:2px
    style Results fill:#90d8ee20,stroke:#333,stroke-width:2px
    style Re-profiling fill:#90a4ee20,stroke:#333,stroke-width:2px
```

### Artificial metagenome reneration
Basic usage described in <a href="./vignettes">**vignettes**</a> and <a href="https://github.com/ctlab/samovar/wiki">**wiki**</a>

You can also try the generator with <a href="https://dsmutin.shinyapps.io/samovaR/">**web** shiny app</a>


#### R algorithm summary
<img src="data/img/additional/algo.png">


#### R generation

<a href="https://github.com/ctlab/samovar/blob/main/samovaR.pdf">See description</a> or <a href="vignettes/samovar-basic.Rmd">source</a> a vignette

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
<a src="https://github.com/ctlab/samovar/blob/main/samovaR_man.pdf">Documentation</a> for the **R package**

## Algorithm details
<img src="data/img/additional/details.png">

## References
- Chechenina А., Vaulin N., Ivanov A., Ulyantsev V. Development of in-silico models of metagenomic communities with given properties and a pipeline for their generation. Bioinformatics Institute 2022/23 URL: https://elibrary.ru/item.asp?id=60029330


## Dependencies

```mermaid
%%{init: {'theme': 'base', 'themeVariables': { 'fontSize': '16px', 'fontFamily': 'arial', 'primaryColor': '#fff', 'primaryTextColor': '#000', 'primaryBorderColor': '#000', 'lineColor': '#000', 'secondaryColor': '#fff', 'tertiaryColor': '#fff'}}}%%
graph LR
    subgraph "R Package Dependencies"
        subgraph "Main"
            direction LR
            tidyverse
            scclust
            Matrix
            methods
        end
        
        subgraph "Visualization"
            direction LR
            ggplot
            plotly
            tsne
            ggnewscale
        end
        
        subgraph "API"
            direction LR
            httr
            jsonlite
            xml2
        end
    end
    
    subgraph "Automated Benchmarking"
        subgraph "Major"
            direction LR
            samova.R
            snakemake
            InSilicoSeq
        end
        
        subgraph "Python packages"
            direction LR
            numpy
            pandas
            requests
            ete3
        end
    end
    
    linkStyle default stroke:#000
```
