# samovaR <a href=""><img src="additional/Logo_SAMOVAR_pprpl.png" align="right" width="150" ></a> 
### R vignette for generating model metagenomes with specified properties

*Metagenomics* serves as a fundamental approach in the analysis of biological communities. 

The field continuously witnesses the emergence of numerous novel tools, which in turn necessitates the validation of these tools to address the crucial challenge at hand. In light of this, we have devised an artificial data generation tool `SAMOVAR`, aimed at enhancing the development of algorithms and expediting scientific discoveries. This addon implement creation of abundance files 


## Installation

To get the tool clone the git repository:
```bash
git clone https://github.com/dsmutin/samovar.git
```

In your R session, source all functions:
```R
source(PATH/TO/SAMOVAR/scripts/source.R)
```


## Usage

Source R script R_samova.R. To specify properties, change inputs in setup

Avialable options:
`default_path` path to save all output plots and files
sample_amount number of samples
minimal_abundance minimal abundance per species to determine as not the noise
mesh_id GMrepo meshID
number_of_clusters number of clusters to split the data
inner_model model of glm connections within clusters
inter_model model of glm connections between clusters

initial initialized sp
initial_level initialized level
generated_amount number of generated samples
