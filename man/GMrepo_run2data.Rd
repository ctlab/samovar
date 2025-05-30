% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GMrepo_run2data.R
\name{GMrepo_run2data}
\alias{GMrepo_run2data}
\title{Get data from GMrepo_run object}
\usage{
GMrepo_run2data(
  run,
  number_to_out = F,
  at_level = "species",
  QC_filter = "QCStatus"
)
}
\arguments{
\item{number_to_out}{False by default, maximum number of obtained data}

\item{at_level}{"species" by default. level to obtain classification from GMrepo}

\item{QC_filter}{QCStatus by default. Perform auto QC filtering based on metadata column, or False for no checking.}

\item{runs}{GMrepo_run object got by GMrepo_type2run or created by user with \code{new('GMrepo_run', metadata = data.frame(), run = run_list)}}
}
\description{
Get data from GMrepo_run object
}
\examples{
library(samovaR)
library(tidyverse)

# get data from GMrepo
data_GMrepo <- GMrepo_type2data(mesh_ids = "D006262", number_to_process = 1000)

# equal to:
run_GMrepo <- GMrepo_type2run(mesh_ids = "D006262", number_to_process = 1000)
data_GMrepo <- GMrepo_run2data(run_GMrepo)

# filter runs before obtaining data (OOP updating data!)
run_GMrepo$filter("checking", 1)

# view
data_GMrepo

# access to metadata
data_GMrepo$run

# access to data
data_GMrepo$data

# access to runs
data_GMrepo$run

# access to taxa
data_GMrepo$species
}
