# exploring-arctic-growing-degree-days
This codebase is primarily intended to support reproducibility of the following manuscript:
 * Hueholt, D.M., E.A. Barnes, J.W. Hurrell, D. Lombardozzi, & A. L. Morrison. "Exploring the Influence of Internal Climate Variability and Forced Change on Arctic Greening" *in prep*

It also supports derivation of growing degree days (GDD) from observations the Global Historical Climatology Network (GHCN) and basic analysis of this dataset.

## Table of Contents
* Setting up environment(#setting-up-environment)
* [Replicating figures in Hueholt et al. 2025](#replicating-figures-in-hueholt-et-al-2025)  
* [Workflow to create GDD datasets](#workflow-to-create-gdd-datasets)
  * [GHCN](#ghcn)  
  * [CESM](#cesm)  
* [Code description](#code-description)   
 
## Setting up environment
The `pixi.toml` file documents necessary packages for this repository. Because OpenCV is used in this environment, **`conda` cannot be reliably used as the package manager** and packages should be installed using another method.

## Replicating figures in Hueholt et al. 2025
To be added.

## Workflow to create GDD datasets
### GHCN
1. Download GHCN dataset [from the National Oceanic and Atmospheric Administration National Centers for Environmental Information (NOAA NCEI)](https://www.ncei.noaa.gov/products/land-based-station/global-historical-climatology-network-daily) or [Amazon Web Services](https://registry.opendata.aws/noaa-ghcn/)
2. Generate a filtered station list with `filter_ghcn` to obtain a subset of interest
3. Create dataset with GDDs and related information using `calc_gdd_info_ghcn`.
4. Make span/coverage dataframe with `check_spancover_ghcn`
5. The GDD dataset and span/coverage dataframe can be used together for plotting and analysis.

### CESM
1. Obtain daily surface temperature output from desired CESM model simulation (see datasheet for links)
2. Create dataset with `calc_gdd`

## Code description
* `cesm_shell`: shell scripts to process and reprocess CESM output before more complex Python-based analysis
* `ghcn`: code related to importing and manage observed GDDs in the GHCN dataset, as well as dedicated plotting code that only involves GHCN data
* `run_*`: shell scripts used to submit jobs to HPC systems (Colorado State University [CASHEW](https://www.engr.colostate.edu/ets/cashew-cluster/), National Science Foundation National Center for Atmospheric Research [casper](https://ncar-hpc-docs.readthedocs.io/en/latest/compute-systems/casper/))
* `fun_*`: modules to be imported by other scripts and functions
* Some code is not used directly in the manuscript figures or analysis, but a core part of basic data analysis:
  * `wrap_map`: plot a basic map for a slice from a CESM dataset
