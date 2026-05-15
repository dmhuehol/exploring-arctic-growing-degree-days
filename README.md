# exploring-arctic-growing-degree-days
This codebase is primarily intended to support reproducibility of the following manuscript:
 * [anonymized in peer review] "Exploring the Influence of Internal Climate Variability and Forced Change on Arctic Greening" *submitted*

It also supports derivation of growing degree days (GDD) from observations the Global Historical Climatology Network (GHCN) and basic analysis of this dataset.

## Table of Contents
* [Setting up environment](#setting-up-environment)
* [Replicating figures in Hueholt et al. 2025](#replicating-figures-in-hueholt-et-al-2025)  
* [Workflow to create GDD datasets](#workflow-to-create-gdd-datasets)
  * [GHCN](#ghcn)  
  * [CESM](#cesm)  
* [Code description](#code-description)   
 
## Setting up environment
The `pixi.toml` file documents necessary packages for this repository. Because OpenCV is used in this environment, **`conda` cannot be reliably used as the package manager** and packages should be installed using another method.

## Replicating figures in Hueholt et al. 2025
### Forced crossover (Figure 1)
From `wrap_crossover_map`, run the `crossover_map` function with input data. Requires gridded file of crossover information obtained from `wrap_calc_exceedance_grid` followed by `wrap_calc_crossover_grid` set to threshold of 80%.

### Exceedance timeseries (Figure 2)
Run `wrap_crossover_ts`; input parameters to obtain individual figures documented in file.

### Composites (Figure 3)
Run wrap_composite_crossover_map for composites based on crossover; or wrap_composite_anom_map for composites based on GDD.

### No-analog crossover (Figure 4)
From `wrap_crossover_map`, run the `crossover_map` function with input data. Requires gridded file of crossover information obtained from `wrap_calc_exceedance_grid` followed by `wrap_calc_crossover_grid` set to threshold of 100%.

### Timeseries from GHCN (Figure 5)
Run `plot_ghcn_ts` from the `ghcn` subdirectory.

### Forced crossover and no-analog crossover without shapefile (Supplemental Figure 1)
Same as Figures 1 and 4, but with `paint_shapefile_bool` set to False.

### Sensitivity test (Supplemental Figure 2)
From `wrap_crossover_sensitivity`, run the `test_sensitivity` function to generate sensitivity tests for crossover maps around 80% with perturbations of 1 percentage point, 5 percentage points, and 10 percentage points.

### Difference between biomass burning subsets (Supplemental Figure 3)
Run `wrap_bmb_crossover_map`

### Composites from differing biomass burning subsets (Supplemental Figure 4 and 5)
Same as standard composites, but pointing the data to the "forcing_cmip6" or "forcing_smoothed" directories.

### Histograms comparing LENS2 and GHCN trends (Supplemental Figure 5 and 6)
Run `wrap_trend_hist` from the `ghcn` subdirectory.

### Tables of trends (Supplemental Table 1 and 2)
Reported by `wrap_trend_hist` from the `ghcn` subdirectory.

### Animation of composites (Supplemental Video 1 and 2)
Create individual images with `wrap_multi_composite_map`; animate with `wrap_animate_frames`.

### Animation of station trends with LENS2 context (Supplemental Video 3 and 4)
Use `wrap_animate_ghcn_trend_map` from the `ghcn` subdirectory.

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
