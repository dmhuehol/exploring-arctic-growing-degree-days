# exploring-arctic-growing-degree-days
This codebase is primarily intended to support reproducibility of the following manuscript:
 * Hueholt, D.M., E.A. Barnes, J.W. Hurrell, D. Lombardozzi, & A. L. Morrison. "Exploring the Influence of Internal Climate Variability and Forced Change on Arctic Greening" *in prep* for submission to One Earth.

It also supports creation and basic analysis of the dataset of observed growing degree days (GDD) above 50N derived from the Global Historical Climatology Network (GHCN). This purpose may be more broadly useful.

# Basic code description
* `cesm_shell`: shell scripts to process and reprocess CESM output before more complex Python-based analysis
* `ghcn`: code related to importing and manage observed GDDs in the GHCN dataset, as well as dedicated plotting code that only involves GHCN data
* `run_*`: shell scripts used to submit jobs to HPC systems (Colorado State University [CASHEW](https://www.engr.colostate.edu/ets/cashew-cluster/), National Science Foundation National Center for Atmospheric Research [casper](https://ncar-hpc-docs.readthedocs.io/en/latest/compute-systems/casper/))
* `fun_*`: modules to be imported by other scripts and functions
* Some code is not used directly in the manuscript figures or analysis, but a core part of basic data analysis:
  * `wrap_map`: plot a basic map for a slice from a CESM dataset
