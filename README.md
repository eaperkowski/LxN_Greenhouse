# Light x Nitrogen greenhouse experiment data repository
Repository for data collected and R scripts created from a light-by-nitrogen greenhouse experiment conducted in the [Smith Ecophys Lab](http://www.smithecophyslab.com/) at Texas Tech University. 

Two manuscripts are being written from this experiment. The first manuscript (Perkowski et al., 2021; [link to paper](https://academic.oup.com/jxb/article/72/15/5766/6296480?login=true)) focuses on the effects of light availability, nitrogen availability, and nitrogen acquisition strategy on plant carbon costs to acquire nitrogen and was recently published in *Journal of Experimental Botany*. The second manuscript (Waring, Perkowski, and Smith, in prep) focuses on the role of soil nitrogen availability and light availability on leaf and whole plant photosynthetic traits.

Currently, this repository contains data relevant to Perkowski et al. (2021). Data and associated metadata files can be located in the `data_sheets` folder. This data sheet does not include the direct metrics used in Perkowski et al. (2021); however, calculations for these metrics can be found in the `LxN_ncost_create_metrics.R` file located in the `R_files` folder. Units for each metric are coded out in this file.

The data frame created from `LxN_ncost_create_metrics.R` is directly loaded into a data analysis R file, which can be found in the `LxN_ncost_analyses.R` file. This data frame is also directly loaded into a file containing manuscript plot code, which can be found in the `LxN_ncost_plots.R` file. Both the `LxN_ncost_analyses.R` and `LxN_ncost_plots.R` files can be accessed in the `R_files` folder.

# UPDATE: August 21, 2022
Data, metadata, and scripts for data analysis and figuremaking are now pushed for the second manuscript investigating impacts of nitrogen fertilization and light availability on leaf and whole plant physiology. Data are included in the `data_sheets` folder as `LxN_physiology_data.csv` and metadata are included as `LxN_physiology_metadata.csv`. Scripts for data analysis and figuremaking are included in the `R_files` folder as `LxN_phys_analyses_plots.R`.

NOTE: Metadata file does not include chlorophyll or proportion of N in photosynthesis units. Calculations and units can be found in the `R_files/LxN_phys_analyses_plots.R` script.

## Zenodo DOI (current release includes data for JXB paper and pre-submission data for PC&E physiology paper)
[![DOI](https://zenodo.org/badge/304118064.svg)](https://zenodo.org/badge/latestdoi/304118064)

## Published papers that include these data
Perkowski EA, Waring EF, Smith NG. 2021. Root mass carbon costs to acquire nitrogen are determined by nitrogen and light availability in two species with different nitrogen acquisition strategies (A Rogers, Ed.). Journal of Experimental Botany 72: 5766–5776. DOI: [10.1093/jxb/erab253](https://doi.org/10.1093/jxb/erab253)
