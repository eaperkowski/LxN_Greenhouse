# Light x Nitrogen greenhouse experiment data repository
Repository for data collected and R scripts created from a light-by-nitrogen greenhouse experiment conducted in the [Smith Ecophys Lab](http://www.smithecophyslab.com/) at Texas Tech University. 

Two manuscripts are being written from this experiment. The first manuscript (Perkowski et al., 2021; [link to paper](https://academic.oup.com/jxb/article/72/15/5766/6296480?login=true)) focuses on the effects of light availability, nitrogen availability, and nitrogen acquisition strategy on plant carbon costs to acquire nitrogen and was recently published in *Journal of Experimental Botany*. The second manuscript (Waring & Smith, in prep) focuses on the relative role of soil nitrogen availability and environmental conditions on leaf and whole plant photosynthetic traits.

Currently, this repository contains data relevant to Perkowski et al. (2021). Data and associated metadata .csv files can be located in the `data_sheets` folder. This data sheet does not include the direct metrics used in Perkowski et al. (2021); however, calculations for these metrics can be found in the `create_ndemand_metrics.R` file located in the `R_files` folder. Units for each metric are coded out in this file.

The data frame created from `create_ndemand_metrics.R` is directly loaded into a data analysis R file, which can be found in the `ndemand_analyses.R` file. This data frame is also directly loaded into a file containing manuscript plot code, which can be found in the `ndemand_plots.R` file. Both the `ndemand_analyses.R` and `ndemand_plots.R` files can be accessed in the `R_files` folder.

NOTE: This repository will be updated with data and metadata from Waring & Smith (in prep) and will be released as a second version on Zenodo.

## Zenodo DOI
[![DOI](https://zenodo.org/badge/304118064.svg)](https://zenodo.org/badge/latestdoi/304118064)

## Published papers that include these data
Perkowski EA, Waring EF, Smith NG. 2021. Root mass carbon costs to acquire nitrogen are determined by nitrogen and light availability in two species with different nitrogen acquisition strategies (A Rogers, Ed.). Journal of Experimental Botany 72: 5766–5776. DOI: [10.1093/jxb/erab253](https://doi.org/10.1093/jxb/erab253)
