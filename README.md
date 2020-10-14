# Light x Nitrogen greenhouse experiment data repository
Repository for data collected and R scripts created from a light-by-nitrogen greenhouse experiment conducted in the [Smith Ecophys Lab](http://www.smithecophyslab.com/) at Texas Tech University. Two manuscripts are being written from this experiment. The first manuscript (Perkowski et al., in prep) focuses on the effects of light availability, nitrogen availability, and nitrogen acquisition strategy on plant carbon costs to acquire nitrogen. The second manuscript (Waring & Smith, in prep) focuses on the relative role of soil nitrogen availability and environmental conditions on leaf- and whole plant photosynthetic traits.

Currently, this repository contains data relevant to the Perkowski et al. (in prep) manuscript. Data and associated metadata .csv files can be located in the `data_sheets` folder. This data sheet does not include the direct metrics used in the Perkowski et al. (in prep) manuscript; however, calculations for these metrics can be found in the `create_ndemand_metrics.R` located in the `R_files` folder.

The data frame created from `create_ndemand_metrics.R` is directly loaded into a data analysis R file, which can be accessed in the same `R_files` as `ndemand_analyses.R`

NOTE: This repository will be updated with data and metadata from Waring & Smith (in prep) and will be released as a second version on Zenodo.
