###########################################################
## Libraries
###########################################################
library(dplyr)

###########################################################
## Load light.nitrogen.datasheet
###########################################################
df <- read.csv("https://raw.githubusercontent.com/eaperkowski/LxN_Greenhouse/main/data_sheets/lightnitrogen.datasheet.csv",
               stringsAsFactors = FALSE,
               na.strings = c("", "NA"))

###########################################################
## Create whole-plant nitrogen mass (n.acq; gN), root 
## carbon mass (root.carbon.mass; gC), and carbon cost to 
## acquire nitrogen (n.cost; gC gN-1) metrics.
###########################################################
df <- df %>%
  mutate(n.acq = ((stem.wt * (n.stem / 100)) + 
                    (leaves.wt * (n.leaf / 100)) + 
                    (roots.wt * (n.root / 100))),
         root.carbon.mass = roots.wt * (c.root / 100),
         n.cost = root.carbon.mass / n.acq,
         nod.root.ratio = nod.wt / roots.wt)