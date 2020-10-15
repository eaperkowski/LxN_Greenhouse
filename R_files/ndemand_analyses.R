###########################################################
# Load libraries
###########################################################
library(lme4)
library(car)
library(emmeans)

###########################################################
# Load data frame produced by "create_ndemand_metrics.R"
###########################################################
source("https://raw.githubusercontent.com/eaperkowski/LxN_Greenhouse/main/R_files/create_ndemand_metrics.R")

###########################################################
# Create categorical fixed effects
###########################################################
df$n.ppm <- as.factor(df$n.ppm)
df$shade.cover <- as.factor(df$shade.cover)
df$spp <- as.factor(df$spp)

###########################################################
# Mixed effect model structure
###########################################################
#   Fixed  effects:
#     - spp           - 2 levels (cotton, soybean)
#     - shade.cover   - 4 levels (0%, 30%, 50%, 80% shadecover)
#     - n.ppm         - 4 levels (0, 70, 210, 630 ppm N added
#                       twice per week)
#
#   Random effects:
#     - block         - three blocks
#
#   Response variables:
#     - n.acq         - whole plant nitrogen mass
#     - root.carbon   - total root carbon mass
#     - n.cost        - root carbon mass / whole plant nitrogen mass
#     - nod.wgt       - root nodule weight

###########################################################
# Carbon cost to acquire nitrogen (g C g-1 N)
###########################################################
ncost <- lmer(log(n.cost) ~ spp * shade.cover * n.ppm + (1 | block), 
              data = df)

# Check normality assumptions
plot(ncost)
qqnorm(residuals(ncost))
qqline(residuals(ncost))

# Model results
summary(ncost)
Anova(ncost)

# Pairwise comparisons
emmeans(ncost, pairwise ~ shade.cover | n.ppm | spp)
emmeans(ncost, pairwise ~ n.ppm | shade.cover | spp)
emmeans(n.acq, pairwise ~ shade.cover | spp) # NGS: n.acq not yet created
emmeans(n.acq, pairwise ~ n.ppm | spp) # NGS: n.acq not yet created
emmeans(ncost, pairwise ~ spp | n.ppm)

###########################################################
# Whole plant nitrogen mass (g N)
###########################################################
n.acq <- lmer(log(n.acq) ~ spp * shade.cover * n.ppm + (1 | block), 
              data = df)

# Check normality assumptions
plot(n.acq)
qqnorm(residuals(n.acq))
qqline(residuals(n.acq))

# Model results
summary(n.acq)
Anova(n.acq)

# Pairwise comparisons
emmeans(n.acq, pairwise ~ n.ppm | shade.cover | spp)
emmeans(n.acq, pairwise ~ shade.cover | n.ppm | spp)
emmeans(n.acq, pairwise ~ shade.cover | spp)
emmeans(n.acq, pairwise ~ n.ppm | spp)

###########################################################
# Total root carbon (g C)
###########################################################
root.carbon <- lmer(sqrt(root.carbon.mass) ~ spp * shade.cover * n.ppm +
                    (1 | block), data = df)

# Check normality assumptions
plot(root.carbon)
qqnorm(residuals(root.carbon))
qqline(residuals(root.carbon))

# Model results
summary(root.carbon)
Anova(root.carbon)

# Pairwise comparisons
emmeans(root.carbon, pairwise ~ n.ppm | shade.cover | spp)
emmeans(root.carbon, pairwise ~ shade.cover | n.ppm | spp)
emmeans(root.carbon, pairwise ~ shade.cover | spp)
emmeans(root.carbon, pairwise ~ n.ppm | spp)

###########################################################
# Root nodule weight
###########################################################
## NOTE: Mixed effect model does not contain species term
## because cotton is not capable of forming root nodules.
## Thus, species comparisons are excluded from this model
soy.df <- subset(df, spp == "Soybean")
soy.df$nod.wt[soy.df$nod.wt == 0] <- NA

nod.wgt <- lmer(log(nod.wt) ~ shade.cover * n.ppm + (1 | block), 
                data = soy.df)

# Check ormality assumptions
plot(nod.wgt)
qqnorm(residuals(nod.wgt))
qqline(residuals(nod.wgt))

# Model results
summary(nod.wgt)
Anova(nod.wgt)

# Pairwise comparisons
emmeans(nod.wgt, pairwise ~ shade.cover | n.ppm)
emmeans(nod.wgt, pairwise ~ shade.cover)
emmeans(nod.wgt, pairwise ~ n.ppm)
emmeans(nod.wgt, pairwise ~ n.ppm | shade.cover)

