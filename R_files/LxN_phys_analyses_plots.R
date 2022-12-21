# Analysis script for the light x N greenhouse experiment
# Note that root directory is the R_files folder of the LxN_Greenhouse repository
# Analysis script for the light x N greenhouse experiment

## load packages
library(tidyverse)
library(lme4)
library(emmeans)
library(multcomp)
library(car)
library(gtable)
library(grid)

## load data
data = read.csv('../data_sheets/LxN_phys_data.csv')
head(data)

## define a function for multi-panel plots
multiplot <- function(..., plotlist=NULL, cols) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # Make the panel
  plotCols = cols                          # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
  
  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }
  
}
## leaf N and SLA
marea_lmer <- lmer(log(marea) ~ spp * shade_cont * fertilizer_cont + (1| block), data = data)
plot(resid(marea_lmer) ~ fitted(marea_lmer))
summary(marea_lmer)
Anova(marea_lmer)
emmeans(marea_lmer, ~spp)
test(emtrends(marea_lmer, ~ spp, var = 'shade_cont'))
cld(emtrends(marea_lmer, ~ spp, var = 'shade_cont'))
### marea higher in cotton
### marea decreases with shading

nmass_lmer <- lmer(log(nmass) ~ spp * shade_cont * fertilizer_cont + (1| block), data = data)
plot(resid(nmass_lmer) ~ fitted(nmass_lmer))
summary(nmass_lmer)
Anova(nmass_lmer)
test(emtrends(nmass_lmer, ~ spp, var = 'fertilizer_cont',
              at = list(shade_cont = c(0))))
test(emtrends(nmass_lmer, ~ spp, var = 'fertilizer_cont',
              at = list(shade_cont = c(30))))
test(emtrends(nmass_lmer, ~ spp, var = 'fertilizer_cont',
              at = list(shade_cont = c(50))))
test(emtrends(nmass_lmer, ~ spp, var = 'fertilizer_cont',
              at = list(shade_cont = c(80))))
test(emtrends(nmass_lmer, ~ spp, var = 'shade_cont',
              at = list(fertilizer_cont = c(0))))
test(emtrends(nmass_lmer, ~ spp, var = 'shade_cont',
              at = list(fertilizer_cont = c(70))))
test(emtrends(nmass_lmer, ~ spp, var = 'shade_cont',
              at = list(fertilizer_cont = c(210))))
test(emtrends(nmass_lmer, ~ spp, var = 'shade_cont',
              at = list(fertilizer_cont = c(630))))
### positive fertilizer effect for cotton at all shade levels
### no fertilizer effect for soybean at any shade level

narea_lmer <- lmer(log(narea) ~ spp * shade_cont * fertilizer_cont + (1| block), data = data)
plot(resid(narea_lmer) ~ fitted(narea_lmer))
summary(narea_lmer)
Anova(narea_lmer)
test(emtrends(narea_lmer, ~ spp, var = 'fertilizer_cont'))
cld(emtrends(narea_lmer, ~ spp, var = 'fertilizer_cont'))
test(emtrends(narea_lmer, ~ spp, var = 'fertilizer_cont',
              at = list(shade_cont = c(0))))
test(emtrends(narea_lmer, ~ spp, var = 'fertilizer_cont',
              at = list(shade_cont = c(30))))
test(emtrends(narea_lmer, ~ spp, var = 'fertilizer_cont',
              at = list(shade_cont = c(50))))
test(emtrends(narea_lmer, ~ spp, var = 'fertilizer_cont',
              at = list(shade_cont = c(80))))
test(emtrends(narea_lmer, ~ 1, var = 'fertilizer_cont',
              at = list(shade_cont = c(0))))
test(emtrends(narea_lmer, ~ 1, var = 'fertilizer_cont',
              at = list(shade_cont = c(30))))
test(emtrends(narea_lmer, ~ 1, var = 'fertilizer_cont',
              at = list(shade_cont = c(50))))
test(emtrends(narea_lmer, ~ 1, var = 'fertilizer_cont',
              at = list(shade_cont = c(80))))
### positive effect of fertilizer for both, stronger in cotton
### positive effect of fertilizer diminishes at lowest light in soybean

## photosynthesis analyses
vcmax25_lmer <- lmer(log(vcmax25) ~ spp * shade_cont * fertilizer_cont + (1| block), data = data)
plot(resid(vcmax25_lmer) ~ fitted(vcmax25_lmer))
summary(vcmax25_lmer)
Anova(vcmax25_lmer)
test(emtrends(vcmax25_lmer, ~ 1, var = 'shade_cont'))
test(emtrends(vcmax25_lmer, ~ spp, var = 'fertilizer_cont'))
cld(emtrends(vcmax25_lmer, ~ spp, var = 'fertilizer_cont'))
### negative effect of shading throughout
### positive effect of fertilizer in cotton, but not soybean

jmax25_lmer <- lmer(log(jmax25) ~ spp * shade_cont * fertilizer_cont + (1| block), data = data)
plot(resid(jmax25_lmer) ~ fitted(jmax25_lmer))
summary(jmax25_lmer)
Anova(jmax25_lmer)
test(emtrends(jmax25_lmer, ~ spp, var = 'shade_cont'))
cld(emtrends(jmax25_lmer, ~ spp, var = 'shade_cont'))
test(emtrends(jmax25_lmer, ~ spp, var = 'fertilizer_cont'))
cld(emtrends(jmax25_lmer, ~ spp, var = 'fertilizer_cont'))
### negative effect of shading throughout, stronger in cotton
### positive effect of fertilizer in cotton, but not soybean

jv25_lmer <- lmer((jv25) ~ spp * shade_cont * fertilizer_cont + (1| block), data = data)
plot(resid(jv25_lmer) ~ fitted(jv25_lmer))
summary(jv25_lmer)
Anova(jv25_lmer)
cld(emmeans(jv25_lmer, ~spp))
### higher in soybean

chlorophyll_mmol.m2_lmer <- lmer((chl_mmol.m2) ~ spp * shade_cont * fertilizer_cont + (1| block), data = data)
plot(resid(chlorophyll_mmol.m2_lmer) ~ fitted(chlorophyll_mmol.m2_lmer))
summary(chlorophyll_mmol.m2_lmer)
Anova(chlorophyll_mmol.m2_lmer)
test(emtrends(chlorophyll_mmol.m2_lmer, ~ spp, var = 'shade_cont'))
test(emtrends(chlorophyll_mmol.m2_lmer, ~ spp, var = 'fertilizer_cont'))
cld(emmeans(chlorophyll_mmol.m2_lmer, ~spp))
# cld(emtrends(chlorophyll_mmol.m2_lmer, ~ spp, var = 'fertilizer_cont'))
### reduced with shade in soy, but not cotton
### positive effect of fertilizer in cotton, but not soybean

## per leaf n metrics
# vcmax25_narea_lmer <- lmer(log(vcmax25_narea) ~ spp * shade_cont * fertilizer_cont + (1| block), data = data)
# plot(resid(vcmax25_narea_lmer) ~ fitted(vcmax25_narea_lmer))
# summary(vcmax25_narea_lmer)
# Anova(vcmax25_narea_lmer)
# test(emtrends(vcmax25_narea_lmer, ~ 1, var = 'shade_cont'))
# test(emtrends(vcmax25_narea_lmer, ~ spp, var = 'fertilizer_cont'))
# cld(emtrends(vcmax25_narea_lmer, ~ spp, var = 'fertilizer_cont'))
# test(emtrends(vcmax25_narea_lmer, ~ 1, var = 'fertilizer_cont',
#               at = list(shade_cont = 0)))
# test(emtrends(vcmax25_narea_lmer, ~ 1, var = 'fertilizer_cont',
#               at = list(shade_cont = 30)))
# test(emtrends(vcmax25_narea_lmer, ~ 1, var = 'fertilizer_cont',
#               at = list(shade_cont = 50)))
# test(emtrends(vcmax25_narea_lmer, ~ 1, var = 'fertilizer_cont',
#               at = list(shade_cont = 80)))
# ### negative effect of shading throughout
# ### negative effect of fertilizer, stronger in cotton
# ### negative fertilizer effect reduced with shade
# 
# jmax25_narea_lmer <- lmer(log(jmax25_narea) ~ spp * shade_cont * fertilizer_cont + (1| block), data = data)
# plot(resid(jmax25_narea_lmer) ~ fitted(jmax25_narea_lmer))
# summary(jmax25_narea_lmer)
# Anova(jmax25_narea_lmer)
# test(emtrends(jmax25_narea_lmer, ~ 1, var = 'shade_cont'))
# test(emtrends(jmax25_narea_lmer, ~ spp, var = 'fertilizer_cont'))
# cld(emtrends(jmax25_narea_lmer, ~ spp, var = 'fertilizer_cont'))
# test(emtrends(jmax25_narea_lmer, ~ 1, var = 'fertilizer_cont',
#               at = list(shade_cont = 0)))
# test(emtrends(jmax25_narea_lmer, ~ 1, var = 'fertilizer_cont',
#               at = list(shade_cont = 30)))
# test(emtrends(jmax25_narea_lmer, ~ 1, var = 'fertilizer_cont',
#               at = list(shade_cont = 50)))
# test(emtrends(jmax25_narea_lmer, ~ 1, var = 'fertilizer_cont',
#               at = list(shade_cont = 80)))
# ### negative effect of shading throughout
# ### negative effect of fertilizer, stronger in cotton
# ### negative fertilizer effect reduced with shade
# 
# chlorophyll_mmol.m2_narea_lmer <- lmer(log(chlorophyll_mmol.m2_narea) ~ spp * shade_cont * fertilizer_cont + (1| block), data = data)
# plot(resid(chlorophyll_mmol.m2_narea_lmer) ~ fitted(chlorophyll_mmol.m2_narea_lmer))
# summary(chlorophyll_mmol.m2_narea_lmer)
# Anova(chlorophyll_mmol.m2_narea_lmer)
# test(emtrends(chlorophyll_mmol.m2_narea_lmer, ~ 1, var = 'shade_cont'))
# test(emtrends(chlorophyll_mmol.m2_narea_lmer, ~ 1, var = 'fertilizer_cont'))
# ### negative effect of shading throughout
# ### negative effect of fertilizer throughout

propN_rubisco_lmer <- lmer(propN_rubisco ~ spp * shade_cont * fertilizer_cont + (1| block), data = data)
plot(resid(propN_rubisco_lmer) ~ fitted(propN_rubisco_lmer))
summary(propN_rubisco_lmer)
Anova(propN_rubisco_lmer)
test(emtrends(propN_rubisco_lmer, ~ 1, var = 'shade_cont'))
test(emtrends(propN_rubisco_lmer, ~ 1, var = 'fertilizer_cont'))
emmeans(propN_rubisco_lmer, ~spp)
### recuction with shade
### reduction with fertilizer

propN_bioenergetics_lmer <- lmer(propN_bioenergetics ~ spp * shade_cont * fertilizer_cont + (1| block), data = data)
plot(resid(propN_bioenergetics_lmer) ~ fitted(propN_bioenergetics_lmer))
summary(propN_bioenergetics_lmer)
Anova(propN_bioenergetics_lmer)
test(emtrends(propN_bioenergetics_lmer, ~ 1, var = 'shade_cont'))
test(emtrends(propN_bioenergetics_lmer, ~ 1, var = 'fertilizer_cont'))
emmeans(propN_bioenergetics_lmer, ~spp)
### recuction with shade
### reduction with fertilizer

propN_lightharvesting_lmer <- lmer(propN_lightharvesting ~ spp * shade_cont * fertilizer_cont + (1| block), data = data)
plot(resid(propN_lightharvesting_lmer) ~ fitted(propN_lightharvesting_lmer))
summary(propN_lightharvesting_lmer)
Anova(propN_lightharvesting_lmer)
test(emtrends(propN_lightharvesting_lmer, ~ 1, var = 'shade_cont'))
test(emtrends(propN_lightharvesting_lmer, ~ 1, var = 'fertilizer_cont'))
emmeans(propN_lightharvesting_lmer, ~spp)
### no effect of shade
### reduction with fertilizer

propN_photosynthesis_lmer <- lmer(propN_photosynthesis ~ spp * shade_cont * fertilizer_cont + (1| block), data = data)
plot(resid(propN_photosynthesis_lmer) ~ fitted(propN_photosynthesis_lmer))
summary(propN_photosynthesis_lmer)
Anova(propN_photosynthesis_lmer)
test(emtrends(propN_photosynthesis_lmer, ~ 1, var = 'shade_cont'))
test(emtrends(propN_photosynthesis_lmer, ~ 1, var = 'fertilizer_cont'))
### reduction with shade
### reduction with fertilizer

## whole-plant metrics
tla_lmer <- lmer(log(tla) ~ spp * shade_cont * fertilizer_cont + (1| block), data = data)
plot(resid(tla_lmer) ~ fitted(tla_lmer))
summary(tla_lmer)
Anova(tla_lmer)
test(emtrends(tla_lmer, ~ spp, var = 'shade_cont'))
cld(emtrends(tla_lmer, ~ spp, var = 'shade_cont'))
test(emtrends(tla_lmer, ~ 1, var = 'fertilizer_cont'))
test(emtrends(tla_lmer, ~ 1, var = 'fertilizer_cont',
              at = list(shade_cont = 0)))
test(emtrends(tla_lmer, ~ 1, var = 'fertilizer_cont',
              at = list(shade_cont = 30)))
test(emtrends(tla_lmer, ~ 1, var = 'fertilizer_cont',
              at = list(shade_cont = 50)))
test(emtrends(tla_lmer, ~ 1, var = 'fertilizer_cont',
              at = list(shade_cont = 80)))
### reduced with shade in both, but to a greater degree in cotton
### increase with fertilizer
### ferilizer effect is slightly reduced with shading

biomass_lmer <- lmer(log(biomass) ~ spp * shade_cont * fertilizer_cont + (1| block), data = data)
plot(resid(biomass_lmer) ~ fitted(biomass_lmer))
summary(biomass_lmer)
Anova(biomass_lmer)
test(emtrends(biomass_lmer, ~ spp, var = 'shade_cont'))
cld(emtrends(biomass_lmer, ~ spp, var = 'shade_cont'))
test(emtrends(biomass_lmer, ~ 1, var = 'fertilizer_cont'))
test(emtrends(biomass_lmer, ~ 1, var = 'fertilizer_cont',
              at = list(shade_cont = 0)))
test(emtrends(biomass_lmer, ~ 1, var = 'fertilizer_cont',
              at = list(shade_cont = 30)))
test(emtrends(biomass_lmer, ~ 1, var = 'fertilizer_cont',
              at = list(shade_cont = 50)))
test(emtrends(biomass_lmer, ~ 1, var = 'fertilizer_cont',
              at = list(shade_cont = 80)))
### reduced with shade in both, but to a greater degree in cotton
### increase with fertilizer
### ferilizer effect is slightly reduced with shading

## some takehomes
### fertilizer and sun increase leaf N
### the leaf N increases are more in support of photosynthetic processes for sun than fertilizer
### proportion of leaf N in photosynthesis (all components!) is reduced with shade and fertilizer!
### the whole-plant effects are strong with fertilizer and sun (both positive)

## make figures and tables

### TABLES
marea_anova_table <- Anova(marea_lmer)
nmass_anova_table <- Anova(nmass_lmer)
narea_anova_table <- Anova(narea_lmer)
marea_nmass_narea_anova_table <- cbind(marea_anova_table, nmass_anova_table, narea_anova_table)
marea_nmass_narea_anova_table_sub <- cbind(marea_nmass_narea_anova_table[,2],
                                           marea_nmass_narea_anova_table[, c(1, 3, 4, 6, 7, 9)])
colnames(marea_nmass_narea_anova_table_sub) <- c('df', 'χ2', 'P-value', 'χ2', 'P-value', 'χ2', 'P-value')
rownames(marea_nmass_narea_anova_table_sub) <- c('Species (Sp)', 'Shading (Sh)', 'Fertilizer (F)',
                                                 'Sp x Sh', 'Sp x F', 'Sh x F', 'Sp x Sh x F')
is.num <- sapply(marea_nmass_narea_anova_table_sub, is.numeric)
marea_nmass_narea_anova_table_sub[is.num] <- lapply(marea_nmass_narea_anova_table_sub[is.num], 
                                                    round, 3)
marea_nmass_narea_anova_table_sub[marea_nmass_narea_anova_table_sub<0.001] <- '<0.001'
write.csv(marea_nmass_narea_anova_table_sub, 'tables/marea_nmass_narea_anova_table.csv')

vcmax25_anova_table <- Anova(vcmax25_lmer)
jmax25_anova_table <- Anova(jmax25_lmer)
chlorophyll_mmol.m2_anova_table <- Anova(chlorophyll_mmol.m2_lmer)
vcmax25_jmax25_chlorophyll_mmol.m2_anova_table <- cbind(vcmax25_anova_table, jmax25_anova_table, chlorophyll_mmol.m2_anova_table)
vcmax25_jmax25_chlorophyll_mmol.m2_anova_table_sub <- cbind(vcmax25_jmax25_chlorophyll_mmol.m2_anova_table[,2],
                                                            vcmax25_jmax25_chlorophyll_mmol.m2_anova_table[, c(1, 3, 4, 6, 7, 9)])
colnames(vcmax25_jmax25_chlorophyll_mmol.m2_anova_table_sub) <- c('df', 'χ2', 'P-value', 'χ2', 'P-value', 'χ2', 'P-value')
rownames(vcmax25_jmax25_chlorophyll_mmol.m2_anova_table_sub) <- c('Species (Sp)', 'Shading (Sh)', 'Fertilizer (F)',
                                                                  'Sp x Sh', 'Sp x F', 'Sh x F', 'Sp x Sh x F')
is.num <- sapply(vcmax25_jmax25_chlorophyll_mmol.m2_anova_table_sub, is.numeric)
vcmax25_jmax25_chlorophyll_mmol.m2_anova_table_sub[is.num] <- lapply(vcmax25_jmax25_chlorophyll_mmol.m2_anova_table_sub[is.num], 
                                                                     round, 3)
vcmax25_jmax25_chlorophyll_mmol.m2_anova_table_sub[vcmax25_jmax25_chlorophyll_mmol.m2_anova_table_sub<0.001] <- '<0.001'
write.csv(vcmax25_jmax25_chlorophyll_mmol.m2_anova_table_sub, 'tables/vcmax25_jmax25_chlorophyll_mmol.m2_anova_table.csv')

propN_rubisco_anova_table <- Anova(propN_rubisco_lmer)
propN_bioenergetics_anova_table <- Anova(propN_bioenergetics_lmer)
propN_lightharvesting_anova_table <- Anova(propN_lightharvesting_lmer)
propN_rubisco_propN_bioenergetics_propN_lightharvesting_anova_table <- cbind(propN_rubisco_anova_table, propN_bioenergetics_anova_table, propN_lightharvesting_anova_table)
propN_rubisco_propN_bioenergetics_propN_lightharvesting_anova_table_sub <- cbind(propN_rubisco_propN_bioenergetics_propN_lightharvesting_anova_table[,2],
                                                                                 propN_rubisco_propN_bioenergetics_propN_lightharvesting_anova_table[, c(1, 3, 4, 6, 7, 9)])
colnames(propN_rubisco_propN_bioenergetics_propN_lightharvesting_anova_table_sub) <- c('df', 'χ2', 'P-value', 'χ2', 'P-value', 'χ2', 'P-value')
rownames(propN_rubisco_propN_bioenergetics_propN_lightharvesting_anova_table_sub) <- c('Species (Sp)', 'Shading (Sh)', 'Fertilizer (F)',
                                                                                       'Sp x Sh', 'Sp x F', 'Sh x F', 'Sp x Sh x F')
is.num <- sapply(propN_rubisco_propN_bioenergetics_propN_lightharvesting_anova_table_sub, is.numeric)
propN_rubisco_propN_bioenergetics_propN_lightharvesting_anova_table_sub[is.num] <- lapply(propN_rubisco_propN_bioenergetics_propN_lightharvesting_anova_table_sub[is.num], 
                                                                                          round, 3)
propN_rubisco_propN_bioenergetics_propN_lightharvesting_anova_table_sub[propN_rubisco_propN_bioenergetics_propN_lightharvesting_anova_table_sub<0.001] <- '<0.001'
write.csv(propN_rubisco_propN_bioenergetics_propN_lightharvesting_anova_table_sub, 'tables/propN_rubisco_propN_bioenergetics_propN_lightharvesting_anova_table.csv')

tla_anova_table <- Anova(tla_lmer)
biomass_anova_table <- Anova(biomass_lmer)
tla_biomass_anova_table <- cbind(tla_anova_table, biomass_anova_table)
tla_biomass_anova_table_sub <- cbind(tla_biomass_anova_table[,2],
                                     tla_biomass_anova_table[, c(1, 3, 4, 6)])
colnames(tla_biomass_anova_table_sub) <- c('df', 'χ2', 'P-value', 'χ2', 'P-value')
rownames(tla_biomass_anova_table_sub) <- c('Species (Sp)', 'Shading (Sh)', 'Fertilizer (F)',
                                           'Sp x Sh', 'Sp x F', 'Sh x F', 'Sp x Sh x F')
is.num <- sapply(tla_biomass_anova_table_sub, is.numeric)
tla_biomass_anova_table_sub[is.num] <- lapply(tla_biomass_anova_table_sub[is.num], 
                                              round, 3)
tla_biomass_anova_table_sub[tla_biomass_anova_table_sub<0.001] <- '<0.001'
write.csv(tla_biomass_anova_table_sub, 'tables/tla_biomass_anova_table.csv')

### FIGURES

### marea
test(emtrends(marea_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 0)))
test(emtrends(marea_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 30)))
test(emtrends(marea_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 50)))
test(emtrends(marea_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 80)))

marea_emtrend_0 <- summary(emtrends(marea_lmer, ~spp, 
                                    var = 'fertilizer_cont', 
                                    at = list(shade_cont = 0)))
marea_emtrend_30 <- summary(emtrends(marea_lmer, ~spp, 
                                     var = 'fertilizer_cont', 
                                     at = list(shade_cont = 30)))
marea_emtrend_50 <- summary(emtrends(marea_lmer, ~spp, 
                                     var = 'fertilizer_cont', 
                                     at = list(shade_cont = 50)))
marea_emtrend_80 <- summary(emtrends(marea_lmer, ~spp, 
                                     var = 'fertilizer_cont', 
                                     at = list(shade_cont = 80)))
marea_intercept_0 <- summary(emmeans(marea_lmer, ~spp, 
                                     at = list(fertilizer_cont = 0, 
                                               shade_cont = 0)))
marea_intercept_30 <- summary(emmeans(marea_lmer, ~spp, 
                                      at = list(fertilizer_cont = 0, 
                                                shade_cont = 30)))
marea_intercept_50 <- summary(emmeans(marea_lmer, ~spp, 
                                      at = list(fertilizer_cont = 0, 
                                                shade_cont = 50)))
marea_intercept_80 <- summary(emmeans(marea_lmer, ~spp, 
                                      at = list(fertilizer_cont = 0, 
                                                shade_cont = 80)))

marea_func_cotton_0 <- function(x){
  exp(marea_emtrend_0[1, 2] * x + marea_intercept_0[1, 2])}
marea_func_cotton_30 <- function(x){
  exp(marea_emtrend_30[1, 2] * x + marea_intercept_30[1, 2])}
marea_func_cotton_50 <- function(x){
  exp(marea_emtrend_50[1, 2] * x + marea_intercept_50[1, 2])}
marea_func_cotton_80 <- function(x){
  exp(marea_emtrend_80[1, 2] * x + marea_intercept_80[1, 2])}

marea_func_soybean_0 <- function(x){
  exp(marea_emtrend_0[2, 2] * x + marea_intercept_0[2, 2])}
marea_func_soybean_30 <- function(x){
  exp(marea_emtrend_30[2, 2] * x + marea_intercept_30[2, 2])}
marea_func_soybean_50 <- function(x){
  exp(marea_emtrend_50[2, 2] * x + marea_intercept_50[2, 2])}
marea_func_soybean_80 <- function(x){
  exp(marea_emtrend_80[2, 2] * x + marea_intercept_80[2, 2])}

#### since there is no significant marea interction between shade and fertilizer, we plot a single line for each species
marea_emtrend <- summary(emtrends(marea_lmer, ~spp, 
                                  var = 'fertilizer_cont'))
marea_intercept <- summary(emmeans(marea_lmer, ~spp, 
                                   at = list(fertilizer_cont = 0)))

marea_func_cotton <- function(x){
  exp(marea_emtrend[1, 2] * x + marea_intercept[1, 2])}

marea_func_soybean <- function(x){
  exp(marea_emtrend[2, 2] * x + marea_intercept[2, 2])}

marea_cotton_plot <- ggplot(aes(y = marea, x = fertilizer_cont, color = shade), 
                            data = subset(data, spp == 'Cotton')) +
  theme(legend.position = "right",
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_jitter(width = 20, pch = 16, alpha = 0.5, size = 2) +
  scale_color_manual(values = c('red', 'orange', 'blue', 'purple')) +
  stat_function(fun = marea_func_cotton_0, color = 'red', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = marea_func_cotton_30, color = 'orange', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = marea_func_cotton_50, color = 'blue', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = marea_func_cotton_80, color = 'purple', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = marea_func_cotton, color = 'black', lwd = 2, lty = 2) +
  labs(color = 'Shade (%)') +
  xlab('Fertilizer (ppm)') +
  ylab(expression('M'[area] * ' (g m' ^ '-2' * ')')) +
  ggtitle(expression(italic('G. hirsutum'))) +
  labs(tag = "(a)") +
  ylim(c(10, 80))

marea_soybean_plot <- ggplot(aes(y = marea, x = fertilizer_cont, color = shade), 
                             data = subset(data, spp == 'Soybean')) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_jitter(width = 20, pch = 16, alpha = 0.5, size = 2) +
  scale_color_manual(values = c('red', 'orange', 'blue', 'purple')) +
  stat_function(fun = marea_func_soybean_0, color = 'red', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = marea_func_soybean_30, color = 'orange', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = marea_func_soybean_50, color = 'blue', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = marea_func_soybean_80, color = 'purple', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = marea_func_soybean, color = 'black', lwd = 2, lty = 2) +
  labs(color = 'Shade (%)') +
  xlab('Fertilizer (ppm)') +
  ylab(expression('M'[area] * ' (g m' ^ '-2' * ')'))+
  ggtitle(expression(italic('G. max'))) +
  labs(tag = "(b)")+
  ylim(c(10, 80))

### nmass
test(emtrends(nmass_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 0)))
test(emtrends(nmass_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 30)))
test(emtrends(nmass_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 50)))
test(emtrends(nmass_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 80)))

nmass_emtrend_0 <- summary(emtrends(nmass_lmer, ~spp, 
                                    var = 'fertilizer_cont', 
                                    at = list(shade_cont = 0)))
nmass_emtrend_30 <- summary(emtrends(nmass_lmer, ~spp, 
                                     var = 'fertilizer_cont', 
                                     at = list(shade_cont = 30)))
nmass_emtrend_50 <- summary(emtrends(nmass_lmer, ~spp, 
                                     var = 'fertilizer_cont', 
                                     at = list(shade_cont = 50)))
nmass_emtrend_80 <- summary(emtrends(nmass_lmer, ~spp, 
                                     var = 'fertilizer_cont', 
                                     at = list(shade_cont = 80)))
nmass_intercept_0 <- summary(emmeans(nmass_lmer, ~spp, 
                                     at = list(fertilizer_cont = 0, 
                                               shade_cont = 0)))
nmass_intercept_30 <- summary(emmeans(nmass_lmer, ~spp, 
                                      at = list(fertilizer_cont = 0, 
                                                shade_cont = 30)))
nmass_intercept_50 <- summary(emmeans(nmass_lmer, ~spp, 
                                      at = list(fertilizer_cont = 0, 
                                                shade_cont = 50)))
nmass_intercept_80 <- summary(emmeans(nmass_lmer, ~spp, 
                                      at = list(fertilizer_cont = 0, 
                                                shade_cont = 80)))

nmass_func_cotton_0 <- function(x){
  exp(nmass_emtrend_0[1, 2] * x + nmass_intercept_0[1, 2])}
nmass_func_cotton_30 <- function(x){
  exp(nmass_emtrend_30[1, 2] * x + nmass_intercept_30[1, 2])}
nmass_func_cotton_50 <- function(x){
  exp(nmass_emtrend_50[1, 2] * x + nmass_intercept_50[1, 2])}
nmass_func_cotton_80 <- function(x){
  exp(nmass_emtrend_80[1, 2] * x + nmass_intercept_80[1, 2])}

nmass_func_soybean_0 <- function(x){
  exp(nmass_emtrend_0[2, 2] * x + nmass_intercept_0[2, 2])}
nmass_func_soybean_30 <- function(x){
  exp(nmass_emtrend_30[2, 2] * x + nmass_intercept_30[2, 2])}
nmass_func_soybean_50 <- function(x){
  exp(nmass_emtrend_50[2, 2] * x + nmass_intercept_50[2, 2])}
nmass_func_soybean_80 <- function(x){
  exp(nmass_emtrend_80[2, 2] * x + nmass_intercept_80[2, 2])}

nmass_cotton_plot <- ggplot(aes(y = nmass, x = fertilizer_cont, color = shade), 
                            data = subset(data, spp == 'Cotton')) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_jitter(width = 20, pch = 16, alpha = 0.5, size = 2) +
  scale_color_manual(values = c('red', 'orange', 'blue', 'purple')) +
  stat_function(fun = nmass_func_cotton_0, color = 'red', lwd = 2, lty = 1) +
  stat_function(fun = nmass_func_cotton_30, color = 'orange', lwd = 2, lty = 1) +
  stat_function(fun = nmass_func_cotton_50, color = 'blue', lwd = 2, lty = 1) +
  stat_function(fun = nmass_func_cotton_80, color = 'purple', lwd = 2, lty = 1) +
  labs(color = 'Shade (%)') +
  xlab('Fertilizer (ppm)') +
  ylab(expression('N'[mass] * ' (g g' ^ '-1' * ')'))+
  labs(tag = "(c)")+
  ylim(c(0, 0.08))

nmass_soybean_plot <- ggplot(aes(y = nmass, x = fertilizer_cont, color = shade), 
                             data = subset(data, spp == 'Soybean')) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_jitter(width = 20, pch = 16, alpha = 0.5, size = 2) +
  scale_color_manual(values = c('red', 'orange', 'blue', 'purple')) +
  stat_function(fun = nmass_func_soybean_0, color = 'red', lwd = 2, lty = 2) +
  stat_function(fun = nmass_func_soybean_30, color = 'orange', lwd = 2, lty = 2) +
  stat_function(fun = nmass_func_soybean_50, color = 'blue', lwd = 2, lty = 2) +
  stat_function(fun = nmass_func_soybean_80, color = 'purple', lwd = 2, lty = 2) +
  labs(color = 'Shade (%)') +
  xlab('Fertilizer (ppm)') +
  ylab(expression('N'[mass] * ' (g g' ^ '-1' * ')'))+
  labs(tag = "(d)")+
  ylim(c(0, 0.08))

### narea

test(emtrends(narea_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 0)))
test(emtrends(narea_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 30)))
test(emtrends(narea_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 50)))
test(emtrends(narea_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 80)))

narea_emtrend_0 <- summary(emtrends(narea_lmer, ~spp, 
                                    var = 'fertilizer_cont', 
                                    at = list(shade_cont = 0)))
narea_emtrend_30 <- summary(emtrends(narea_lmer, ~spp, 
                                     var = 'fertilizer_cont', 
                                     at = list(shade_cont = 30)))
narea_emtrend_50 <- summary(emtrends(narea_lmer, ~spp, 
                                     var = 'fertilizer_cont', 
                                     at = list(shade_cont = 50)))
narea_emtrend_80 <- summary(emtrends(narea_lmer, ~spp, 
                                     var = 'fertilizer_cont', 
                                     at = list(shade_cont = 80)))
narea_intercept_0 <- summary(emmeans(narea_lmer, ~spp, 
                                     at = list(fertilizer_cont = 0, 
                                               shade_cont = 0)))
narea_intercept_30 <- summary(emmeans(narea_lmer, ~spp, 
                                      at = list(fertilizer_cont = 0, 
                                                shade_cont = 30)))
narea_intercept_50 <- summary(emmeans(narea_lmer, ~spp, 
                                      at = list(fertilizer_cont = 0, 
                                                shade_cont = 50)))
narea_intercept_80 <- summary(emmeans(narea_lmer, ~spp, 
                                      at = list(fertilizer_cont = 0, 
                                                shade_cont = 80)))

narea_func_cotton_0 <- function(x){
  exp(narea_emtrend_0[1, 2] * x + narea_intercept_0[1, 2])}
narea_func_cotton_30 <- function(x){
  exp(narea_emtrend_30[1, 2] * x + narea_intercept_30[1, 2])}
narea_func_cotton_50 <- function(x){
  exp(narea_emtrend_50[1, 2] * x + narea_intercept_50[1, 2])}
narea_func_cotton_80 <- function(x){
  exp(narea_emtrend_80[1, 2] * x + narea_intercept_80[1, 2])}

narea_func_soybean_0 <- function(x){
  exp(narea_emtrend_0[2, 2] * x + narea_intercept_0[2, 2])}
narea_func_soybean_30 <- function(x){
  exp(narea_emtrend_30[2, 2] * x + narea_intercept_30[2, 2])}
narea_func_soybean_50 <- function(x){
  exp(narea_emtrend_50[2, 2] * x + narea_intercept_50[2, 2])}
narea_func_soybean_80 <- function(x){
  exp(narea_emtrend_80[2, 2] * x + narea_intercept_80[2, 2])}

narea_cotton_plot <- ggplot(aes(y = narea, x = fertilizer_cont, color = shade), 
                            data = subset(data, spp == 'Cotton')) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_jitter(width = 20, pch = 16, alpha = 0.5, size = 2) +
  scale_color_manual(values = c('red', 'orange', 'blue', 'purple')) +
  stat_function(fun = narea_func_cotton_0, color = 'red', lwd = 2, lty = 1) +
  stat_function(fun = narea_func_cotton_30, color = 'orange', lwd = 2, lty = 1) +
  stat_function(fun = narea_func_cotton_50, color = 'blue', lwd = 2, lty = 1) +
  stat_function(fun = narea_func_cotton_80, color = 'purple', lwd = 2, lty = 1) +
  labs(color = 'Shade (%)') +
  xlab('Fertilizer (ppm)') +
  ylab(expression('N'[area] * ' (g m' ^ '-2' * ')'))+
  labs(tag = "(e)")+
  ylim(c(0.5, 3))

narea_soybean_plot <- ggplot(aes(y = narea, x = fertilizer_cont, color = shade), 
                             data = subset(data, spp == 'Soybean')) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_jitter(width = 20, pch = 16, alpha = 0.5, size = 2) +
  scale_color_manual(values = c('red', 'orange', 'blue', 'purple')) +
  stat_function(fun = narea_func_soybean_0, color = 'red', lwd = 2, lty = 2) +
  stat_function(fun = narea_func_soybean_30, color = 'orange', lwd = 2, lty = 1) +
  stat_function(fun = narea_func_soybean_50, color = 'blue', lwd = 2, lty = 1) +
  stat_function(fun = narea_func_soybean_80, color = 'purple', lwd = 2, lty = 2) +
  labs(color = 'Shade (%)') +
  xlab('Fertilizer (ppm)') +
  ylab(expression('N'[area] * ' (g m' ^ '-2' * ')'))+
  labs(tag = "(f)")+
  ylim(c(0.5, 3))

marea_cotton_plot_g <- ggplotGrob(marea_cotton_plot)
nmass_cotton_plot_g <- ggplotGrob(nmass_cotton_plot)
narea_cotton_plot_g <- ggplotGrob(narea_cotton_plot)
marea_nmass_narea_cotton_g <- rbind(marea_cotton_plot_g, 
                                    nmass_cotton_plot_g, 
                                    narea_cotton_plot_g,
                                    size = "max")
marea_soybean_plot_g <- ggplotGrob(marea_soybean_plot)
nmass_soybean_plot_g <- ggplotGrob(nmass_soybean_plot)
narea_soybean_plot_g <- ggplotGrob(narea_soybean_plot)
marea_nmass_narea_soybean_g <- rbind(marea_soybean_plot_g, 
                                     nmass_soybean_plot_g, 
                                     narea_soybean_plot_g,
                                     size = "max")
marea_nmass_narea_g <- cbind(marea_nmass_narea_cotton_g, 
                             marea_nmass_narea_soybean_g,
                             size = 'max')

jpeg(filename = "plots/marea_nmass_narea.jpeg", 
     width = 14, height = 14, units = 'in', res = 600)
grid.newpage()
grid.draw(marea_nmass_narea_g)
dev.off()

### vcmax25
test(emtrends(vcmax25_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 0)))
test(emtrends(vcmax25_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 30)))
test(emtrends(vcmax25_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 50)))
test(emtrends(vcmax25_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 80)))

vcmax25_emtrend_0 <- summary(emtrends(vcmax25_lmer, ~spp, 
                                      var = 'fertilizer_cont', 
                                      at = list(shade_cont = 0)))
vcmax25_emtrend_30 <- summary(emtrends(vcmax25_lmer, ~spp, 
                                       var = 'fertilizer_cont', 
                                       at = list(shade_cont = 30)))
vcmax25_emtrend_50 <- summary(emtrends(vcmax25_lmer, ~spp, 
                                       var = 'fertilizer_cont', 
                                       at = list(shade_cont = 50)))
vcmax25_emtrend_80 <- summary(emtrends(vcmax25_lmer, ~spp, 
                                       var = 'fertilizer_cont', 
                                       at = list(shade_cont = 80)))
vcmax25_intercept_0 <- summary(emmeans(vcmax25_lmer, ~spp, 
                                       at = list(fertilizer_cont = 0, 
                                                 shade_cont = 0)))
vcmax25_intercept_30 <- summary(emmeans(vcmax25_lmer, ~spp, 
                                        at = list(fertilizer_cont = 0, 
                                                  shade_cont = 30)))
vcmax25_intercept_50 <- summary(emmeans(vcmax25_lmer, ~spp, 
                                        at = list(fertilizer_cont = 0, 
                                                  shade_cont = 50)))
vcmax25_intercept_80 <- summary(emmeans(vcmax25_lmer, ~spp, 
                                        at = list(fertilizer_cont = 0, 
                                                  shade_cont = 80)))

vcmax25_func_cotton_0 <- function(x){
  exp(vcmax25_emtrend_0[1, 2] * x + vcmax25_intercept_0[1, 2])}
vcmax25_func_cotton_30 <- function(x){
  exp(vcmax25_emtrend_30[1, 2] * x + vcmax25_intercept_30[1, 2])}
vcmax25_func_cotton_50 <- function(x){
  exp(vcmax25_emtrend_50[1, 2] * x + vcmax25_intercept_50[1, 2])}
vcmax25_func_cotton_80 <- function(x){
  exp(vcmax25_emtrend_80[1, 2] * x + vcmax25_intercept_80[1, 2])}

vcmax25_func_soybean_0 <- function(x){
  exp(vcmax25_emtrend_0[2, 2] * x + vcmax25_intercept_0[2, 2])}
vcmax25_func_soybean_30 <- function(x){
  exp(vcmax25_emtrend_30[2, 2] * x + vcmax25_intercept_30[2, 2])}
vcmax25_func_soybean_50 <- function(x){
  exp(vcmax25_emtrend_50[2, 2] * x + vcmax25_intercept_50[2, 2])}
vcmax25_func_soybean_80 <- function(x){
  exp(vcmax25_emtrend_80[2, 2] * x + vcmax25_intercept_80[2, 2])}

#### since there is no significant vcmax25 interction between shade and fertilizer, we plot a single line for each species
vcmax25_emtrend <- summary(emtrends(vcmax25_lmer, ~spp, 
                                    var = 'fertilizer_cont'))
vcmax25_intercept <- summary(emmeans(vcmax25_lmer, ~spp, 
                                     at = list(fertilizer_cont = 0)))

vcmax25_func_cotton <- function(x){
  exp(vcmax25_emtrend[1, 2] * x + vcmax25_intercept[1, 2])}

vcmax25_func_soybean <- function(x){
  exp(vcmax25_emtrend[2, 2] * x + vcmax25_intercept[2, 2])}

vcmax25_cotton_plot <- ggplot(aes(y = vcmax25, x = fertilizer_cont, color = shade), 
                              data = subset(data, spp == 'Cotton')) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_jitter(width = 20, pch = 16, alpha = 0.5, size = 2) +
  scale_color_manual(values = c('red', 'orange', 'blue', 'purple')) +
  stat_function(fun = vcmax25_func_cotton_0, color = 'red', lwd = 2, lty = 1, alpha = 0.3) +
  stat_function(fun = vcmax25_func_cotton_30, color = 'orange', lwd = 2, lty = 1, alpha = 0.3) +
  stat_function(fun = vcmax25_func_cotton_50, color = 'blue', lwd = 2, lty = 1, alpha = 0.3) +
  stat_function(fun = vcmax25_func_cotton_80, color = 'purple', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = vcmax25_func_cotton, color = 'black', lwd = 2, lty = 1) +
  labs(color = 'Shade (%)') +
  xlab('Fertilizer (ppm)') +
  ylab(expression('V'[cmax25] * ' (µmol m' ^ '-2 ' * 's' ^ '-1' * ')')) +
  ggtitle(expression(italic('G. hirsutum'))) +
  labs(tag = "(a)") +
  ylim(c(10, 175))

vcmax25_soybean_plot <- ggplot(aes(y = vcmax25, x = fertilizer_cont, color = shade), 
                               data = subset(data, spp == 'Soybean')) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_jitter(width = 20, pch = 16, alpha = 0.5, size = 2) +
  scale_color_manual(values = c('red', 'orange', 'blue', 'purple')) +
  stat_function(fun = vcmax25_func_soybean_0, color = 'red', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = vcmax25_func_soybean_30, color = 'orange', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = vcmax25_func_soybean_50, color = 'blue', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = vcmax25_func_soybean_80, color = 'purple', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = vcmax25_func_soybean, color = 'black', lwd = 2, lty = 2) +
  labs(color = 'Shade (%)') +
  xlab('Fertilizer (ppm)') +
  ylab(expression('V'[cmax25] * ' (µmol m' ^ '-2 ' * 's' ^ '-1' * ')')) +
  ggtitle(expression(italic('G. max'))) +
  labs(tag = "(b)") +
  ylim(c(10, 175))

### jmax25
test(emtrends(jmax25_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 0)))
test(emtrends(jmax25_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 30)))
test(emtrends(jmax25_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 50)))
test(emtrends(jmax25_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 80)))

jmax25_emtrend_0 <- summary(emtrends(jmax25_lmer, ~spp, 
                                     var = 'fertilizer_cont', 
                                     at = list(shade_cont = 0)))
jmax25_emtrend_30 <- summary(emtrends(jmax25_lmer, ~spp, 
                                      var = 'fertilizer_cont', 
                                      at = list(shade_cont = 30)))
jmax25_emtrend_50 <- summary(emtrends(jmax25_lmer, ~spp, 
                                      var = 'fertilizer_cont', 
                                      at = list(shade_cont = 50)))
jmax25_emtrend_80 <- summary(emtrends(jmax25_lmer, ~spp, 
                                      var = 'fertilizer_cont', 
                                      at = list(shade_cont = 80)))
jmax25_intercept_0 <- summary(emmeans(jmax25_lmer, ~spp, 
                                      at = list(fertilizer_cont = 0, 
                                                shade_cont = 0)))
jmax25_intercept_30 <- summary(emmeans(jmax25_lmer, ~spp, 
                                       at = list(fertilizer_cont = 0, 
                                                 shade_cont = 30)))
jmax25_intercept_50 <- summary(emmeans(jmax25_lmer, ~spp, 
                                       at = list(fertilizer_cont = 0, 
                                                 shade_cont = 50)))
jmax25_intercept_80 <- summary(emmeans(jmax25_lmer, ~spp, 
                                       at = list(fertilizer_cont = 0, 
                                                 shade_cont = 80)))

jmax25_func_cotton_0 <- function(x){
  exp(jmax25_emtrend_0[1, 2] * x + jmax25_intercept_0[1, 2])}
jmax25_func_cotton_30 <- function(x){
  exp(jmax25_emtrend_30[1, 2] * x + jmax25_intercept_30[1, 2])}
jmax25_func_cotton_50 <- function(x){
  exp(jmax25_emtrend_50[1, 2] * x + jmax25_intercept_50[1, 2])}
jmax25_func_cotton_80 <- function(x){
  exp(jmax25_emtrend_80[1, 2] * x + jmax25_intercept_80[1, 2])}

jmax25_func_soybean_0 <- function(x){
  exp(jmax25_emtrend_0[2, 2] * x + jmax25_intercept_0[2, 2])}
jmax25_func_soybean_30 <- function(x){
  exp(jmax25_emtrend_30[2, 2] * x + jmax25_intercept_30[2, 2])}
jmax25_func_soybean_50 <- function(x){
  exp(jmax25_emtrend_50[2, 2] * x + jmax25_intercept_50[2, 2])}
jmax25_func_soybean_80 <- function(x){
  exp(jmax25_emtrend_80[2, 2] * x + jmax25_intercept_80[2, 2])}

#### since there is no significant jmax25 interction between shade and fertilizer, we plot a single line for each species
jmax25_emtrend <- summary(emtrends(jmax25_lmer, ~spp, 
                                   var = 'fertilizer_cont'))
jmax25_intercept <- summary(emmeans(jmax25_lmer, ~spp, 
                                    at = list(fertilizer_cont = 0)))

jmax25_func_cotton <- function(x){
  exp(jmax25_emtrend[1, 2] * x + jmax25_intercept[1, 2])}

jmax25_func_soybean <- function(x){
  exp(jmax25_emtrend[2, 2] * x + jmax25_intercept[2, 2])}

jmax25_cotton_plot <- ggplot(aes(y = jmax25, x = fertilizer_cont, color = shade), 
                             data = subset(data, spp == 'Cotton')) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_jitter(width = 20, pch = 16, alpha = 0.5, size = 2) +
  scale_color_manual(values = c('red', 'orange', 'blue', 'purple')) +
  stat_function(fun = jmax25_func_cotton_0, color = 'red', lwd = 2, lty = 1, alpha = 0.3) +
  stat_function(fun = jmax25_func_cotton_30, color = 'orange', lwd = 2, lty = 1, alpha = 0.3) +
  stat_function(fun = jmax25_func_cotton_50, color = 'blue', lwd = 2, lty = 1, alpha = 0.3) +
  stat_function(fun = jmax25_func_cotton_80, color = 'purple', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = jmax25_func_cotton, color = 'black', lwd = 2, lty = 1) +
  labs(color = 'Shade (%)') +
  xlab('Fertilizer (ppm)') +
  ylab(expression('J'[max25] * ' (µmol m' ^ '-2 ' * 's' ^ '-1' * ')')) +
  labs(tag = "(c)") +
  ylim(c(30, 210))

jmax25_soybean_plot <- ggplot(aes(y = jmax25, x = fertilizer_cont, color = shade), 
                              data = subset(data, spp == 'Soybean')) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_jitter(width = 20, pch = 16, alpha = 0.5, size = 2) +
  scale_color_manual(values = c('red', 'orange', 'blue', 'purple')) +
  stat_function(fun = jmax25_func_soybean_0, color = 'red', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = jmax25_func_soybean_30, color = 'orange', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = jmax25_func_soybean_50, color = 'blue', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = jmax25_func_soybean_80, color = 'purple', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = jmax25_func_soybean, color = 'black', lwd = 2, lty = 2) +
  labs(color = 'Shade (%)') +
  xlab('Fertilizer (ppm)') +
  ylab(expression('J'[max25] * ' (µmol m' ^ '-2 ' * 's' ^ '-1' * ')')) +
  labs(tag = "(d)") +
  ylim(c(30, 210))

### chlorophyll_mmol.m2

test(emtrends(chlorophyll_mmol.m2_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 0)))
test(emtrends(chlorophyll_mmol.m2_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 30)))
test(emtrends(chlorophyll_mmol.m2_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 50)))
test(emtrends(chlorophyll_mmol.m2_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 80)))

chlorophyll_mmol.m2_emtrend_0 <- summary(emtrends(chlorophyll_mmol.m2_lmer, ~spp, 
                                                  var = 'fertilizer_cont', 
                                                  at = list(shade_cont = 0)))
chlorophyll_mmol.m2_emtrend_30 <- summary(emtrends(chlorophyll_mmol.m2_lmer, ~spp, 
                                                   var = 'fertilizer_cont', 
                                                   at = list(shade_cont = 30)))
chlorophyll_mmol.m2_emtrend_50 <- summary(emtrends(chlorophyll_mmol.m2_lmer, ~spp, 
                                                   var = 'fertilizer_cont', 
                                                   at = list(shade_cont = 50)))
chlorophyll_mmol.m2_emtrend_80 <- summary(emtrends(chlorophyll_mmol.m2_lmer, ~spp, 
                                                   var = 'fertilizer_cont', 
                                                   at = list(shade_cont = 80)))
chlorophyll_mmol.m2_intercept_0 <- summary(emmeans(chlorophyll_mmol.m2_lmer, ~spp, 
                                                   at = list(fertilizer_cont = 0, 
                                                             shade_cont = 0)))
chlorophyll_mmol.m2_intercept_30 <- summary(emmeans(chlorophyll_mmol.m2_lmer, ~spp, 
                                                    at = list(fertilizer_cont = 0, 
                                                              shade_cont = 30)))
chlorophyll_mmol.m2_intercept_50 <- summary(emmeans(chlorophyll_mmol.m2_lmer, ~spp, 
                                                    at = list(fertilizer_cont = 0, 
                                                              shade_cont = 50)))
chlorophyll_mmol.m2_intercept_80 <- summary(emmeans(chlorophyll_mmol.m2_lmer, ~spp, 
                                                    at = list(fertilizer_cont = 0, 
                                                              shade_cont = 80)))

chlorophyll_mmol.m2_func_cotton_0 <- function(x){
  (chlorophyll_mmol.m2_emtrend_0[1, 2] * x + chlorophyll_mmol.m2_intercept_0[1, 2])}
chlorophyll_mmol.m2_func_cotton_30 <- function(x){
  (chlorophyll_mmol.m2_emtrend_30[1, 2] * x + chlorophyll_mmol.m2_intercept_30[1, 2])}
chlorophyll_mmol.m2_func_cotton_50 <- function(x){
  (chlorophyll_mmol.m2_emtrend_50[1, 2] * x + chlorophyll_mmol.m2_intercept_50[1, 2])}
chlorophyll_mmol.m2_func_cotton_80 <- function(x){
  (chlorophyll_mmol.m2_emtrend_80[1, 2] * x + chlorophyll_mmol.m2_intercept_80[1, 2])}

chlorophyll_mmol.m2_func_soybean_0 <- function(x){
  (chlorophyll_mmol.m2_emtrend_0[2, 2] * x + chlorophyll_mmol.m2_intercept_0[2, 2])}
chlorophyll_mmol.m2_func_soybean_30 <- function(x){
  (chlorophyll_mmol.m2_emtrend_30[2, 2] * x + chlorophyll_mmol.m2_intercept_30[2, 2])}
chlorophyll_mmol.m2_func_soybean_50 <- function(x){
  (chlorophyll_mmol.m2_emtrend_50[2, 2] * x + chlorophyll_mmol.m2_intercept_50[2, 2])}
chlorophyll_mmol.m2_func_soybean_80 <- function(x){
  (chlorophyll_mmol.m2_emtrend_80[2, 2] * x + chlorophyll_mmol.m2_intercept_80[2, 2])}

#### since there is no significant chlorophyll_mmol.m2 interction between shade and fertilizer, we plot a single line for each species
chlorophyll_mmol.m2_emtrend <- summary(emtrends(chlorophyll_mmol.m2_lmer, ~spp, 
                                                var = 'fertilizer_cont'))
chlorophyll_mmol.m2_intercept <- summary(emmeans(chlorophyll_mmol.m2_lmer, ~spp, 
                                                 at = list(fertilizer_cont = 0)))

chlorophyll_mmol.m2_func_cotton <- function(x){
  (chlorophyll_mmol.m2_emtrend[1, 2] * x + chlorophyll_mmol.m2_intercept[1, 2])}

chlorophyll_mmol.m2_func_soybean <- function(x){
  (chlorophyll_mmol.m2_emtrend[2, 2] * x + chlorophyll_mmol.m2_intercept[2, 2])}

chlorophyll_mmol.m2_cotton_plot <- ggplot(aes(y = chl_mmol.m2, x = fertilizer_cont, color = shade), 
                                          data = subset(data, spp == 'Cotton')) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_jitter(width = 20, pch = 16, alpha = 0.5, size = 2) +
  scale_color_manual(values = c('red', 'orange', 'blue', 'purple')) +
  stat_function(fun = chlorophyll_mmol.m2_func_cotton_0, color = 'red', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = chlorophyll_mmol.m2_func_cotton_30, color = 'orange', lwd = 2, lty = 1, alpha = 0.3) +
  stat_function(fun = chlorophyll_mmol.m2_func_cotton_50, color = 'blue', lwd = 2, lty = 1, alpha = 0.3) +
  stat_function(fun = chlorophyll_mmol.m2_func_cotton_80, color = 'purple', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = chlorophyll_mmol.m2_func_cotton, color = 'black', lwd = 2, lty = 1) +
  labs(color = 'Shade (%)') +
  xlab('Fertilizer (ppm)') +
  ylab(expression('Chl'[area] * ' (mmol m' ^ '-2 ' * ')')) +
  labs(tag = "(e)") +
  ylim(c(0, 1))

chlorophyll_mmol.m2_soybean_plot <- ggplot(aes(y = chl_mmol.m2, x = fertilizer_cont, color = shade), 
                                           data = subset(data, spp == 'Soybean')) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_jitter(width = 20, pch = 16, alpha = 0.5, size = 2) +
  scale_color_manual(values = c('red', 'orange', 'blue', 'purple')) +
  stat_function(fun = chlorophyll_mmol.m2_func_soybean_0, color = 'red', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = chlorophyll_mmol.m2_func_soybean_30, color = 'orange', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = chlorophyll_mmol.m2_func_soybean_50, color = 'blue', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = chlorophyll_mmol.m2_func_soybean_80, color = 'purple', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = chlorophyll_mmol.m2_func_soybean, color = 'black', lwd = 2, lty = 2) +
  labs(color = 'Shade (%)') +
  xlab('Fertilizer (ppm)') +
  ylab(expression('Chl'[area] * ' (mmol m' ^ '-2 ' * ')')) +
  labs(tag = "(f)") +
  ylim(c(0, 1))

vcmax25_cotton_plot_g <- ggplotGrob(vcmax25_cotton_plot)
jmax25_cotton_plot_g <- ggplotGrob(jmax25_cotton_plot)
chlorophyll_mmol.m2_cotton_plot_g <- ggplotGrob(chlorophyll_mmol.m2_cotton_plot)
vcmax25_jmax25_chlorophyll_mmol.m2_cotton_g <- rbind(vcmax25_cotton_plot_g, 
                                                     jmax25_cotton_plot_g, 
                                                     chlorophyll_mmol.m2_cotton_plot_g,
                                                     size = "max")
vcmax25_soybean_plot_g <- ggplotGrob(vcmax25_soybean_plot)
jmax25_soybean_plot_g <- ggplotGrob(jmax25_soybean_plot)
chlorophyll_mmol.m2_soybean_plot_g <- ggplotGrob(chlorophyll_mmol.m2_soybean_plot)
vcmax25_jmax25_chlorophyll_mmol.m2_soybean_g <- rbind(vcmax25_soybean_plot_g, 
                                                      jmax25_soybean_plot_g, 
                                                      chlorophyll_mmol.m2_soybean_plot_g,
                                                      size = "max")
vcmax25_jmax25_chlorophyll_mmol.m2_g <- cbind(vcmax25_jmax25_chlorophyll_mmol.m2_cotton_g, 
                                              vcmax25_jmax25_chlorophyll_mmol.m2_soybean_g,
                                              size = 'max')

jpeg(filename = "plots/vcmax25_jmax25_chlorophyll_mmol.m2.jpeg", 
     width = 14, height = 14, units = 'in', res = 600)
grid.newpage()
grid.draw(vcmax25_jmax25_chlorophyll_mmol.m2_g)
dev.off()

### propN_rubisco
test(emtrends(propN_rubisco_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 0)))
test(emtrends(propN_rubisco_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 30)))
test(emtrends(propN_rubisco_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 50)))
test(emtrends(propN_rubisco_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 80)))

propN_rubisco_emtrend_0 <- summary(emtrends(propN_rubisco_lmer, ~spp, 
                                            var = 'fertilizer_cont', 
                                            at = list(shade_cont = 0)))
propN_rubisco_emtrend_30 <- summary(emtrends(propN_rubisco_lmer, ~spp, 
                                             var = 'fertilizer_cont', 
                                             at = list(shade_cont = 30)))
propN_rubisco_emtrend_50 <- summary(emtrends(propN_rubisco_lmer, ~spp, 
                                             var = 'fertilizer_cont', 
                                             at = list(shade_cont = 50)))
propN_rubisco_emtrend_80 <- summary(emtrends(propN_rubisco_lmer, ~spp, 
                                             var = 'fertilizer_cont', 
                                             at = list(shade_cont = 80)))
propN_rubisco_intercept_0 <- summary(emmeans(propN_rubisco_lmer, ~spp, 
                                             at = list(fertilizer_cont = 0, 
                                                       shade_cont = 0)))
propN_rubisco_intercept_30 <- summary(emmeans(propN_rubisco_lmer, ~spp, 
                                              at = list(fertilizer_cont = 0, 
                                                        shade_cont = 30)))
propN_rubisco_intercept_50 <- summary(emmeans(propN_rubisco_lmer, ~spp, 
                                              at = list(fertilizer_cont = 0, 
                                                        shade_cont = 50)))
propN_rubisco_intercept_80 <- summary(emmeans(propN_rubisco_lmer, ~spp, 
                                              at = list(fertilizer_cont = 0, 
                                                        shade_cont = 80)))

propN_rubisco_func_cotton_0 <- function(x){
  (propN_rubisco_emtrend_0[1, 2] * x + propN_rubisco_intercept_0[1, 2])}
propN_rubisco_func_cotton_30 <- function(x){
  (propN_rubisco_emtrend_30[1, 2] * x + propN_rubisco_intercept_30[1, 2])}
propN_rubisco_func_cotton_50 <- function(x){
  (propN_rubisco_emtrend_50[1, 2] * x + propN_rubisco_intercept_50[1, 2])}
propN_rubisco_func_cotton_80 <- function(x){
  (propN_rubisco_emtrend_80[1, 2] * x + propN_rubisco_intercept_80[1, 2])}

propN_rubisco_func_soybean_0 <- function(x){
  (propN_rubisco_emtrend_0[2, 2] * x + propN_rubisco_intercept_0[2, 2])}
propN_rubisco_func_soybean_30 <- function(x){
  (propN_rubisco_emtrend_30[2, 2] * x + propN_rubisco_intercept_30[2, 2])}
propN_rubisco_func_soybean_50 <- function(x){
  (propN_rubisco_emtrend_50[2, 2] * x + propN_rubisco_intercept_50[2, 2])}
propN_rubisco_func_soybean_80 <- function(x){
  (propN_rubisco_emtrend_80[2, 2] * x + propN_rubisco_intercept_80[2, 2])}

#### since there is no significant propN_rubisco interction between shade and fertilizer, we plot a single line for each species
propN_rubisco_emtrend <- summary(emtrends(propN_rubisco_lmer, ~spp, 
                                          var = 'fertilizer_cont'))
propN_rubisco_intercept <- summary(emmeans(propN_rubisco_lmer, ~spp, 
                                           at = list(fertilizer_cont = 0)))

propN_rubisco_func_cotton <- function(x){
  (propN_rubisco_emtrend[1, 2] * x + propN_rubisco_intercept[1, 2])}

propN_rubisco_func_soybean <- function(x){
  (propN_rubisco_emtrend[2, 2] * x + propN_rubisco_intercept[2, 2])}

propN_rubisco_cotton_plot <- ggplot(aes(y = propN_rubisco, x = fertilizer_cont, color = shade), 
                                    data = subset(data, spp == 'Cotton')) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_jitter(width = 20, pch = 16, alpha = 0.5, size = 2) +
  scale_color_manual(values = c('red', 'orange', 'blue', 'purple')) +
  stat_function(fun = propN_rubisco_func_cotton_0, color = 'red', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = propN_rubisco_func_cotton_30, color = 'orange', lwd = 2, lty = 1, alpha = 0.3) +
  stat_function(fun = propN_rubisco_func_cotton_50, color = 'blue', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = propN_rubisco_func_cotton_80, color = 'purple', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = propN_rubisco_func_cotton, color = 'black', lwd = 2, lty = 1) +
  labs(color = 'Shade (%)') +
  xlab('Fertilizer (ppm)') +
  ylab(expression('ρ'[rubisco] * ' (gN gN' ^ '-1' * ')')) +
  ggtitle(expression(italic('G. hirsutum'))) +
  labs(tag = "(a)") +
  ylim(c(0, 0.9))

propN_rubisco_soybean_plot <- ggplot(aes(y = propN_rubisco, x = fertilizer_cont, color = shade), 
                                     data = subset(data, spp == 'Soybean')) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_jitter(width = 20, pch = 16, alpha = 0.5, size = 2) +
  scale_color_manual(values = c('red', 'orange', 'blue', 'purple')) +
  stat_function(fun = propN_rubisco_func_soybean_0, color = 'red', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = propN_rubisco_func_soybean_30, color = 'orange', lwd = 2, lty = 1, alpha = 0.3) +
  stat_function(fun = propN_rubisco_func_soybean_50, color = 'blue', lwd = 2, lty = 1, alpha = 0.3) +
  stat_function(fun = propN_rubisco_func_soybean_80, color = 'purple', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = propN_rubisco_func_soybean, color = 'black', lwd = 2, lty = 1) +
  labs(color = 'Shade (%)') +
  xlab('Fertilizer (ppm)') +
  ylab(expression('ρ'[rubisco] * ' (gN gN' ^ '-1' * ')')) +
  ggtitle(expression(italic('G. max'))) +
  labs(tag = "(b)") +
  ylim(c(0, 0.9))

### propN_bioenergetics
test(emtrends(propN_bioenergetics_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 0)))
test(emtrends(propN_bioenergetics_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 30)))
test(emtrends(propN_bioenergetics_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 50)))
test(emtrends(propN_bioenergetics_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 80)))

propN_bioenergetics_emtrend_0 <- summary(emtrends(propN_bioenergetics_lmer, ~spp, 
                                                  var = 'fertilizer_cont', 
                                                  at = list(shade_cont = 0)))
propN_bioenergetics_emtrend_30 <- summary(emtrends(propN_bioenergetics_lmer, ~spp, 
                                                   var = 'fertilizer_cont', 
                                                   at = list(shade_cont = 30)))
propN_bioenergetics_emtrend_50 <- summary(emtrends(propN_bioenergetics_lmer, ~spp, 
                                                   var = 'fertilizer_cont', 
                                                   at = list(shade_cont = 50)))
propN_bioenergetics_emtrend_80 <- summary(emtrends(propN_bioenergetics_lmer, ~spp, 
                                                   var = 'fertilizer_cont', 
                                                   at = list(shade_cont = 80)))
propN_bioenergetics_intercept_0 <- summary(emmeans(propN_bioenergetics_lmer, ~spp, 
                                                   at = list(fertilizer_cont = 0, 
                                                             shade_cont = 0)))
propN_bioenergetics_intercept_30 <- summary(emmeans(propN_bioenergetics_lmer, ~spp, 
                                                    at = list(fertilizer_cont = 0, 
                                                              shade_cont = 30)))
propN_bioenergetics_intercept_50 <- summary(emmeans(propN_bioenergetics_lmer, ~spp, 
                                                    at = list(fertilizer_cont = 0, 
                                                              shade_cont = 50)))
propN_bioenergetics_intercept_80 <- summary(emmeans(propN_bioenergetics_lmer, ~spp, 
                                                    at = list(fertilizer_cont = 0, 
                                                              shade_cont = 80)))

propN_bioenergetics_func_cotton_0 <- function(x){
  (propN_bioenergetics_emtrend_0[1, 2] * x + propN_bioenergetics_intercept_0[1, 2])}
propN_bioenergetics_func_cotton_30 <- function(x){
  (propN_bioenergetics_emtrend_30[1, 2] * x + propN_bioenergetics_intercept_30[1, 2])}
propN_bioenergetics_func_cotton_50 <- function(x){
  (propN_bioenergetics_emtrend_50[1, 2] * x + propN_bioenergetics_intercept_50[1, 2])}
propN_bioenergetics_func_cotton_80 <- function(x){
  (propN_bioenergetics_emtrend_80[1, 2] * x + propN_bioenergetics_intercept_80[1, 2])}

propN_bioenergetics_func_soybean_0 <- function(x){
  (propN_bioenergetics_emtrend_0[2, 2] * x + propN_bioenergetics_intercept_0[2, 2])}
propN_bioenergetics_func_soybean_30 <- function(x){
  (propN_bioenergetics_emtrend_30[2, 2] * x + propN_bioenergetics_intercept_30[2, 2])}
propN_bioenergetics_func_soybean_50 <- function(x){
  (propN_bioenergetics_emtrend_50[2, 2] * x + propN_bioenergetics_intercept_50[2, 2])}
propN_bioenergetics_func_soybean_80 <- function(x){
  (propN_bioenergetics_emtrend_80[2, 2] * x + propN_bioenergetics_intercept_80[2, 2])}

#### since there is no significant propN_bioenergetics interction between shade and fertilizer, we plot a single line for each species
propN_bioenergetics_emtrend <- summary(emtrends(propN_bioenergetics_lmer, ~spp, 
                                                var = 'fertilizer_cont'))
propN_bioenergetics_intercept <- summary(emmeans(propN_bioenergetics_lmer, ~spp, 
                                                 at = list(fertilizer_cont = 0)))

propN_bioenergetics_func_cotton <- function(x){
  (propN_bioenergetics_emtrend[1, 2] * x + propN_bioenergetics_intercept[1, 2])}

propN_bioenergetics_func_soybean <- function(x){
  (propN_bioenergetics_emtrend[2, 2] * x + propN_bioenergetics_intercept[2, 2])}

propN_bioenergetics_cotton_plot <- ggplot(aes(y = propN_bioenergetics, x = fertilizer_cont, color = shade), 
                                          data = subset(data, spp == 'Cotton')) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_jitter(width = 20, pch = 16, alpha = 0.5, size = 2) +
  scale_color_manual(values = c('red', 'orange', 'blue', 'purple')) +
  stat_function(fun = propN_bioenergetics_func_cotton_0, color = 'red', lwd = 2, lty = 1, alpha = 0.3) +
  stat_function(fun = propN_bioenergetics_func_cotton_30, color = 'orange', lwd = 2, lty = 1, alpha = 0.3) +
  stat_function(fun = propN_bioenergetics_func_cotton_50, color = 'blue', lwd = 2, lty = 1, alpha = 0.3) +
  stat_function(fun = propN_bioenergetics_func_cotton_80, color = 'purple', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = propN_bioenergetics_func_cotton, color = 'black', lwd = 2, lty = 1) +
  labs(color = 'Shade (%)') +
  xlab('Fertilizer (ppm)') +
  ylab(expression('ρ'[bioenergetics] * ' (gN gN' ^ '-1' * ')')) +
  labs(tag = "(c)") +
  ylim(c(0.02, 0.13))

propN_bioenergetics_soybean_plot <- ggplot(aes(y = propN_bioenergetics, x = fertilizer_cont, color = shade), 
                                           data = subset(data, spp == 'Soybean')) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_jitter(width = 20, pch = 16, alpha = 0.5, size = 2) +
  scale_color_manual(values = c('red', 'orange', 'blue', 'purple')) +
  stat_function(fun = propN_bioenergetics_func_soybean_0, color = 'red', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = propN_bioenergetics_func_soybean_30, color = 'orange', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = propN_bioenergetics_func_soybean_50, color = 'blue', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = propN_bioenergetics_func_soybean_80, color = 'purple', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = propN_bioenergetics_func_soybean, color = 'black', lwd = 2, lty = 1) +
  labs(color = 'Shade (%)') +
  xlab('Fertilizer (ppm)') +
  ylab(expression('ρ'[bioenergetics] * ' (gN gN' ^ '-1' * ')')) +
  labs(tag = "(d)") +
  ylim(c(0.02, 0.13))

### propN_lightharvesting

test(emtrends(propN_lightharvesting_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 0)))
test(emtrends(propN_lightharvesting_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 30)))
test(emtrends(propN_lightharvesting_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 50)))
test(emtrends(propN_lightharvesting_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 80)))

propN_lightharvesting_emtrend_0 <- summary(emtrends(propN_lightharvesting_lmer, ~spp, 
                                                    var = 'fertilizer_cont', 
                                                    at = list(shade_cont = 0)))
propN_lightharvesting_emtrend_30 <- summary(emtrends(propN_lightharvesting_lmer, ~spp, 
                                                     var = 'fertilizer_cont', 
                                                     at = list(shade_cont = 30)))
propN_lightharvesting_emtrend_50 <- summary(emtrends(propN_lightharvesting_lmer, ~spp, 
                                                     var = 'fertilizer_cont', 
                                                     at = list(shade_cont = 50)))
propN_lightharvesting_emtrend_80 <- summary(emtrends(propN_lightharvesting_lmer, ~spp, 
                                                     var = 'fertilizer_cont', 
                                                     at = list(shade_cont = 80)))
propN_lightharvesting_intercept_0 <- summary(emmeans(propN_lightharvesting_lmer, ~spp, 
                                                     at = list(fertilizer_cont = 0, 
                                                               shade_cont = 0)))
propN_lightharvesting_intercept_30 <- summary(emmeans(propN_lightharvesting_lmer, ~spp, 
                                                      at = list(fertilizer_cont = 0, 
                                                                shade_cont = 30)))
propN_lightharvesting_intercept_50 <- summary(emmeans(propN_lightharvesting_lmer, ~spp, 
                                                      at = list(fertilizer_cont = 0, 
                                                                shade_cont = 50)))
propN_lightharvesting_intercept_80 <- summary(emmeans(propN_lightharvesting_lmer, ~spp, 
                                                      at = list(fertilizer_cont = 0, 
                                                                shade_cont = 80)))

propN_lightharvesting_func_cotton_0 <- function(x){
  (propN_lightharvesting_emtrend_0[1, 2] * x + propN_lightharvesting_intercept_0[1, 2])}
propN_lightharvesting_func_cotton_30 <- function(x){
  (propN_lightharvesting_emtrend_30[1, 2] * x + propN_lightharvesting_intercept_30[1, 2])}
propN_lightharvesting_func_cotton_50 <- function(x){
  (propN_lightharvesting_emtrend_50[1, 2] * x + propN_lightharvesting_intercept_50[1, 2])}
propN_lightharvesting_func_cotton_80 <- function(x){
  (propN_lightharvesting_emtrend_80[1, 2] * x + propN_lightharvesting_intercept_80[1, 2])}

propN_lightharvesting_func_soybean_0 <- function(x){
  (propN_lightharvesting_emtrend_0[2, 2] * x + propN_lightharvesting_intercept_0[2, 2])}
propN_lightharvesting_func_soybean_30 <- function(x){
  (propN_lightharvesting_emtrend_30[2, 2] * x + propN_lightharvesting_intercept_30[2, 2])}
propN_lightharvesting_func_soybean_50 <- function(x){
  (propN_lightharvesting_emtrend_50[2, 2] * x + propN_lightharvesting_intercept_50[2, 2])}
propN_lightharvesting_func_soybean_80 <- function(x){
  (propN_lightharvesting_emtrend_80[2, 2] * x + propN_lightharvesting_intercept_80[2, 2])}

#### since there is no significant propN_lightharvesting interction between shade and fertilizer, we plot a single line for each species
propN_lightharvesting_emtrend <- summary(emtrends(propN_lightharvesting_lmer, ~spp, 
                                                  var = 'fertilizer_cont'))
propN_lightharvesting_intercept <- summary(emmeans(propN_lightharvesting_lmer, ~spp, 
                                                   at = list(fertilizer_cont = 0)))

propN_lightharvesting_func_cotton <- function(x){
  (propN_lightharvesting_emtrend[1, 2] * x + propN_lightharvesting_intercept[1, 2])}

propN_lightharvesting_func_soybean <- function(x){
  (propN_lightharvesting_emtrend[2, 2] * x + propN_lightharvesting_intercept[2, 2])}

propN_lightharvesting_cotton_plot <- ggplot(aes(y = propN_lightharvesting, x = fertilizer_cont, color = shade), 
                                            data = subset(data, spp == 'Cotton')) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_jitter(width = 20, pch = 16, alpha = 0.5, size = 2) +
  scale_color_manual(values = c('red', 'orange', 'blue', 'purple')) +
  stat_function(fun = propN_lightharvesting_func_cotton_0, color = 'red', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = propN_lightharvesting_func_cotton_30, color = 'orange', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = propN_lightharvesting_func_cotton_50, color = 'blue', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = propN_lightharvesting_func_cotton_80, color = 'purple', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = propN_lightharvesting_func_cotton, color = 'black', lwd = 2, lty = 1) +
  labs(color = 'Shade (%)') +
  xlab('Fertilizer (ppm)') +
  ylab(expression('ρ'[lightharvesting] * ' (gN gN' ^ '-1' * ')')) +
  labs(tag = "(e)") +
  ylim(c(0, 0.1))

propN_lightharvesting_soybean_plot <- ggplot(aes(y = propN_lightharvesting, x = fertilizer_cont, color = shade), 
                                             data = subset(data, spp == 'Soybean')) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_jitter(width = 20, pch = 16, alpha = 0.5, size = 2) +
  scale_color_manual(values = c('red', 'orange', 'blue', 'purple')) +
  stat_function(fun = propN_lightharvesting_func_soybean_0, color = 'red', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = propN_lightharvesting_func_soybean_30, color = 'orange', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = propN_lightharvesting_func_soybean_50, color = 'blue', lwd = 2, lty = 1, alpha = 0.3) +
  stat_function(fun = propN_lightharvesting_func_soybean_80, color = 'purple', lwd = 2, lty = 2, alpha = 0.3) +
  stat_function(fun = propN_lightharvesting_func_soybean, color = 'black', lwd = 2, lty = 1) +
  labs(color = 'Shade (%)') +
  xlab('Fertilizer (ppm)') +
  ylab(expression('ρ'[lightharvesting] * ' (gN gN' ^ '-1' * ')')) +
  labs(tag = "(f)") +
  ylim(c(0, 0.1))

propN_rubisco_cotton_plot_g <- ggplotGrob(propN_rubisco_cotton_plot)
propN_bioenergetics_cotton_plot_g <- ggplotGrob(propN_bioenergetics_cotton_plot)
propN_lightharvesting_cotton_plot_g <- ggplotGrob(propN_lightharvesting_cotton_plot)
propN_rubisco_propN_bioenergetics_propN_lightharvesting_cotton_g <- rbind(propN_rubisco_cotton_plot_g, 
                                                                          propN_bioenergetics_cotton_plot_g, 
                                                                          propN_lightharvesting_cotton_plot_g,
                                                                          size = "max")
propN_rubisco_soybean_plot_g <- ggplotGrob(propN_rubisco_soybean_plot)
propN_bioenergetics_soybean_plot_g <- ggplotGrob(propN_bioenergetics_soybean_plot)
propN_lightharvesting_soybean_plot_g <- ggplotGrob(propN_lightharvesting_soybean_plot)
propN_rubisco_propN_bioenergetics_propN_lightharvesting_soybean_g <- rbind(propN_rubisco_soybean_plot_g, 
                                                                           propN_bioenergetics_soybean_plot_g, 
                                                                           propN_lightharvesting_soybean_plot_g,
                                                                           size = "max")
propN_rubisco_propN_bioenergetics_propN_lightharvesting_g <- cbind(propN_rubisco_propN_bioenergetics_propN_lightharvesting_cotton_g, 
                                                                   propN_rubisco_propN_bioenergetics_propN_lightharvesting_soybean_g,
                                                                   size = 'max')

jpeg(filename = "plots/propN_rubisco_propN_bioenergetics_propN_lightharvesting.jpeg", 
     width = 14, height = 14, units = 'in', res = 600)
grid.newpage()
grid.draw(propN_rubisco_propN_bioenergetics_propN_lightharvesting_g)
dev.off()

### tla
test(emtrends(tla_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 0)))
test(emtrends(tla_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 30)))
test(emtrends(tla_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 50)))
test(emtrends(tla_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 80)))

tla_emtrend_0 <- summary(emtrends(tla_lmer, ~spp, 
                                  var = 'fertilizer_cont', 
                                  at = list(shade_cont = 0)))
tla_emtrend_30 <- summary(emtrends(tla_lmer, ~spp, 
                                   var = 'fertilizer_cont', 
                                   at = list(shade_cont = 30)))
tla_emtrend_50 <- summary(emtrends(tla_lmer, ~spp, 
                                   var = 'fertilizer_cont', 
                                   at = list(shade_cont = 50)))
tla_emtrend_80 <- summary(emtrends(tla_lmer, ~spp, 
                                   var = 'fertilizer_cont', 
                                   at = list(shade_cont = 80)))
tla_intercept_0 <- summary(emmeans(tla_lmer, ~spp, 
                                   at = list(fertilizer_cont = 0, 
                                             shade_cont = 0)))
tla_intercept_30 <- summary(emmeans(tla_lmer, ~spp, 
                                    at = list(fertilizer_cont = 0, 
                                              shade_cont = 30)))
tla_intercept_50 <- summary(emmeans(tla_lmer, ~spp, 
                                    at = list(fertilizer_cont = 0, 
                                              shade_cont = 50)))
tla_intercept_80 <- summary(emmeans(tla_lmer, ~spp, 
                                    at = list(fertilizer_cont = 0, 
                                              shade_cont = 80)))

tla_func_cotton_0 <- function(x){
  exp(tla_emtrend_0[1, 2] * x + tla_intercept_0[1, 2])}
tla_func_cotton_30 <- function(x){
  exp(tla_emtrend_30[1, 2] * x + tla_intercept_30[1, 2])}
tla_func_cotton_50 <- function(x){
  exp(tla_emtrend_50[1, 2] * x + tla_intercept_50[1, 2])}
tla_func_cotton_80 <- function(x){
  exp(tla_emtrend_80[1, 2] * x + tla_intercept_80[1, 2])}

tla_func_soybean_0 <- function(x){
  exp(tla_emtrend_0[2, 2] * x + tla_intercept_0[2, 2])}
tla_func_soybean_30 <- function(x){
  exp(tla_emtrend_30[2, 2] * x + tla_intercept_30[2, 2])}
tla_func_soybean_50 <- function(x){
  exp(tla_emtrend_50[2, 2] * x + tla_intercept_50[2, 2])}
tla_func_soybean_80 <- function(x){
  exp(tla_emtrend_80[2, 2] * x + tla_intercept_80[2, 2])}

tla_cotton_plot <- ggplot(aes(y = tla, x = fertilizer_cont, color = shade), 
                          data = subset(data, spp == 'Cotton')) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_jitter(width = 20, pch = 16, alpha = 0.5, size = 2) +
  scale_color_manual(values = c('red', 'orange', 'blue', 'purple')) +
  stat_function(fun = tla_func_cotton_0, color = 'red', lwd = 2, lty = 1) +
  stat_function(fun = tla_func_cotton_30, color = 'orange', lwd = 2, lty = 1) +
  stat_function(fun = tla_func_cotton_50, color = 'blue', lwd = 2, lty = 1) +
  stat_function(fun = tla_func_cotton_80, color = 'purple', lwd = 2, lty = 2) +
  labs(color = 'Shade (%)') +
  xlab('Fertilizer (ppm)') +
  ylab(expression('Total Leaf Area (cm' ^ '2' * ')')) +
  ggtitle(expression(italic('G. hirsutum'))) +
  labs(tag = "(a)") +
  ylim(c(0, 800))

tla_soybean_plot <- ggplot(aes(y = tla, x = fertilizer_cont, color = shade), 
                           data = subset(data, spp == 'Soybean')) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_jitter(width = 20, pch = 16, alpha = 0.5, size = 2) +
  scale_color_manual(values = c('red', 'orange', 'blue', 'purple')) +
  stat_function(fun = tla_func_soybean_0, color = 'red', lwd = 2, lty = 1) +
  stat_function(fun = tla_func_soybean_30, color = 'orange', lwd = 2, lty = 1) +
  stat_function(fun = tla_func_soybean_50, color = 'blue', lwd = 2, lty = 1) +
  stat_function(fun = tla_func_soybean_80, color = 'purple', lwd = 2, lty = 1) +
  labs(color = 'Shade (%)') +
  xlab('Fertilizer (ppm)') +
  ylab(expression('Total Leaf Area (cm' ^ '2' * ')')) +
  ggtitle(expression(italic('G. max'))) +
  labs(tag = "(b)") +
  ylim(c(0, 800))

### biomass

test(emtrends(biomass_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 0)))
test(emtrends(biomass_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 30)))
test(emtrends(biomass_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 50)))
test(emtrends(biomass_lmer, ~spp, 
              var = 'fertilizer_cont', 
              at = list(shade_cont = 80)))

biomass_emtrend_0 <- summary(emtrends(biomass_lmer, ~spp, 
                                      var = 'fertilizer_cont', 
                                      at = list(shade_cont = 0)))
biomass_emtrend_30 <- summary(emtrends(biomass_lmer, ~spp, 
                                       var = 'fertilizer_cont', 
                                       at = list(shade_cont = 30)))
biomass_emtrend_50 <- summary(emtrends(biomass_lmer, ~spp, 
                                       var = 'fertilizer_cont', 
                                       at = list(shade_cont = 50)))
biomass_emtrend_80 <- summary(emtrends(biomass_lmer, ~spp, 
                                       var = 'fertilizer_cont', 
                                       at = list(shade_cont = 80)))
biomass_intercept_0 <- summary(emmeans(biomass_lmer, ~spp, 
                                       at = list(fertilizer_cont = 0, 
                                                 shade_cont = 0)))
biomass_intercept_30 <- summary(emmeans(biomass_lmer, ~spp, 
                                        at = list(fertilizer_cont = 0, 
                                                  shade_cont = 30)))
biomass_intercept_50 <- summary(emmeans(biomass_lmer, ~spp, 
                                        at = list(fertilizer_cont = 0, 
                                                  shade_cont = 50)))
biomass_intercept_80 <- summary(emmeans(biomass_lmer, ~spp, 
                                        at = list(fertilizer_cont = 0, 
                                                  shade_cont = 80)))

biomass_func_cotton_0 <- function(x){
  exp(biomass_emtrend_0[1, 2] * x + biomass_intercept_0[1, 2])}
biomass_func_cotton_30 <- function(x){
  exp(biomass_emtrend_30[1, 2] * x + biomass_intercept_30[1, 2])}
biomass_func_cotton_50 <- function(x){
  exp(biomass_emtrend_50[1, 2] * x + biomass_intercept_50[1, 2])}
biomass_func_cotton_80 <- function(x){
  exp(biomass_emtrend_80[1, 2] * x + biomass_intercept_80[1, 2])}

biomass_func_soybean_0 <- function(x){
  exp(biomass_emtrend_0[2, 2] * x + biomass_intercept_0[2, 2])}
biomass_func_soybean_30 <- function(x){
  exp(biomass_emtrend_30[2, 2] * x + biomass_intercept_30[2, 2])}
biomass_func_soybean_50 <- function(x){
  exp(biomass_emtrend_50[2, 2] * x + biomass_intercept_50[2, 2])}
biomass_func_soybean_80 <- function(x){
  exp(biomass_emtrend_80[2, 2] * x + biomass_intercept_80[2, 2])}

biomass_cotton_plot <- ggplot(aes(y = biomass, x = fertilizer_cont, color = shade), 
                              data = subset(data, spp == 'Cotton')) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_jitter(width = 20, pch = 16, alpha = 0.5, size = 2) +
  scale_color_manual(values = c('red', 'orange', 'blue', 'purple')) +
  stat_function(fun = biomass_func_cotton_0, color = 'red', lwd = 2, lty = 1) +
  stat_function(fun = biomass_func_cotton_30, color = 'orange', lwd = 2, lty = 1) +
  stat_function(fun = biomass_func_cotton_50, color = 'blue', lwd = 2, lty = 1) +
  stat_function(fun = biomass_func_cotton_80, color = 'purple', lwd = 2, lty = 2) +
  labs(color = 'Shade (%)') +
  xlab('Fertilizer (ppm)') +
  ylab(expression('Biomass (g)'))+
  labs(tag = "(c)") +
  ylim(c(0, 6))

biomass_soybean_plot <- ggplot(aes(y = biomass, x = fertilizer_cont, color = shade), 
                               data = subset(data, spp == 'Soybean')) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_jitter(width = 20, pch = 16, alpha = 0.5, size = 2) +
  scale_color_manual(values = c('red', 'orange', 'blue', 'purple')) +
  stat_function(fun = biomass_func_soybean_0, color = 'red', lwd = 2, lty = 1) +
  stat_function(fun = biomass_func_soybean_30, color = 'orange', lwd = 2, lty = 1) +
  stat_function(fun = biomass_func_soybean_50, color = 'blue', lwd = 2, lty = 1) +
  stat_function(fun = biomass_func_soybean_80, color = 'purple', lwd = 2, lty = 1) +
  labs(color = 'Shade (%)') +
  xlab('Fertilizer (ppm)') +
  ylab(expression('Biomass (g)'))+
  labs(tag = "(d)") +
  ylim(c(0, 6))

tla_cotton_plot_g <- ggplotGrob(tla_cotton_plot)
biomass_cotton_plot_g <- ggplotGrob(biomass_cotton_plot)
tla_biomass_cotton_g <- rbind(tla_cotton_plot_g, 
                              biomass_cotton_plot_g, 
                              size = "max")
tla_soybean_plot_g <- ggplotGrob(tla_soybean_plot)
biomass_soybean_plot_g <- ggplotGrob(biomass_soybean_plot)
tla_biomass_soybean_g <- rbind(tla_soybean_plot_g,
                               biomass_soybean_plot_g, 
                               size = "max")
tla_biomass_g <- cbind(tla_biomass_cotton_g, 
                       tla_biomass_soybean_g,
                       size = 'max')

jpeg(filename = "plots/tla_biomass.jpeg", 
     width = 14, height = 11, units = 'in', res = 600)
grid.newpage()
grid.draw(tla_biomass_g)
dev.off()


## some useful calculations

### change in total photosynthetic N
propN_all_0ppm <- summary(emmeans(propN_rubisco_lmer, ~1, at = list(fertilizer_cont = 0)))[1, 2] +
  summary(emmeans(propN_bioenergetics_lmer, ~1, at = list(fertilizer_cont = 0)))[1, 2] +
  summary(emmeans(propN_lightharvesting_lmer, ~1, at = list(fertilizer_cont = 0)))[1, 2]

propN_all_630ppm <- summary(emmeans(propN_rubisco_lmer, ~1, at = list(fertilizer_cont = 630)))[1, 2] +
  summary(emmeans(propN_bioenergetics_lmer, ~1, at = list(fertilizer_cont = 630)))[1, 2] +
  summary(emmeans(propN_lightharvesting_lmer, ~1, at = list(fertilizer_cont = 630)))[1, 2]

propN_all_630ppm - propN_all_0ppm

(propN_all_630ppm - propN_all_0ppm)/propN_all_0ppm


vcmax25_shade0 <- summary(emmeans(vcmax25_lmer, ~1, at = list(shade_cont = 0)))[1, 2]
vcmax25_shade80 <- summary(emmeans(vcmax25_lmer, ~1, at = list(shade_cont = 80)))[1, 2]

(vcmax25_shade0 - vcmax25_shade80)/ vcmax25_shade80

(exp(vcmax25_shade0) - exp(vcmax25_shade80)) / 800

jmax25_shade0 <- summary(emmeans(jmax25_lmer, ~1, at = list(shade_cont = 0)))[1, 2]
jmax25_shade80 <- summary(emmeans(jmax25_lmer, ~1, at = list(shade_cont = 80)))[1, 2]

(jmax25_shade0 - jmax25_shade80)/ jmax25_shade80


vcmax25_fertilizer0 <- summary(emmeans(vcmax25_lmer, ~spp, at = list(fertilizer_cont = 0)))[1, 2]
vcmax25_fertilizer630 <- summary(emmeans(vcmax25_lmer, ~spp, at = list(fertilizer_cont = 630)))[1, 2]

(vcmax25_fertilizer630 - vcmax25_fertilizer0)/ vcmax25_fertilizer0

jmax25_fertilizer0 <- summary(emmeans(jmax25_lmer, ~spp, at = list(fertilizer_cont = 0)))[1, 2]
jmax25_fertilizer630 <- summary(emmeans(jmax25_lmer, ~spp, at = list(fertilizer_cont = 630)))[1, 2]

(jmax25_fertilizer630 - jmax25_fertilizer0)/ jmax25_fertilizer0



