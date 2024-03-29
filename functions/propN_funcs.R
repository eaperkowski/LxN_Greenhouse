###############################################################################
# p_rubisco(vcmax25, narea):
###############################################################################
#
# Calculates proportion of leaf nitrogen in Rubisco following equations from
# Niinemets et al. (1997) and Niinements et al. (1998)
#
# Function arguments:
#   - vcmax25     = maximum Rubisco carboxylation rate, standardized to 25degC
#                  (μmol m^-2 s^-1)
#   - narea       = leaf nitrogen per leaf area (gN m^-2)
#
# Returns:
# Vector with proportion of leaf N to Rubisco (p_rubisco; gN Rubisco gN^-1)
p_rubisco <- function(vcmax25, narea){
  
  vcr = 20.5 #umol CO2 (g Rubisco)-1 s-1 at 25C
  # vcmax25 in umol m-2 s-1
  # marea is g m-2
  # nmass is gN g-1
  # 6.25 converts nitrogen to protein (rubisco)
  
  p_rubisco = vcmax25 / (vcr * narea * 6.25)
  
  p_rubisco # g N Rubisco / g N
  
}

###############################################################################
# p_bioenergetics(jmax25, narea):
###############################################################################
#
# Calculates proportion of leaf nitrogen allocated to bioenergetics. Function
# follows equations from Niinemets et al. (1997) and Niinements et al. (1998)
#
# Function arguments:
#   - jmax25     = maximum RuBP regeneration rate, standardized to 25degC
#                  (μmol m^-2 s^-1)
#   - narea      = leaf nitrogen per leaf area (gN m^-2)
#
#
# Returns:
# Vector with proportion of leaf N to bioenergetics (p_bioenergetics; 
# g N bioenergetics g N^-1)
p_bioenergetics <- function(jmax25, narea){
  
  jmc = 156 # capacity of electron transport per unit of cytochrome f (umol electrons (umol cyt f)-1 s-1)
  # jmax25 in umol m-2 s-1
  # marea is g m-2
  # nmass is gN g-1
  # 8.06 converts nitrogen to protein (cytf)
  
  p_bioenergetics = jmax25 / (jmc * narea * 8.06)
  
  p_bioenergetics # g N bioenergetics / g N
  
}

###############################################################################
# p_lightharvesting(chlorophyll, nmass):
###############################################################################
#
# Calculates proportion of leaf nitrogen allocated to light harvesting. Function
# follows equations from Niinemets et al. (1997) and Niinements et al. (1998)
#
# Function arguments:
#   - chlorophyll = chlorophyll content (mmol g^-1)
#   - nmass       = leaf nitrogen per leaf mass (g N g^-1)
#
#
# Returns:
# Vector with proportion of leaf N to light harvesting (p_lightharvesting; 
# g N light harvesting g N^-1)
p_lightharvesting <- function(chlorophyll, nmass){
  
  # chlorophyll content in mmol g-1
  # nmass in gN g-1
  # marea in g m-2
  chlorophyll_mol = chlorophyll / 1000
  cb = 2.75 / 1000 # chlorophyll binding in mol chlorophyll (g N chlorophyll)-1 assuming most is in LHCII
  
  p_lightharvesting = chlorophyll_mol *  (1 / nmass) * (1 / cb)
  
  p_lightharvesting # g N light harvesting / g N
  
  # cytf_m2 = (jmax25 / 156) / 1000 # mmol cyt f m-2
  # cytf_chlor = cytf_m2 * (1/marea) * (1/chlorophyll_mol) # mmol cyt f (mol chlorophyll)-1
  # psii_chlor = 1.98 * cytf_chlor - (0.365 * (cytf_chlor^2)) # mmol (mol chlorophyll)-1
  # psi_chlor = 1.7 # mmol (mol chlorophyll)-1
  # psii_g = psii_chlor * (chlorophyll / 1000) # mmol g-1
  # psi_g = psi_chlor * (chlorophyll / 1000) # mmol g-1
  # 
  # lhcii_g = chlorophyll - (psii_g + psi_g)
  # 
  # psii_prop = psii_g / chlorophyll
  # psi_prop = psi_g / chlorophyll
  # lhcii_prop = lhcii_g / chlorophyll
  # 
  # cb_psii = 0.858 # mmol chl (gN)-1
  # cb_psi = 2.18 # mmol chl (gN)-1
  # cb_lhcii = 2.75 # mmol chl (gN)-1
  # 
  # p_lightharvesting = (chlorophyll / nmass) * (1/((psii_prop*cb_psii) + (psi_prop*cb_psi) +(lhcii_prop*cb_lhcii)))
}

###############################################################################
# References
###############################################################################

# Niinemets Ü., Kull O. & Tenhunen J.D. (1998) An analysis of light effects on 
# foliar morphology, physiology, and light interception in temperate deciduous 
# woody species of contrasting shade tolerance. Tree Physiology 18, 681–696.

# Niinemets Ü. & Tenhunen J.D. (1997) A model separating leaf structural and 
# physiological effects on carbon gain along light gradients for the 
# shade-tolerant species Acer saccharum. Plant, Cell & Environment 20, 845–866.
