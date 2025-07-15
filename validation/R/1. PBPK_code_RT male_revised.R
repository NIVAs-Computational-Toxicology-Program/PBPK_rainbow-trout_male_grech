### Model PBPK (rainbow trout / zebrafish / fathead minnows / stickleback) generic ###
# 10/2017
# A. Grech, C.Tebby, C. Brochot, R. Beaudouin.

## clean environment
rm(list=ls())


##############################
### Variable -(Notation)- Units
# Quantity        -(Q_x) - microg 
# Volumes:      -(V_x) - mL
# Time:            -(t)   - d
# Flows:           -(F_x) - mL/d
# Concentrations: -(C_x) - microg/mL
# Masses:         -(BW)  - g
# Lenght:          -(L)   - mm
# Temperature: -(TC_c)  - Celsius
#     ''                 -(TC_k)  - Kelvin
# Density of each tissue is considered equal to 1
###############################



############
## Libraries
library("deSolve") # for ODE solving


# Arrhenius temperatures function 
KT <- function(T, TR, TA){ exp ( (TA / TR) - (TA / T) ) } # Implements the Arrhenius equation to scale physiological parameters by temperature (in Kelvin).


# PBPK function -  Defines the PBPK model ODE system. Inputs: t: time, Initial.Values: initial values of state variables (quantities in compartments), Parameters: model parameters, input.T: temperature function over time
PBPK.model <- function(t, Initial.Values, Parameters, input.T) { 
  
  with( as.list(c(Initial.Values, Parameters)),{
    
    # Temperature (Converts input temperature to Kelvin and Celsius)
    TC_k = input.T(t) + 273.15 # (degree K) 
    TC_c = input.T(t)          # (degree C) #
    
    # body weight : DEB growth model anisomorphic (Calculates temperature-scaled growth rate (dL) and resulting body weight (BW) using Dynamic Energy Budget (DEB) theory.)
    DEB_v_t   = DEB_v  * KT(T=TC_k , TR=TR_DEB ,TA=TA)  # mm/d 
    DEB_Lm    = DEB_v / (DEB_KM * DEB_g)                # mm  
    DEB_M     = (DEB_EHm / DEB_EHb)^(1/3)               # EHm and EHb = J
    dL        = (DEB_v_t / (3 * (f_cst + DEB_g))) * (f_cst * DEB_M - (L/DEB_Lm)) # Length increment based on DEB growth
    BW        = a_BW_L *( (L/10)/DEB_shape ) ^ b_BW_L   # BW = g; a = g/cm;  L = mm --> /10 = cm
    
    #  volumes (mL)- Defines volumes of each compartment as fractions of body weight.
    V_art    = sc_blood  * BW  * frac_art_ven
    V_ven	 = sc_blood  * BW  - V_art
    V_liver  = sc_liver  * BW
    V_gonads = sc_gonads * BW
    V_brain  = sc_brain  * BW
    V_fat    = sc_fat    * BW
    V_skin   = sc_skin   * BW
    V_GIT    = sc_GIT    * BW
    V_kidney = sc_kidney * BW
    V_rp	 = sc_rp     * BW
    V_pp     = BW - (V_art + V_ven + V_liver + V_gonads + V_brain + V_fat + V_skin + V_GIT + V_kidney + V_rp)
    
    V_bal    = BW - (V_art + V_ven + V_liver + V_gonads + V_brain + V_fat + V_skin + V_GIT + V_kidney + V_rp + V_pp) # balance
    
    # Blood flow (mL/d) - Calculates organ-specific blood flows using allometric scaling.
    F_card_g = F_card_ref * KT(T=TC_k , TR=TR_Fcard,TA=TA)  * (BW/Bw_Fcard_ref)^(-0.1) # cardiac output = mL/d/g
    F_card   = F_card_g * BW # mL/d
    
    F_liver  = frac_liver  * F_card
    F_gonads = frac_gonads * F_card
    F_brain  = frac_brain  * F_card
    F_fat    = frac_fat    * F_card
    F_skin   = frac_skin   * F_card
    F_GIT    = frac_GIT    * F_card
    F_kidney = frac_kidney * F_card
    F_rp	   = frac_rp     * F_card
    F_pp     = F_card - (F_liver + F_gonads + F_brain + F_fat + F_skin + F_GIT + F_kidney + F_rp)
    
    F_bal = F_card - (F_liver + F_gonads + F_brain + F_fat + F_skin + F_GIT + F_kidney + F_rp + F_pp) # balance
    
    # Effective respiratory volume (mL/d)- Computes water intake rate based on oxygen demand and respiratory uptake rate.
    V_O2_g = V_O2_ref * KT(T=TC_k, TR=TR_VO2,TA=TA)  * (BW/Bw_VO2_ref)^(-0.1) #mg O2/d/g 
    V_O2   = V_O2_g * BW   # mg O2/d
    
    C_O2_water = ((-0.24 * TC_c + 14.04) * Sat )/10^3    # mg O2/mL
    F_water    = V_O2/ (O2_EE * C_O2_water)              # mL/d
    
    Kx         = min (F_water, F_card * PC_blood_water)  # mL/d
    
    # Concentrations in tissues (microg/g = microg/mL)- Calculates chemical concentrations in blood and tissues.
    C_art    = Q_art    / V_art
    C_plasma = C_art    / Ratio_blood_plasma
    C_ven    = Q_ven    / V_ven
    C_liver  = Q_liver  / V_liver
    C_gonads = Q_gonads / V_gonads
    C_brain  = Q_brain  / V_brain
    C_fat    = Q_fat    / V_fat
    C_skin   = Q_skin   / V_skin
    C_GIT    = Q_GIT    / V_GIT
    C_kidney = Q_kidney / V_kidney
    C_rp     = Q_rp     / V_rp
    C_pp     = Q_pp     / V_pp
    C_tot    = ((Q_art + Q_ven + Q_liver + Q_gonads + Q_brain + Q_fat + Q_skin + Q_GIT + Q_kidney + Q_pp + Q_rp)
                / (V_art + V_ven + V_liver + V_gonads + V_brain + V_fat + V_skin + V_GIT + V_kidney + V_pp + V_rp))
    
    # concentration in the water : exposure --> microg/g = microg/mL - Models exposure concentration and its elimination from water.
    C_water       = Q_water / V_water
    dQ_elim_water = Ke_water * Q_water
    
    # concentration in urine - Urine volume and concentration calculation
    dV_urine        = urine_rate * BW
    C_urine         = Q_urine  / V_urine
    
    # Scaling clearence and excretion constant - Scale metabolic clearance to tissue volume or body weight
    if ( !is.na(Vmax)) { Vmax_sc        =  Vmax * V_liver  } 
    if ( !is.na(Cl_liver)) { Cl_liver_sc    = Cl_liver * V_liver }
    rate_plasma_sc = rate_plasma * BW
    
    # Excretion - Calculate loss of chemical via feces, gills, bile, urine.
    dQ_feces        = Q_lumen_GIT * Ke_feces
    dQ_excret_gills = Kx * (C_ven * Unbound_fraction / PC_blood_water)
    dQ_bile         = (Ke_bile * Q_liver) - (Q_bile * K_BG)
    dQ_urine        = Ke_urine * Q_kidney               # initialized by scenario for urine discharge
    dQ_urine_cum    = Ke_urine * Q_kidney               # not initialized 
    dQ_excret       = dQ_excret_gills + dQ_urine_cum + dQ_feces
    
    
    # Quantity absorbed : differentials in microg/d - Models chemical uptake from water and gastrointestinal tract
    dQ_admin_gills = Kx * C_water         ###Gills
    dQ_lumen_GIT   = ( - Ku * Q_lumen_GIT
                       - Q_lumen_GIT * Ke_feces
                       + Q_bile * K_BG)  ###GIT
    dQ_admin_GIT   = 0
    
    # Metabolism - First-order or Michaelis-Menten metabolism in plasma and liver
    dQ_met_plasma  = rate_plasma_sc * Q_ven
    
    dQ_met_liver    = ifelse ( is.na(Cl_liver),
                               ( ( Vmax_sc *(C_liver/PC_liver) )/(Km + C_liver/PC_liver) ), 
                               ( Cl_liver_sc * C_liver / PC_liver) )
    
    dQ_met = dQ_met_plasma + dQ_met_liver
    
    
    # Blood quantity - Mass balance in arterial and venous compartments
    dQ_art    = ( F_card * C_ven - F_fat * C_art - F_rp * C_art - F_pp * C_art
                  - F_liver  * C_art - F_kidney * C_art   - F_GIT * C_art
                  - F_gonads * C_art  - F_skin   * C_art - F_brain  * C_art )
    
    dQ_ven    = (F_brain * C_brain/PC_brain
                 + (F_liver + F_rp + F_gonads + F_GIT ) * C_liver/PC_liver
                 +  F_fat   * C_fat/PC_fat
                 + (F_kidney + ((1 - a_Fpp) * F_pp) + ((1 - a_Fs) * F_skin)) * (C_kidney / PC_kidney)
                 + a_Fpp * F_pp  * (C_pp / PC_pp) 
                 + a_Fs * F_skin * (C_skin / PC_skin)
                 - F_card * C_ven  
                 + dQ_admin_gills - dQ_excret_gills
                 - dQ_met_plasma)
    
    # Quantity in tissues - Net exchange of chemical between blood and each tissue
    dQ_gonads = F_gonads*(C_art - C_gonads/PC_gonads)
    dQ_fat    = F_fat   *(C_art - C_fat/PC_fat)
    dQ_rp     = F_rp    *(C_art - C_rp/PC_rp)
    dQ_pp     = F_pp    *(C_art - C_pp/PC_pp)
    dQ_brain  = F_brain *(C_art - C_brain/PC_brain)
    dQ_skin   = F_skin  *(C_art - C_skin/PC_skin)
    
    dQ_liver  =  ( F_liver * C_art  
                   + F_rp *(C_rp/PC_rp)
                   + F_GIT *(C_GIT/PC_GIT)
                   + F_gonads*(C_gonads/PC_gonads)
                   - ( F_liver + F_rp + F_GIT + F_gonads  ) *(C_liver/PC_liver)  
                   - (Ke_bile * Q_liver) - dQ_met_liver )
    
    dQ_kidney = ( F_kidney * C_art 
                  + (1-a_Fpp) * F_pp * (C_pp / PC_pp) 
                  + (1-a_Fs) * F_skin * (C_skin / PC_skin) 
                  - (F_kidney + (1-a_Fpp) * F_pp + (1-a_Fs) * F_skin) * (C_kidney / PC_kidney) 
                  - dQ_urine)
    
    dQ_GIT    = ( Ku * Q_lumen_GIT +  F_GIT*(C_art - C_GIT/PC_GIT))
    
    # Chemical kinetic in  aquarium water - Net change of chemical in aquarium water: Models the chemical concentration in water, accounting for administration, elimination, and excretion.
    dQ_water = (dQ_excret - dQ_elim_water - dQ_admin_gills)
    
    
    # Mass-balance (Ensures conservation of mass and tracks system balance). List of:Derivatives (dQ_x), Concentrations (C_x), Volumes (V_x), Mass balances, Scalars (BW, F_card, Kx, etc.)
    Q_Body      = (Q_art + Q_ven + Q_bile + Q_liver + Q_gonads + Q_brain + Q_fat + Q_skin + Q_kidney + Q_GIT + Q_pp + Q_rp) 
    Q_admin_tot = Q_admin_gills + ivQuantity - Q_lumen_GIT + Q_admin_GIT # amount entering body
    Q_elim_tot  = Q_met + Q_excret
    Mass_Bal    = Q_admin_tot - Q_Body - Q_elim_tot
    Mass_Bal_Sys = (Q_admin_tot - Q_admin_gills) - Q_Body - ( Q_elim_tot - Q_excret) - Q_water - Q_elim_water																										 
    
    
    list("Q"=c(dQ_water, dL, dV_urine,
               dQ_art,  dQ_ven,  
               dQ_brain, dQ_fat,  dQ_GIT,  dQ_lumen_GIT, dQ_gonads, 
               dQ_kidney, dQ_liver, dQ_skin, dQ_rp, dQ_pp,
               dQ_admin_gills, dQ_admin_GIT,
               dQ_met_plasma, dQ_met_liver, dQ_met, dQ_elim_water,
               dQ_excret, dQ_feces, dQ_excret_gills, dQ_bile, dQ_urine, dQ_urine_cum), # derivatives
         
         
         c( C_water = C_water,
            C_art = C_art,       C_ven = C_ven,      C_plasma = C_plasma,
            C_brain = C_brain,   C_fat = C_fat,      C_GIT = C_GIT, C_gonads = C_gonads,
            C_kidney = C_kidney, C_liver = C_liver,  C_skin = C_skin, 
            C_rp    = C_rp,      C_pp = C_pp, C_urine = C_urine, C_tot = C_tot),
         
         c( V_water = V_water,
            V_art = V_art,  V_ven = V_ven, 
            V_brain = V_brain, V_fat = V_fat, V_GIT = V_GIT,  V_gonads= V_gonads, 
            V_kidney = V_kidney, V_liver = V_liver, V_skin = V_skin, 
            V_rp  = V_rp,  V_pp = V_pp),
         
         c(V_bal = V_bal, 
           F_bal = F_bal, 
           Q_admin_tot = Q_admin_tot, 
           Q_elim_tot  = Q_elim_tot, 
           Q_Body = Q_Body, 
           Mass_Bal = Mass_Bal,
           Mass_Bal_Sys =Mass_Bal_Sys),	   
         
         
         "BW"      = BW,
         "F_card"  = F_card,
         "F_water" = F_water,
         "Kx"      = Kx)
  }) }


#################################################
### Partition coefficient function: QSAR model ##
### Partition coefficient function #####
# 10/2017
# A. Grech, C.Tebby, C. Brochot, R. Beaudouin.
##############################

##############################
### Units:
# Volumes:          mL
# Masses:           g
###############################

# Estimate tissue:blood partition coefficients using QSARs from:a) Nichols et al. (2006) for high logKow, b) Bertelsen et al. (1998) for lower logKow
PC_qsar_model = function(a_PC  = 0.78,  # PC between blood and water and logKow (Nichols et al. 2006)
                         b_PC  = 0.82,  # Empirical intercept for blood:water PC vs. logKow (Nichols et al. 2006)
                         
                         e_Bar = 0.05, # Parameter 1 of Nichols et al. 2006
                         
                         a_Bar = 0.74, # Parameter 1 of Bertelsen et 1998
                         b_Bar = 0.72, # Parameter 2 of Bertelsen et 1998
                         c_Bar = 1.00, # Parameter 3 of Bertelsen et 1998
                         
                         log_Kow, # parameter from xlogKOW for BaP (accounts for non-lipid organic matter)
                         
                         # Water content of tissues (fractions of tissue volume)
                         water_liver,  water_gonads, water_brain,
                         water_fat,    water_skin,   water_GIT,
                         water_kidney, water_rp,     water_pp,
                         
                         # Lipid content of tissues (fractions of tissue volume)
                         lipids_liver, lipids_gonads, lipids_brain,  lipids_fat,
                         lipids_skin,  lipids_GIT,    lipids_kidney, lipids_rp,
                         lipids_pp,
                         
                         # Optional user-supplied PC overrides (if not NULL, these override the QSAR estimates)
                         PC_bw, PC_l, PC_go, PC_b, PC_f, PC_s, PC_gi, PC_k, PC_r, PC_p   ){
  
  
  
  #===============================================
  # (I.) QSAR Model : partition coefficient between blood and water, using Nichols et al. regression
  #=============================================== 
  
  PC_QSAR_blood_water = 10^(a_PC * log_Kow - b_PC)
  
  PC_blood_water = ifelse (is.null(PC_bw), PC_QSAR_blood_water, PC_bw)  
  
  if( log_Kow > 4.04) {
    
    #===============================================
    # (I.) QSAR Model
    # Nichols, J W.; Schultz, I R.; Fitzsimmons, P N.
    # Aquat. Toxicol. 2006, 78, 74-90.
    #=============================================== 
    
    #PC = lipidcontent * Kow + 0.05 * nonlipidcontent * Kow + watercontent
    
    # Calcul des PCs # SU 
    PC_QSAR_liver       = (water_liver   +  lipids_liver  * 10^( log_Kow) +  e_Bar * (1 - lipids_liver  - water_liver ) )/PC_blood_water
    PC_QSAR_gonads      = (water_gonads  +  lipids_gonads * 10^( log_Kow) +  e_Bar * (1 - lipids_gonads - water_gonads) )/PC_blood_water
    PC_QSAR_brain       = (water_brain   +  lipids_brain  * 10^( log_Kow) +  e_Bar * (1 - lipids_brain  - water_brain)  )/PC_blood_water
    PC_QSAR_fat         = (water_fat     +  lipids_fat    * 10^( log_Kow) +  e_Bar * (1 - lipids_fat    - water_fat )   )/PC_blood_water
    PC_QSAR_skin        = (water_skin    +  lipids_skin   * 10^( log_Kow) +  e_Bar * (1 - lipids_skin   - water_skin )  )/PC_blood_water
    PC_QSAR_GIT         = (water_GIT     +  lipids_GIT    * 10^( log_Kow) +  e_Bar * (1 - lipids_GIT    - water_GIT )   )/PC_blood_water
    PC_QSAR_kidney      = (water_kidney  +  lipids_kidney * 10^( log_Kow) +  e_Bar * (1 - lipids_kidney - water_kidney) )/PC_blood_water
    PC_QSAR_rp          = (water_rp      +  lipids_rp     * 10^( log_Kow) +  e_Bar * (1 - lipids_rp     - water_rp )    )/PC_blood_water
    PC_QSAR_pp          = (water_pp      +  lipids_pp     * 10^( log_Kow) +  e_Bar * (1 - lipids_pp     - water_pp )    )/PC_blood_water
    
  } else { 
    
    #===============================================
    # (II.) QSAR Model
    # Bertelsen, S. L.; Hoffman, A. D.; Gallinat, C. A.; Elonen, C. M.; Nichols, J. W.
    # Environ. Toxicol. Chem. 1998, 17, 1447-1455.
    #=============================================== 
    
    # for lower hydrophobicity (logKow ≤ 4.04), use Bertelsen et al. (1998)
    
    PC_QSAR_liver       = (water_liver   + 10^(a_Bar*log_Kow + b_Bar + c_Bar * log10(lipids_liver )  )) /PC_blood_water
    PC_QSAR_gonads      = (water_gonads  + 10^(a_Bar*log_Kow + b_Bar + c_Bar * log10(lipids_gonads)  )) /PC_blood_water
    PC_QSAR_brain       = (water_brain   + 10^(a_Bar*log_Kow + b_Bar + c_Bar * log10(lipids_brain)   )) /PC_blood_water
    PC_QSAR_fat         = (water_fat     + 10^(a_Bar*log_Kow + b_Bar + c_Bar * log10(lipids_fat)     )) /PC_blood_water
    PC_QSAR_skin        = (water_skin    + 10^(a_Bar*log_Kow + b_Bar + c_Bar * log10(lipids_skin)    )) /PC_blood_water
    PC_QSAR_GIT         = (water_GIT     + 10^(a_Bar*log_Kow + b_Bar + c_Bar * log10(lipids_GIT)     )) /PC_blood_water
    PC_QSAR_kidney      = (water_kidney  + 10^(a_Bar*log_Kow + b_Bar + c_Bar * log10(lipids_kidney)  )) /PC_blood_water
    PC_QSAR_rp          = (water_rp      + 10^(a_Bar*log_Kow + b_Bar + c_Bar * log10(lipids_rp)      )) /PC_blood_water
    PC_QSAR_pp          = (water_pp      + 10^(a_Bar*log_Kow + b_Bar + c_Bar * log10(lipids_pp)      )) /PC_blood_water
  }
  
  #===============================================
  
  # (III.) Override QSAR partition coefficients with user-supplied values if provided
  PC_liver   = ifelse (is.null(PC_l), PC_QSAR_liver, PC_l)
  PC_gonads  = ifelse (is.null(PC_go), PC_QSAR_gonads, PC_go)  
  PC_brain   = ifelse (is.null(PC_b), PC_QSAR_brain, PC_b)
  PC_fat     = ifelse (is.null(PC_f), PC_QSAR_fat, PC_f) 
  PC_skin    = ifelse (is.null(PC_s), PC_QSAR_skin, PC_s)
  PC_GIT     = ifelse (is.null(PC_gi), PC_QSAR_GIT, PC_gi)   
  PC_kidney  = ifelse (is.null(PC_k), PC_QSAR_kidney, PC_k) 
  PC_rp      = ifelse (is.null(PC_r), PC_QSAR_rp, PC_r)
  PC_pp      = ifelse (is.null(PC_p), PC_QSAR_pp, PC_p) 
  
  # return a named vector of final tissue:blood partition coefficients
  return(   c("PC_blood_water" = PC_blood_water,  "PC_liver"   = PC_liver,  "PC_gonads"      = PC_gonads,
              "PC_brain"       = PC_brain,        "PC_fat"     = PC_fat ,   "PC_skin"        = PC_skin,
              "PC_GIT"         = PC_GIT,    "PC_kidney"      = PC_kidney,
              "PC_rp"          = PC_rp,       "PC_pp"          = PC_pp ))
}

############################################################
####  Parameterisation for PBTK model male rainbow trout  ###
### Parameterisation for PBPK model male/female rainbow trout #####
## Note: main differences are indicated in the composition of tissues
# 10/2017
# A. Grech, C.Tebby, C. Brochot, R. Beaudouin.
# 
##############################

##############################
### Units:
###
# Quantity          microg
# Volumes:          mL
# Time:             d
# Flows:            mL/d
# Concentrations:   microg/mL
# Vmax:             microg/d/mL Liver
# Km:               microg/mL
# Masses:           g
# Lenght:           mm
# Temperature:      Celsius
# Ventilation rate: mL/d

# Density of each tissue is considered equal to 1

###############################

library("deSolve")


#===============================================
# source("PBTK_model.r")
# source("PBTK_PC_QSAR_function.r")
#===============================================

#===============================================
# (II.) Model integration
#===============================================

PBPK.RT.M = function(
    #===============================================
    # (II.1) Parameters
    #===============================================
    
    # Physiological parameters (species- and life-stage-specific constants, derived from experiments, databases (find my pet) or literature (e.g., Elliott 1969, Curtis 1991).)
    Bw_i       = NULL    ,  # Body weight initial (g fresh weight) - need to be filled in
    Bw_Fcard_ref= 270.1  ,  # Body weight of reference for F_card (from Elliott (1969)?)
    Bw_VO2_ref = 1000    ,  # Body weight of reference for VO2 from Elliot experiments 
    DEB_v      = 0.38    ,  # Energy conductance (mm/d) (DEB model parameter) - from https://www.bio.vu.nl/thb/deb/deblab/add_my_pet/index.html  
    DEB_g      = 0.196   ,  # Energy investment ratio (SU) (DEB model parameter) - from https://www.bio.vu.nl/thb/deb/deblab/add_my_pet/index.html
    DEB_KM     = 0.027   ,  # Somatic maintenance rate coefficient (1/d) (DEB model parameter) - from https://www.bio.vu.nl/thb/deb/deblab/add_my_pet/index.html
    DEB_EHm    = 452     ,  # Energy at State of maturity at metamorphosis (J) - from https://www.bio.vu.nl/thb/deb/deblab/add_my_pet/index.html
    DEB_EHb    = 52      ,  # Energy at State of maturity at birth (J) - from https://www.bio.vu.nl/thb/deb/deblab/add_my_pet/index.html
    DEB_shape  = 0.113   ,  # Shape factor (dimensionless) for surface/volume scaling - from https://www.bio.vu.nl/thb/deb/deblab/add_my_pet/index.html
    a_BW_L     = 0.0051  ,  # Allometric coefficient to convert length (cm) to weight (g) - a relation BW(g)=F(L(cm)) 
    b_BW_L     = 3.3513  ,  # # Allometric exponent for length-weight relationship - b relation BW(g)=F(L(cm))
    
    # Environmental condition
    Temperature  = 10   ,    # ambient water temperature (celsius)
    TA           = 6930 ,    # Arrhenius temperature  in Kelvin, for scaling rates
    TR_DEB       = 293.65,   # Reference temperature for DEB scaling (K)
    TR_Fcard     = 279.15,   # Reference T for cardiac output (K)
    TR_VO2       = 283.15,   # Reference T for VO2 scaling (K)
    f_cst        = 1    ,    # food level 1 = ad-libitum, 0= starvation
    V_water      = 1000 ,    # Volume of exposure water (aquarium) (mL)
    
    # Effective respiratory volume & cardiac output  
    F_card_ref  = 28.5,  # reference cardiac output (mL/min/g): F_card_ref = mL/min/g --> F_card_ref * (60*24) = mL/d/g
    V_O2_ref    = 3.26, # reference oxygen comsuption rate (mg O2/d/g), determined by visual fitting from Elliot, 1969 experiments
    O2_EE       = 0.71, # Oxygen extraction efficiency of (fraction of extracted water) - default: 71% proposed by Erickson, 1990
    Sat         = 0.90, # dissolved oxygen saturation in water - default: 90% proposed by Erickson, 1990
    frac_art_ven = 1/3, # fraction of arterial blood vs. veinous blood - default: 1/3 proposed by XXX.
    
    sc_blood  = 0.0527, # volume scaling factor: fraction of BW (%)
    sc_gonads = 0.018 , # volume scaling factor: fraction of BW (%), default: male = 0.018, female: 0.199
    sc_brain  = 0.0049, # volume scaling factor: fraction of BW (%)
    sc_liver  = 0.0146, # volume scaling factor: fraction of BW (%), male: 0.0146, female:0.0165 
    sc_fat    = 0.0797, # volume scaling factor: fraction of BW (%)
    sc_skin   = 0.10  , # volume scaling factor: fraction of BW (%)
    sc_GIT    = 0.099 , # volume scaling factor: fraction of BW (%)
    sc_kidney = 0.0076, # volume scaling factor: fraction of BW (%)
    sc_rp     = 0.0257, # volume scaling factor: fraction of BW (%), male: 0.0257, female: 0.0189 
    
    sc_pp     = (1 -  sc_blood - sc_gonads - sc_brain - sc_liver - sc_fat  # Volue sclaing factor, poorly perfused tissue (residual)
                 -  sc_skin - sc_GIT - sc_kidney -sc_rp),
    
    frac_gonads =0.0138  , # Fraction of arterial blood flow
    frac_brain  =0.018   , # Fraction of arterial blood flow
    frac_liver  =0.0158  , # Fraction of arterial blood flow
    frac_fat    =0.0059  , # Fraction of arterial blood flow
    frac_skin   =0.0570  , # Fraction of arterial blood flow
    frac_GIT    =0.1759  , # Fraction of arterial blood flow
    frac_kidney =0.0817  , # Fraction of arterial blood flow
    frac_rp     =0.0952  , # Fraction of arterial blood flow
    
    frac_pp     = (1 - frac_gonads - frac_brain - frac_liver - frac_fat # Remaining flow to poorly perfused tissues
                   - frac_skin - frac_GIT -  frac_kidney - frac_rp),
    
    # Exposure quantity (microg)
    WaterQuantity  = NULL, # Amount of chemical introduced via water - need to be filled in.
    IngestQuantity = NULL, # Amount ingested orally (feed)- need to be filled in.
    ivQuantity     = NULL, # Intravenous dosing (µg)- need to be filled in.
    
    # Chemical parameters 
    log_Kow        = NULL , # Octanol-water partition coefficient (used for PC estimation) - need to be filled in
    Unbound_fraction = 1, # Fraction of unbound chemical in blood (bioavailable). between 0 and 1. assumption: Limits excretion by gills.
    Ratio_blood_plasma = 1, # blood/plasma global equilibrium concentration ratio (accounts for hemtocrite, fu, ...)
    
    # Organs composition (Relative to total volume, % of tissue volume)
    water_liver = 0.746   , # relative value (%)
    water_brain = 0.75    , # male: 0.75, female: 0.76 
    water_gonads= 0.66    , # female value
    water_fat   = 0.05    ,
    water_skin  = 0.667   ,
    water_GIT   = 0.582   , 
    water_kidney= 0.789   ,
    water_rp    = 0.5293  , # value of zebrafish model
    water_pp    = 0.769   , # muscle value
    
    
    lipids_liver = 0.045 , # relative value (%)
    lipids_brain = 0.073 , # male: 0.073, female: 0.076 
    lipids_gonads= 0.0563 , # female value 
    lipids_fat   = 0.942 , 
    lipids_skin  = 0.029 ,
    lipids_GIT   = 0.288 ,
    lipids_kidney= 0.052 ,
    lipids_rp    = 0.0677 , # value of zebrafish model
    lipids_pp    = 0.03  , # muscle value
    
    # PARTITION COEFFICIENT OVERRIDES (optional manual values)
    PC_bw = NULL,  # Blood:water PC
    PC_l  = NULL,  # Liver:blood
    PC_go = NULL,  # Gonads:blood
    PC_b  = NULL,  # Brain:blood
    PC_f  = NULL,  # Fat:blood
    PC_s  = NULL,  # Skin:blood
    PC_gi = NULL,  # GIT:blood
    PC_k  = NULL,  # Kidney:blood
    PC_r  = NULL,  # Richly perfused:blood
    PC_p  = NULL,  # Poorly perfused:blood
    
    # VENOUS RETURN FRACTIONS
    a_Fpp = 0.4, # Fraction of PPT blood going to venous
    a_Fs  = 0.1, # Fraction of skin blood going to venous
    
    # Oral absorption  parameters
    Ku       = 0,   # Uptake rate constant/diffusion coefficient (1/d) from lumen to GIT wall 
    frac_abs = 0,   # Absorption fraction.  /!\ Parametrisation with frac_abs implies Ke_feces null ! # Total absorbed fraction from gut (bypasses Ku if defined). # Note: if frac_abs is used, Ke_feces should be set to 0 to avoid double-counting losses
    
    # Metabolism
    Km   = NULL , # Michaelis constant (microg/ml)
    Vmax = NA, # Max metabolic rate (microg/d/mL liver)
    Cl_liver   = NA , # Liver clearance (mL/d/mL liver) — alternative to Vmax/Km (mL/d/mL liver)
    rate_plasma  = 0, # Plasma metabolism rate (mg/d/g fish)
    
    # Excretion
    Ke_bile  = 0,     # First-order bile excretion rate (1/d) - need to be fille din
    Ke_urine = 0,     # Urine elimination rate (1/d)- need to be fille din
    Ke_water = 0,     # Water elimination rate (gills) (1/d)- need to be fille din
    K_BG     = 0,     # Transfer rate from bile to gut lumen (1/d). If 0 → storage; if large → fast elimination. 0 : no bile excretion: accumulation in bilary vesicule # 1E12 : Complete bile excretion in GIT_lumen
    
    # Feces and urination
    Ke_feces   = 0.83 , # Fecal excretion rate (1/d) estiamted from Nichols  et al. 2004
    urine_rate = 1.2e-03*60*24/29.82, # # Daily urine rate (mL/g/d), derived from burst frequency and volumeV_burst = 1.2 mL.kg-1 every 29.82 minutes proposed by Curtis 1991 --> 1.2e-03 mL.g BW-1
    urination_interval = NULL  ,    # Optional fixed interval between urination events. Default: 30min proposed by Curtis 1991 --> 30/60/24hr
    
    # temps
    Times  = NULL,       # number of days simulation (days) - neeed to be filled in
    Time_end_exposure = max(Times), # Last day of exposure days- neeed to be filled in
    period = NULL,       # days between two repeated doses- neeed to be filled in
    frac_renewed    = NULL, # fraction of the water of the aquaria renewed- neeed to be filled in
    time_final_dose = NULL, # Time of last administered dose- neeed to be filled in
    time_first_dose = NULL # Time of first administered dose- neeed to be filled in
){
  
  
  #===============================================
  # (I.) Vectors of parameters
  #=============================================== 
  
  # Physiological parameters
  parms_physio <- c (# DEB model parameters
    "Bw_i" = Bw_i,                     # Initial body weight of the fish (g)
    "Bw_VO2_ref" = Bw_VO2_ref,         # Reference BW for scaling VO2 (oxygen consumption)
    "Bw_Fcard_ref" = Bw_Fcard_ref,     # Reference BW for cardiac output scaling
    "DEB_v" = DEB_v,                   # Energy conductance (growth rate driver)
    "DEB_g" = DEB_g,                   # Energy investment ratio
    "DEB_KM" = DEB_KM,                 # Somatic maintenance rate (1/d)
    "DEB_EHm" = DEB_EHm,               # Maturity at metamorphosis (J)
    "DEB_EHb" = DEB_EHb,               # Maturity at birth (J)
    "DEB_shape"  = DEB_shape,          # Shape factor (scales surface:volume)
    "a_BW_L" = a_BW_L, "b_BW_L" = b_BW_L, # Length-weight scaling parameters
    
    # cardiac ouptut & respiratory parameters
    "F_card_ref" = F_card_ref,         # Reference cardiac output (mL/min/g)
    "O2_EE" = O2_EE,                   # Oxygen extraction efficiency (fraction)
    "Sat" = Sat,                       # Dissolved oxygen saturation in water
    "V_O2_ref" = V_O2_ref,             # Reference VO2 (mg O2/d/g)
    "frac_art_ven" = frac_art_ven,     # Arterial to venous blood volume ratio
    
    # relative flux to cardiac output (lood Flow Fractions (as % of cardiac output)) 
    "frac_liver" = frac_liver, # Liver - blood flow fraction
    "frac_gonads"= frac_gonads, # gonads - blood flow fraction
    "frac_brain"  = frac_brain, # brain - blood flow fraction
    "frac_fat"   = frac_fat ,   # fat - blood flow fraction
    "frac_skin"   = frac_skin, # skin - blood flow fraction
    "frac_GIT" = frac_GIT , # GIT - blood flow fraction
    "frac_kidney" = frac_kidney, # kidney - blood flow fraction
    "frac_rp"    = frac_rp, # Richly perfused tissue- blood flow fraction
    "frac_pp"    = frac_pp, # Poorly perfused tissue- blood flow fraction
    
    # relative weight to BW. Tissue Volume Fractions (as % of BW)
    "sc_GIT"    = sc_GIT,    
    "sc_skin"  = sc_skin,   
    "sc_blood"   = sc_blood,  
    "sc_liver"  = sc_liver,  
    "sc_gonads"= sc_gonads, 
    "sc_fat"   = sc_fat,
    "sc_brain"  = sc_brain,  
    "sc_rp"    = sc_rp, 
    "sc_pp"    = sc_pp,     
    "sc_kidney" = sc_kidney,
    
    # urination
    "urine_rate" = urine_rate, # Urine excretion rate (mL/g/d)
    "urination_interval" = urination_interval # Time between urinations (if pulsed excretion)
  )
  
  # Partition coefficients : QSAR model 
  parms_PC  = PC_qsar_model(log_Kow =log_Kow, 
                            PC_bw = PC_bw,
                            water_liver = water_liver  , 
                            lipids_liver  = lipids_liver ,   
                            PC_l  = PC_l,
                            water_gonads= water_gonads , 
                            lipids_gonads = lipids_gonads,   
                            PC_go = PC_go,
                            water_brain = water_brain  , 
                            lipids_brain  = lipids_brain, 
                            PC_b  = PC_b,
                            water_fat   = water_fat    , 
                            lipids_fat    = lipids_fat,      
                            PC_f  = PC_f,
                            water_skin  = water_skin   , 
                            lipids_skin   = lipids_skin,     
                            PC_s  = PC_s,
                            water_GIT   = water_GIT    , 
                            lipids_GIT    = lipids_GIT,      
                            PC_gi = PC_gi,
                            water_kidney= water_kidney , 
                            lipids_kidney = lipids_kidney,   
                            PC_k  = PC_k,
                            water_rp    = water_rp     , 
                            lipids_rp     = lipids_rp,  
                            PC_r  = PC_r,
                            water_pp    = water_pp     ,  
                            lipids_pp     = lipids_pp,
                            PC_p  = PC_p)
  
  
  # chemical-dependent - chemical-specific kinetic constants and distribution modifiers. Non relevant parameters can be set to NA
  parms_chemical <- c(
    
    # Metabolism rates
    "Cl_liver" = Cl_liver, # Liver clearance (if used instead of Vmax/Km)
    "rate_plasma" = rate_plasma, # Plasma clearance (mg/d/g)
    "Km" = Km, # Michaelis-Menten constant (µg/mL)
    "Vmax" = Vmax, # Maximum metabolism rate (µg/d/mL liver)
    
    # Binding & distribution
    "log_Kow" = log_Kow, # Octanol-water partition coefficient
    "Unbound_fraction" = Unbound_fraction, # Free fraction in blood
    "Ratio_blood_plasma" = Ratio_blood_plasma, # Blood:plasma concentration ratio
    
    # Excretion rates
    "Ke_bile" = Ke_bile,               # Biliary excretion rate (1/d)
    "Ke_urine" = Ke_urine,             # Urinary excretion rate (1/d)
    "Ke_water" = Ke_water,             # Gills/water excretion rate (1/d)
    "Ke_feces" = Ke_feces,             # Fecal excretion rate (1/d)
    "K_BG" = K_BG,                     # Bile-to-GIT transfer rate (1/d)
    
    # Venous blood flow distribution
    "a_Fpp" = a_Fpp,  # Venous return fraction from poorly perfused tissue
    "a_Fs"   = a_Fs, # Venous return fraction from skin
    
    
    # Oral absorption 
    "Ku" = Ku, # Gut uptake rate constant
    "frac_abs" = frac_abs)	# Total absorbed fraction (bypasses Ku if set)				  
  
  
  # exposure parameters - Environmental and dosing-related settings that control model dynamics over time
  parms_exposure <- c(# environmental parameters
    "TA" = TA,                        # Arrhenius temperature (K)
    "TR_DEB" = TR_DEB,                # Reference temp for DEB growth (K)
    "TR_Fcard" = TR_Fcard,            # Reference temp for cardiac output (K)
    "TR_VO2" = TR_VO2,                # Reference temp for VO2 (K)
    "f_cst" = f_cst,                  # Feeding level: 1 = full feeding
    "V_water" = V_water,              # Volume of the aquarium (mL)
    
    # Exposure quantity
    "IngestQuantity" = IngestQuantity,   # Oral dosing (µg)
    "WaterQuantity" = WaterQuantity,     # Waterborne dose (µg)
    "ivQuantity" = ivQuantity,           # IV bolus (µg)
    
    # Exposure time scenario
    "Time_end_exposure" = Time_end_exposure, # End of exposure (d)
    "period" = period,                       # Interval between doses (d)
    "time_final_dose" = time_final_dose,     # Time of final dose
    "time_first_dose" = time_first_dose      # Time of first dose
  )
  
  # Final Parameter Collection to be passed to the ODE solver
  Parameters <-c(parms_physio, parms_PC, parms_chemical, parms_exposure) 
  
  
  #===============================================
  #(II.1) Initial values
  # This section defines the initial state of the system at time t = 0. It includes chemical quantities (Q_*) in each compartment, initial body structure (L), and other model state variables like urine volume or metabolite accumulation.
  #===============================================
  
  # Initial Structural Length
  L0 = ((Bw_i/ a_BW_L) ^ (1/b_BW_L)) * DEB_shape * 10  # g --> cm Ltotale --> Lstruc --> mm. Initial body weight (Bw_i), Empirical allometric relationship (a_BW_L, b_BW_L), DEB_shape: shape correction length into structural volume, Multiplied by 10 to convert cm to mm.
  
  # State Variable Initialization (holds the initial values of all state variables, typically quantiaty of chemical, ug, mg etc.)
  Initial.Values <- c(
                    # Exposure Medium & Body Size
                      "Q_water" = WaterQuantity,   # Initial chemical quantity in water (µg)
                      "L"       = L0,              # Initial fish length (mm)
                      "V_urine" = 0,               # Initial urine volume (mL)
                      
                      # Circulatory System
                      "Q_art"  = 0,  # Arterial blood chemical load (µg) - Q_art is 0, consistent with bolus delivery via veins
                      "Q_ven"   = ivQuantity, # Venous blood: IV dose goes here at t = 0 - initialized with chemical if using IV administration
                      
                      # Tissue Compartments (Q-quantity in tissue at start of experiment, here with oral administration, where all other compartments starts empty)
                      "Q_brain"= 0, 
                      "Q_fat"=0, 
                      "Q_GIT" =0, 
                      "Q_lumen_GIT" = (IngestQuantity * frac_abs),  # Q_lumen_GIT contains the orally administered chemical, assuming instantaneous dosing and immediate availability.
                      "Q_gonads"= 0, 
                      "Q_kidney" =0,  
                      "Q_liver"= 0,  
                      "Q_skin" =0, 
                      "Q_rp"=0, 
                      "Q_pp"=0, 
                      
                      # Uptake & Administration Tracking. bookkeeping variables for tracking uptake routes. Q_admin_GIT is duplicated to align with mass balance terms
                      "Q_admin_gills"=0, # Cumulative uptake via gills (µg)
                      "Q_admin_GIT"=(IngestQuantity * frac_abs), # Redundant with lumen for tracking total administered via GIT
                      
                      # Metabolism
                      "Q_met_plasma"=0, # Metabolized in plasma (µg)
                      "Q_met_liver"=0, # Metabolized in liver (µg)
                      "Q_met"=0, # Total metabolized chemical (plasma + liver)
                      
                      # Excretion. rack different elimination pathways of the chemical.Start at 0 to measure cumulative amounts over time
                      "Q_elim_water"     = 0,  # Eliminated via gills into water
                      "Q_excret"         = 0,  # Total excreted (gills + urine + feces)
                      "Q_feces"          = 0,  # Fecal excretion
                      "Q_excret_gills"   = 0,  # Excreted via gills
                      "Q_bile"           = 0,  # Stored in bile (not necessarily excreted yet)
                      "Q_urine"          = 0,  # Currently excreted urine
                      "Q_urine_cum"      = 0   # Cumulative urine excretion
                      )
  
  #===============================================
  #(II.2) Environmental and biological scenario
  ## objective: external time-varying inputs for the model simulation, including: 1) Temperature (can be constant or time-varying), 2) Optional urine discharge events (if intermittent urination is modeled)
  #===============================================
  
  # Temperature scenario. 
  if( length (Temperature) > 1) { # If the Temperature input is a vector (i.e., varies over time)...
    times = seq(0, length(Temperature)-1, 1) # Create a time vector (0 to N-1 days)
    import = Temperature # Use the provided time-varying temperature values
  } else { # If temperature is constant...
    times = Times # Use the simulation time vector
    import = rep(Temperature, length(times))} # Repeat the same temperature for all time points
  
  # Combine time and temperature values into a data frame. objective: 1)Create an interpolation function for T  over time, 2) return the nearest value if 't' is outside the range
  Scenario = data.frame(times = times, import = import)
  Temperature.Scenario = approxfun(Scenario$times, Scenario$import, rule = 2)
  
  # scenario for urine discharge (optinonal)
  if (!is.null(urination_interval)){ # If a urination interval is defined by the user (e.g., every 30 minutes)...
    events_urine <- list(data = rbind(data.frame(var = c("V_urine"), # Create a sequence of times from the beginning to the end of simulation. At each of these times, we reset urine volume and chemical amount to 0
                                                 time = seq(Times[1], rev(Times)[1] , by=urination_interval), 
                                                 value = 0, 
                                                 method = c("replace")),
                                      data.frame(var = c("Q_urine"), 
                                                 time = seq(Times[1], rev(Times)[1] , by=urination_interval), 
                                                 value = 0, 
                                                 method = c("replace"))))
  }else{ events_urine <- NULL } # If no urination interval is defined, do not include any urine-resetting events
  
  #===============================================
  #(II.3) Exposure scenario :   events_repeated
  #===============================================
  events_repeated_f <- NULL
  
  # repeated dose scenario for oral bolus
  if (!is.null(period) & IngestQuantity!=0 ) {
    
    # no stop before the end of the experiement (If no end time is specified for dosing, set it to the last multiple of 'period' before the simulation ends)
    if (is.null(time_final_dose)){ time_final_dose <- floor(max(Times)/period)*period}
    
    # First dose NOT at t0 of the experiment (If first dose is NOT at time 0, reset initial oral dose amounts to 0)
    if (!is.null(time_first_dose)) { Initial.Values["Q_admin_GIT"]<-0 ; Initial.Values["Q_lumen_GIT"]<-0 }
    
    # First dose at t0 of the experiment (If first dose time is not defined, assume it happens at 'period' (i.e., first interval))
    if (is.null(time_first_dose)){ time_first_dose = period }
    
    # events : a list with a data frame for each compartment --> Cpt name / Time / Dose / method
    events_repeated_f <- list(data = rbind(data.frame(var = c("Q_lumen_GIT"), 
                                                      time = seq(time_first_dose, time_final_dose , by=period), 
                                                      value = as.numeric(c(parms_chemical["frac_abs"]*parms_exposure["IngestQuantity"])), 
                                                      method = c("add")),
                                           data.frame(var = c("Q_admin_GIT"), 
                                                      time = seq(time_first_dose, time_final_dose , by=period), 
                                                      value = as.numeric(c(parms_chemical["frac_abs"]*parms_exposure["IngestQuantity"])), 
                                                      method = c("add"))
    )) 
  } 
  
  # repeated dose scenario for water exposure (repeated bolus adminstration - experimental)
  
  events_repeated_w <- NULL # Initialize an empty event list for repeated oral dosing
  
  
  if (!is.null(period) & WaterQuantity!=0 ) { # If a dosing period is defined AND ingestion amount is non-zero, then set up repeated oral doses
    
    # no stop before the end of the experiment
    if (is.null(time_final_dose)){ time_final_dose <- floor(max(Times)/period)*period}
    
    # First dose NOT at t0 of the experiment
    if (!is.null(time_first_dose)) { Initial.Values["Q_water"]<-0 }
    
    # First dose at t0 of the experiment
    if (is.null(time_first_dose)){ time_first_dose = period }
    
    # events : a list with a data frame for each compartment --> Cpt name / Time / Dose / method
    events_repeated_w <- list(data = rbind(data.frame(var = c("Q_water"), # 1. Add chemical to gut lumen (Q_lumen_GIT) at every dose time
                                                      time  = seq(time_first_dose, time_final_dose , by=period), 
                                                      value = (1-frac_renewed), 
                                                      method = c("multiply")),
                                           data.frame(var = c("Q_water"), # 2. Add to cumulative dose tracker (Q_admin_GIT) at same times
                                                      time =  0.0000001+seq(time_first_dose, time_final_dose , by=period), 
                                                      value = as.numeric( parms_exposure["WaterQuantity"]* frac_renewed), 
                                                      method = c("add")),
                                           data.frame(var = c("Q_water"),
                                                      time  = Time_end_exposure,
                                                      value = 0,
                                                      method = c("multiply"))
                                           
                                           #                                            # In depuration phase, water renewal
                                           #                                             ,	 data.frame(var = c("Q_water"),
                                           #                                                        time  = 0.0000001+seq(Time_end_exposure, max(Times) , by=period),
                                           #                                                        value = (1-frac_renewed),
                                           #                                                        method = c("multiply"))
    ))    
  }
  
  if ((is.null(period))& (WaterQuantity!=0) ){
    events_repeated_w <- list(data = rbind(data.frame(var = c("Q_water"), 
                                                      time  = Time_end_exposure, 
                                                      value = 0, 
                                                      method = c("multiply"))))
  }
  
  
  all_event_data <- Filter(Negate(is.null), list(
    if (!is.null(events_repeated_f)) events_repeated_f$data else NULL,
    if (!is.null(events_repeated_w)) events_repeated_w$data else NULL,
    if (!is.null(events_urine)) events_urine$data else NULL
  ))
  
  # Slå sammen og sorter eventene etter tid
  events_combined <- do.call(rbind, all_event_data)
  events_combined <- events_combined[order(events_combined$time), ]
  
  # Lag event-strukturen til ode()
  events <- list(data = events_combined)
  
  
  # events<-list(data = rbind(events_repeated_f[[1]], events_repeated_w[[1]], events_urine[[1]]))
  
  
  #===============================================
  # (III) Numerical resolution and post-calculation
  " objective: uns the PBPK model simulation, solves the ODEs, and returns the simulation output along with partition coefficients and fish body length over time"
  #===============================================
  
out <- ode(
  times  = Times,              # Time points at which output is desired (e.g. seq(0, 20, 0.1))
  func   = PBPK.model,         # The model function that calculates derivatives (defined earlier)
  y      = Initial.Values,     # Vector of initial state values (Q_water, Q_ven, Q_liver, etc.)
  parms  = Parameters,         # All parameters needed by the model (combined in earlier steps)
  input.T = Temperature.Scenario, # Temperature interpolation function passed to the model
  method = "lsodes",           # Solver method: suitable for stiff systems with many compartments
  events = events              # List of dosing and urination events (if defined earlier)
)

    PC = parms_PC
  
  return( list("Sim"=out, # Simulation results (time, state variable values)
               "PC"=PC, # Partition coefficients used in the model
               "Length" = out[,"L"] / DEB_shape * 1/10 ) ) # Convert structural length (L, mm) to total length (cm)
}