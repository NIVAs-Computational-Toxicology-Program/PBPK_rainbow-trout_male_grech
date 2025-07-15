#################
## Example TCDD, rt
## simualtions: pulsed exposure, continous exposure
###################
###############################
####### 2,3,7,8-Tetrachlorodibenzo-p-dioxin (TCDD)

###########################################################
# APPLICATION RAINOW TROUT
##### Data from Branson, 1985
##### Water exposure : 6 h of exposure + 139 days depuration
########################################################### 

#===============================================
# Observed data
#===============================================

Res_times = c(0.083, 0.25, 7, 22, 42, 64, 78, 118, 134, 136, 139) # days

Times_water = c(0.25, 0.30, 0.75, 1, 1.25,  1.5, 1.75, 2, 3, 4.17, 5, 6) / 24 # h/24 = days
Mean_water = c(362, 330, 337, 320, 321, 354, 316, 327, 307, 296, 310, 291) * 0.63

Mean_value_fish = (c(1.01, 2.58, 2.06, 1.97, 0.98, 0.78, 0.81, 0.583, 0.65, 0.82, 0.49 ) * 10^(-3)) # ng/g * 10^(-3) = microg/g
SD_value_fish = (c(0.09, 0.06, 0.15, 0.07, 0.04, 0.05, NA, 0.1, NA, NA, NA) * 10^(-3)) # ng/g * 10^(-3) = microg/g

#===============================================
# Simulation - example puled exposure, water
#===============================================

Quantity_ext = 107 * 10^(-6) # ng/L * 10^(-6) = microg/mL

# tps <- round(sort(unique(c(seq(0, 140, 0.1), Res_times, seq(0, 140, 1)))), 8)
tps = sort( unique( c(seq(0,140,0.1), Res_times)))  # pulse: sort( unique( c(seq(0,140,0.1), Res_times))), continous: pulse: sort( unique( c(seq(0,140,1), Res_times)))

res_branson <- PBPK.RT.M( WaterQuantity = Quantity_ext * (145000/30), 
                          IngestQuantity    = 0,
                          ivQuantity        = 0,
                          Bw_i          = 35,
                          log_Kow       = 6.8,
                          Temperature   = 13,  
                          V_water       = (145000/30), # mL
                          Cl_liver      = 0,     # mL/g/d
                          K_BG          = 0,                        
                          #Ke_water     = (2 /145) * (24 * 60) * 0.035 , # que dans depuration phase
                          Ke_bile       = 1.5, # d-1
                          Ku            = (0.58*0.83)/(1-0.58), # d-1 
                          frac_abs      = 1,
                          
                          Times         = tps,
                          time_first_dose = NULL, # pulse: NULL, continous: 0
                          time_final_dose = NULL, # pulse: NULL, continous: 140
                          Time_end_exposure= 6/24, # pulse: 6/24, continous: 140 
                          period        = 140, # pulse: 140, continous: 1
                          frac_renewed  = 1)


##############
## data wrangling 
######

result <- res_branson  # eller hva du har kalt objektet

sim_data_branson <- as.data.frame(result$Sim)   # Hele simuleringen som data.frame
pc_values_branson <- as.data.frame(result$PC)                  # Partisjonskoeffisienter (named vector)
fish_length_branson <- as.data.frame(result$Length)            # Lengde over tid

# explore data - Simulation
print(sim_data_branson)
str(sim_data_branson)
names(sim_data_branson) # names of the columns
head(sim_data_branson) # first rows of the simulation data

# explore data - partioning tissues
print(pc_values_branson)
str(pc_values_branson)
names(pc_values_branson) # names of the columns
head(pc_values_branson) # first rows of the partitioning data

# explore data - partioning tissues
print(fish_length_branson)
str(fish_length_branson)
names(fish_length_branson) # names of the columns
head(fish_length_branson) # first rows of the partitioning data

