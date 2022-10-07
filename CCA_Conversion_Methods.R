#### Set up environment ----
rm(list = ls()) ; options(cores = parallel::detectCores())

# Packages
library(readxl)
library(tidyverse)

# Dataset
CCA_CC_SC <- read_excel("Data/CCA_CC_SC.xlsx")

# Units to convert
mol_to_umol      = 1e+06
nmol_to_umol     = 1e+03
g_to_mg          = 1e+03
molar_mass_CaCO3 = 100.1 
day_to_hour      = 24.00
hour_to_minute   = 60.00

#### Reshape methodologies ----
# Rename Methods from Steeve, Ben and Chris paper – Global Change Biology
CCA_CC_SC$Method_family = ifelse(CCA_CC_SC$Method %in% unique(CCA_CC_SC$Method)[c(3,12,25)], 
                                 "TA anomaly", CCA_CC_SC$Method)
CCA_CC_SC$Method_family = ifelse(CCA_CC_SC$Method %in% unique(CCA_CC_SC$Method)[c(4,10,13)], 
                                 "Staining"  , CCA_CC_SC$Method_family)
CCA_CC_SC$Method_family = ifelse(CCA_CC_SC$Method %in% unique(CCA_CC_SC$Method)[c(5,07,21)], 
                                 "Isotopes"  , CCA_CC_SC$Method_family)
CCA_CC_SC$Method_family = ifelse(CCA_CC_SC$Method %in% unique(CCA_CC_SC$Method)[c(1,09,11,14,20,22,23,26)], 
                                 "BW or RGR" , CCA_CC_SC$Method_family)
CCA_CC_SC = CCA_CC_SC %>% mutate(., Method = replace_na(Method, "NA")) %>% 
  dplyr::filter(., Method != "Net photosynthesis", Method != "Photo", Method != "Net respiration", Method != "SOPHIE TO GET ?")

#### Biomass ----
biomass_corrected_unit <- unique(CCA_CC_SC$Unit)[c(3,5,9,13,15,17,21:24,29:30,33,37,43,45:46,53:54,56:57)]
CCA_biom <- CCA_CC_SC %>% dplyr::filter(., Unit %in% biomass_corrected_unit)
CCA_biom$Unit = ifelse(CCA_biom$Unit %in% unique(CCA_biom$Unit)[c(1,6:7,9:10,13:14,18:19,21)], "umol_g_h", CCA_biom$Unit)
CCA_biom <- CCA_biom %>% 
  mutate(std_biom_corrected_unit = case_when(Unit == "umol_g_h" ~ Rate,
                                             Unit == "nmoles g-1 dry weight min-l" ~ Rate / nmol_to_umol * hour_to_minute,
                                             Unit == "mgCaCO3 g-1 d-1" ~ Rate / g_to_mg / molar_mass_CaCO3 / day_to_hour * mol_to_umol,
                                             Unit == "g CaCO3/g/h" ~ Rate / molar_mass_CaCO3 * mol_to_umol,
                                             Unit == "umol/g(dwt)/min" ~ Rate * hour_to_minute,
                                             Unit == "mg CaCO3 mg-1 day-1." ~ Rate / g_to_mg / molar_mass_CaCO3 / day_to_hour * mol_to_umol,
                                             Unit == "mg wet weight mg-1 day-1" ~ Rate / molar_mass_CaCO3 / day_to_hour * mol_to_umol,
                                             Unit == "g per g (day-1)" ~ Rate / molar_mass_CaCO3 / day_to_hour * mol_to_umol,
                                             Unit == "µmolC/g day" ~ Rate / day_to_hour,
                                             Unit == "µmol CaCO3/g/day" ~ Rate / day_to_hour,
                                             Unit == "mg g–1 d–1)" ~ Rate / g_to_mg /day_to_hour * mol_to_umol,
                                             Unit == "mgCaCO3 g-1 DW h-1" ~ Rate / g_to_mg * mol_to_umol,
                                             Unit == "µmol CaCO3/g/day" ~ Rate / day_to_hour)) 

## Mistakes
# Too much incertitude when it's wet ––– Remove.
CCA_biom <- CCA_biom %>% dplyr::filter(., !grepl("wet",`Unit`), !grepl("Wet",`Method`))
# Error unit in Graba-Landry et al. 2018 ––– in day and not in hours.
CCA_biom$std_biom_corrected_unit[c(32,33)] <- CCA_biom$std_biom_corrected_unit[c(32,33)] / 24
# Error unit in Donham et al 2022 ––– Percentage, useless - Remove.
CCA_biom <- CCA_biom %>% dplyr::filter(., `Paper name` != "Donham et al 2022")
# Suspicious w/ Barner et al. 2018 ––– Remove (waiting for this where they found the values)
CCA_biom <- CCA_biom %>% dplyr::filter(., `Paper name` != "Barner et al. 2018")
# Remove negative values and recode NA values
CCA_biom <- CCA_biom %>% dplyr::filter(., std_biom_corrected_unit >= 0)
CCA_biom$Method_family[CCA_biom$Method_family == "NA"] <- NA
CCA_biom$Climate[CCA_biom$Climate == "NA"] <- NA
# Add Genus informations
CCA_biom$Genus   <- stringr::word(CCA_biom$Species, 1)
CCA_biom$Species <- paste(stringr::word(CCA_biom$Species, 1), stringr::word(CCA_biom$Species, 2), sep = "_")

## Fill the gaps – Methods
# Legrand et al. 2019
CCA_biom$Method_family[c(33:36)] <- "TA anomaly"
# Navarte et al. 2019, 2020
CCA_biom$Method[c(39:41)] <- "RGR" ; CCA_biom$Method_family[c(39:41)] <- "BW or RGR"
# Qui-Minet et al. 2019
CCA_biom$Method[c(44:54)] <- "Total alkalinity" ; CCA_biom$Method_family[c(44:54)] <- "TA anomaly"
# Sordo et al. 2018
CCA_biom$Method[c(59:62)] <- "BW" ; CCA_biom$Method_family[c(59:62)] <- "BW or RGR"

## Fill the gaps – Climate
CCA_biom$Climate[which(CCA_biom$`Paper name` %in% c("Navarte et al 2020", "Noisette et al 2013 JEMBE", 
                                                   "Ragazzola et al 2020", "Schubert et al 2021", 
                                                   "Schubert et al. 2019", "Sordo et al 2020"))] = "Cool Temperate"
CCA_biom$Climate[which(CCA_biom$`Paper name` == "McGraw et al. 2010")] = "Warm temperate"
CCA_biom$Climate[which(CCA_biom$`Paper name` == "Wei et al 2020")] = "Tropical"

## Some values are really high... Discuss about it w/ Steeve, Chris and Ben
table(CCA_biom$Method_family, CCA_biom$Genus)
Biom_stat = CCA_biom %>% group_by(Method_family, Climate) %>% 
  summarise(growth_rate = mean(std_biom_corrected_unit), SD = sd(std_biom_corrected_unit))

CCA_biom %>% dplyr::filter(., std_biom_corrected_unit <= 10) %>% 
ggplot() + geom_boxplot(aes(x = Method_family, y = std_biom_corrected_unit)) 


#### Surface ----
surface_corrected_unit <- unique(CCA_CC_SC$Unit)[c(1,8,11,13:14,16,20,30,33,49,57,61)]
CCA_surf <- CCA_CC_SC %>% dplyr::filter(., Unit %in% surface_corrected_unit)

table(CCA_surf$Unit)

table(CCA_surf$Method, CCA_surf$`Paper name`) %>% View()

CCA_surf <- CCA_surf %>% drop_na(Rate) %>% dplyr::filter(., !grepl("\"CCA\"", Species)) %>% 
  dplyr::filter(., !grepl("CCA", Species))
unique(CCA_surf$Species)
