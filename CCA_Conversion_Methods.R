library(readxl)
library(tidyverse)

CCA_CC_SC <- read_excel("Data/CCA_CC_SC.xlsx")

# Units to convert
nmol_to_umol     = 1e+03
g_to_mg          = 1e+03
molar_mass_CaCO3 = 100.1 
day_to_hour      = 24.00
hour_to_minute   = 60.00

# Rename Methods
CCA_CC_SC$Method_family = ifelse(CCA_CC_SC$Method %in% unique(CCA_CC_SC$Method)[c(3,12,25)]               , "TA anomaly", CCA_CC_SC$Method)
CCA_CC_SC$Method_family = ifelse(CCA_CC_SC$Method %in% unique(CCA_CC_SC$Method)[c(4,10,13)]               , "Staining"  , CCA_CC_SC$Method_family)
CCA_CC_SC$Method_family = ifelse(CCA_CC_SC$Method %in% unique(CCA_CC_SC$Method)[c(5,07,21)]               , "Isotopes"  , CCA_CC_SC$Method_family)
CCA_CC_SC$Method_family = ifelse(CCA_CC_SC$Method %in% unique(CCA_CC_SC$Method)[c(1,09,11,14,20,22,23,26)], "BW or RGR" , CCA_CC_SC$Method_family)
CCA_CC_SC = CCA_CC_SC %>% mutate(., Method = replace_na(Method, "NA")) %>% 
  dplyr::filter(., Method != "Net photosynthesis", Method != "Photo", Method != "Net respiration", Method != "SOPHIE TO GET ?")

# Biomass
biomass_corrected_unit <- unique(CCA_CC_SC$Unit)[c(4,6,10,17,19,23:26,31:32,35,40,45,47:48,55:56,58:59)]
CCA_biom <- CCA_CC_SC %>% dplyr::filter(., Unit %in% biomass_corrected_unit)
CCA_biom$Unit = ifelse(CCA_biom$Unit %in% unique(CCA_biom$Unit)[c(1,5:6,8:9,12:13,17:18,20)], "umol_g_h", CCA_biom$Unit)
CCA_biom = CCA_biom %>% 
  mutate(std_biom_corrected_unit = case_when(Unit == "umol_g_h" ~ Rate,
                                             Unit == "nmoles g-1 dry weight min-l" ~ Rate / nmol_to_umol * hour_to_minute,
                                             Unit == "mgCaCO3 g-1 d-1" ~ Rate * g_to_mg / day_to_hour,
                                             Unit == "g CaCO3/g/h" ~ Rate / molar_mass_CaCO3_g_per_mol,
                                             Unit == "umol/g(dwt)/min" ~ Rate * hour_to_minute,
                                             Unit == "mg CaCO3 mg-1 day-1." ~ Rate / molar_mass_CaCO3_g_per_mol / day_to_hour,
                                             Unit == "mg wet weight mg-1 day-1" ~ Rate / molar_mass_CaCO3_g_per_mol / day_to_hour,
                                             Unit == "g per g (day-1)" ~ Rate / molar_mass_CaCO3_g_per_mol / day_to_hour,
                                             Unit == "µmolC/g day" ~ Rate / day_to_hour,
                                             Unit == "µmol CaCO3/g/day" ~ Rate / day_to_hour,
                                             Unit == "µmol CaCO3/g/day" ~ Rate / g_to_mg / molar_mass_CaCO3_g_per_mol))

CCA_CC_SC %>% filter(Method == "NA") %>% View()
unique(CCA_CC_SC$Method)

# Surface
surface_corrected_unit <- unique(CCA_CC_SC$Unit)[c(1,8,11,13:14,16,20,30,33,49,57,61)]
CCA_surf <- CCA_CC_SC %>% dplyr::filter(., Unit %in% surface_corrected_unit)