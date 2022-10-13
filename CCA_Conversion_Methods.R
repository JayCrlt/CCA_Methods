#### Set up environment ----
rm(list = ls()) ; options(cores = parallel::detectCores())

# Packages
library(readxl)
library(tidyverse)

# Useful functions
'%notin%' <- function(x,y)!('%in%'(x,y))

# Dataset
CCA_CC_SC <- read_excel("Data/CCA_CC_SC.xlsx")

# Units to convert
mol_to_umol      = 1e+006
nmol_to_umol     = 1e+003
g_to_mg          = 1e+003
molar_mass_CaCO3 = 100.10 
day_to_hour      = 24e+00
hour_to_minute   = 60e+00
year_to_day      = 365.25
m2_to_cm2        = 1e+004
calcite_density  = 2.7100
year_to_month    = 12e+00
mm_to_um         = 1e+003
cm_to_mm         = 1e+001

# Biomass and surface units vectors
biomass_corrected_unit <- c("µmol/g/h", "nmoles g-1 dry weight min-l", "mgCaCO3 g-1 d-1", "mg g–1 d–1)", "g CaCO3/g/h", "ulmol CaCO3 g-1 h-", 
                            "umol C h -1 g -1\002","umol/g(dwt)/min", "umol_g_h", "µmol CaCO3/g/h", "mg CaCO3 mg-1 day-1.", 
                            "mg wet weight mg-1 day-1", "[myu]molCaCO3g^-1hr^-1", "umol CaCO3 gDW-1 h-1", "g per g (day-1)", "µmolC/g day", 
                            "µmol CaCO3/g/day", "Net Calcification (μmol CaCO3 g−1 DW  hr−1)", "Dark Net Calcification (μmol CaCO3 g−1 DW  hr−1)", 
                            "mgCaCO3 g-1 DW h-1", "[myu]molCaCO3gFW^-1hr^-1")

surface_corrected_unit <- c("mg/cm2/d", "Calcif in umol m-2 h-1", "mgCaCO3 cm-2 d-1", "CaCO3 [mg/cm**2/day]", "mg CaCO3 cm–2 d–1)", "umol CaCO3 cm-2 h-1", 
                            "mmol/cm**2/day", "mg cm-2 day-1", "micro mol CaCO3 cm-2 h-1", "µmol CaCO3/cm2/h", "mg/cm-2/year", "CaCO3 mg cm–2 d–1")

biomass_umol_g_h       <- c("µmol/g/h", "ulmol CaCO3 g-1 h-", "umol C h -1 g -1\002", "umol_g_h", "µmol CaCO3/g/h", "[myu]molCaCO3g^-1hr^-1", 
                            "umol CaCO3 gDW-1 h-1", "Net Calcification (μmol CaCO3 g−1 DW  hr−1)", "Dark Net Calcification (μmol CaCO3 g−1 DW  hr−1)", 
                            "[myu]molCaCO3gFW^-1hr^-1")

surface_mg_cm2_day     <- c("mg/cm2/d", "mgCaCO3 cm-2 d-1", "CaCO3 [mg/cm**2/day]", "mg CaCO3 cm–2 d–1)", "mg cm-2 day-1", "CaCO3 mg cm–2 d–1")
extension_mm_year      <- c("mm y-1", "mm-year", "mm/year")

#### Reshape methodologies ----
# Rename Methods from Steeve, Ben and Chris paper – Global Change Biology
CCA_CC_SC$Method_family = NA
CCA_CC_SC$Method_family[which(CCA_CC_SC$Method %in% c("TA", "TA anomaly", "TA anamoly Light calcification"))] = "TA anomaly"
CCA_CC_SC$Method_family[which(CCA_CC_SC$Method %in% c("Staining", "Calcofluor white", "alizarin red"))] = "Staining"
CCA_CC_SC$Method_family[which(CCA_CC_SC$Method %in% c("45Ca", "Isotope/growth band", "13C"))] = "Isotopes"
CCA_CC_SC$Method_family[which(CCA_CC_SC$Method %in% c("BW", "RGR Wet weight", "change in total wieght", "Linear extension", 
                                                      "RGR (SA)", "Relative growth rate", "length extension", 
                                                      "Area extension"))] = "BW or RGR"
CCA_CC_SC = CCA_CC_SC %>% mutate(., Method = replace_na(Method, "NA"), Method_family = replace_na(Method_family, "NA")) %>% 
  dplyr::filter(., Method != "Net photosynthesis", Method != "Photo", Method != "Net respiration", Method != "SOPHIE TO GET ?")

#### Biomass ----
CCA_biom <- CCA_CC_SC %>% dplyr::filter(., Unit %in% biomass_corrected_unit)
CCA_biom$Unit = ifelse(CCA_biom$Unit %in% biomass_umol_g_h, "umol_g_h", CCA_biom$Unit)
CCA_biom <- CCA_biom %>% 
  mutate(std_biom_corrected_unit = case_when(Unit == "umol_g_h" ~                    Rate,
                                             Unit == "µmolC/g day" ~                 Rate                              / day_to_hour,
                                             Unit == "µmol CaCO3/g/day" ~            Rate                              / day_to_hour,
                                             Unit == "umol/g(dwt)/min" ~             Rate                              * hour_to_minute,
                                             Unit == "nmoles g-1 dry weight min-l" ~ Rate                              * hour_to_minute   / nmol_to_umol,
                                             Unit == "g CaCO3/g/h" ~                 Rate / molar_mass_CaCO3                              * mol_to_umol,
                                             Unit == "mgCaCO3 g-1 DW h-1" ~          Rate / molar_mass_CaCO3 / g_to_mg                    * mol_to_umol,
                                             Unit == "mgCaCO3 g-1 d-1" ~             Rate / molar_mass_CaCO3 / g_to_mg / day_to_hour      * mol_to_umol,
                                             Unit == "mg g–1 d–1)" ~                 Rate / molar_mass_CaCO3 / g_to_mg / day_to_hour      * mol_to_umol,
                                             Unit == "g per g (day-1)" ~             Rate / molar_mass_CaCO3           / day_to_hour      * mol_to_umol,
                                             Unit == "mg CaCO3 mg-1 day-1." ~        Rate / molar_mass_CaCO3           / day_to_hour      * mol_to_umol,
                                             Unit == "mg wet weight mg-1 day-1" ~    Rate / molar_mass_CaCO3           / day_to_hour      * mol_to_umol)) 

# Check where there is a lack of information
table(CCA_biom$Method_family, CCA_biom$Method)

## Mistakes or doubts
# Too much incertitude when it's wet ––– Remove.
CCA_biom <- CCA_biom %>% dplyr::filter(., !grepl("wet",`Unit`))
# Error unit in Graba-Landry et al. 2018 ––– in day and not in hours.
CCA_biom$std_biom_corrected_unit[which(CCA_biom$`Paper name` == "Graba-Landry et al. 2018")] = 
  CCA_biom$std_biom_corrected_unit[which(CCA_biom$`Paper name` == "Graba-Landry et al. 2018")] / 24
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
# BW or RGR
CCA_biom$Method[which(CCA_biom$`Paper name` %in% c("Narvarte et al. 2019", "Navarte et al 2020", "Sordo et al. 2018"))] <- "RGR"
CCA_biom$Method_family[which(CCA_biom$`Paper name` %in% c("Narvarte et al. 2019", "Navarte et al 2020", "Sordo et al. 2018"))] <- "BW or RGR"
# TA Anomaly
CCA_biom$Method[which(CCA_biom$`Paper name` %in% c("Legrand et al 2019", "Qui-Minet et al. 2019"))] <- "Total alkalinity"
CCA_biom$Method_family[which(CCA_biom$`Paper name` %in% c("Legrand et al 2019", "Qui-Minet et al. 2019"))] <- "TA anomaly"

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
CCA_surf <- CCA_CC_SC %>% dplyr::filter(., Unit %in% surface_corrected_unit) %>% drop_na(Rate) %>% 
  dplyr::filter(., !grepl("\"CCA\"", Species), !grepl("CCA", Species))
CCA_surf$Unit = ifelse(CCA_surf$Unit %in% surface_mg_cm2_day, "mg_cm2_day", CCA_surf$Unit)

# Mistake pre-analysis with Johnson & Carpenter 2012
CCA_surf$Rate[which(CCA_surf$`Paper name` == "Johnson and Carpenter 2012")] <- c(0.0059, 0.0062)
CCA_surf$Error[which(CCA_surf$`Paper name` == "Johnson and Carpenter 2012")] <- c(0.0003, 0.00025)
CCA_surf$Species[which(CCA_surf$`Paper name` == "Johnson and Carpenter 2012")] <- "Hydrolithon onkodes"

# Std the units
CCA_surf <- CCA_surf %>% 
  mutate(std_surf_corrected_unit = case_when(Unit == "mg_cm2_day" ~               Rate,
                                             Unit == "mg/cm-2/year" ~             Rate                                  / year_to_day,
                                             Unit == "mmol/cm**2/day" ~           Rate * molar_mass_CaCO3,
                                             Unit == "µmol CaCO3/cm2/h" ~         Rate * molar_mass_CaCO3 / mol_to_umol * day_to_hour * g_to_mg,
                                             Unit == "umol CaCO3 cm-2 h-1" ~      Rate * molar_mass_CaCO3 / mol_to_umol * day_to_hour * g_to_mg,
                                             Unit == "micro mol CaCO3 cm-2 h-1" ~ Rate * molar_mass_CaCO3 / mol_to_umol * day_to_hour * g_to_mg,
                                             Unit == "Calcif in umol m-2 h-1" ~   Rate * molar_mass_CaCO3 / mol_to_umol * day_to_hour * g_to_mg / m2_to_cm2)) 

## Mistakes
# This study is not working with surface std.
CCA_surf <- CCA_surf %>% dplyr::filter(., `Paper name` != "Noisette et al. 2013 J phycol")
# there is a serious doubt to include this study
CCA_surf <- CCA_surf %>% dplyr::filter(., `Paper name` != "Burdett et al. 2018")
# Remove negative values and recode NA values
CCA_surf <- CCA_surf %>% dplyr::filter(., std_surf_corrected_unit >= 0)
CCA_surf$Method_family[CCA_surf$Method_family == "NA"] <- NA
CCA_surf$Climate[CCA_surf$Climate == "NA"] <- NA

## Fill the gaps
CCA_surf$Method_family[which(CCA_surf$`Paper name` == "Schubert et al 2021")] = "X-ray CT Scan"
CCA_surf$Method_family[which(CCA_surf$`Paper name` == "Tanaka et al. 2016")] = "BW or RGR"
CCA_surf$Climate[which(CCA_surf$`Paper name`       %in% c("Schubert et al 2021", "Schubert et al 2022"))] = "Cool Temperate"
CCA_surf$Climate[which(CCA_surf$`Paper name`       == "Tanaka et al. 2016")] = "Tropical"
CCA_surf$Climate[which(CCA_surf$`Paper name`       == "Westfield et al 2022")] = "Polar"
CCA_surf$Method[which(CCA_surf$`Paper name`        == "Schubert et al 2021")] = "X-ray CT Scan"
CCA_surf$Method[which(CCA_surf$`Paper name`        == "Tanaka et al. 2016")] = "RGR"

# Add Genus informations
CCA_surf$Genus   <- stringr::word(CCA_surf$Species, 1)
CCA_surf$Species <- paste(stringr::word(CCA_surf$Species, 1), stringr::word(CCA_surf$Species, 2), sep = "_")

## Some values are really high... Discuss about it w/ Steeve, Chris and Ben
table(CCA_surf$Method_family, CCA_surf$Genus)
Surfstat = CCA_surf %>% group_by(Method_family, Climate) %>% 
  summarise(growth_rate = mean(std_surf_corrected_unit), SD = sd(std_surf_corrected_unit))

CCA_surf %>% ggplot() + geom_boxplot(aes(x = Method_family, y = std_surf_corrected_unit)) 

#### Volume ----
CCA_vol <- CCA_CC_SC %>% dplyr::filter(., Unit == "g CaCO3/cm3/year")
CCA_vol <- CCA_vol %>% mutate(std_biom_corrected_unit = Rate * calcite_density / molar_mass_CaCO3 * mol_to_umol / 
                                year_to_day / day_to_hour,
                              Method_family = rep("X-ray CT Scan", length(CCA_vol$Method)))
CCA_vol$Climate[which(CCA_vol$`Paper name` == "Williams et al 2020")] = "Polar"

# Add Genus informations
CCA_vol$Genus   <- stringr::word(CCA_vol$Species, 1)
CCA_vol$Species <- paste(stringr::word(CCA_vol$Species, 1), stringr::word(CCA_vol$Species, 2), sep = "_")

## Add Volume dataset to Biomass dataset
CCA_biom = rbind(CCA_biom, CCA_vol) %>% data.frame() ; rm(CCA_vol)

#### Non-standardized datasets
CCA_non_std <- CCA_CC_SC %>% dplyr::filter(., Unit %notin% c(biomass_corrected_unit, surface_corrected_unit, "g CaCO3/cm3/year"),
                                           !grepl("%", Unit), Unit %notin% c(NA, "um", "ugCaCO3", "1 hour incubations"),
                                           Unit %notin% c("umol CacO3 h-1","g buoyant wt. yr-1", "mm2 42days-1", "mm per mm (day-1)")) 
                                                                                                   # non-usable (3) | already used (1)

CCA_non_std$Unit = ifelse(CCA_non_std$Unit %in% extension_mm_year, "mm_year", CCA_non_std$Unit)
CCA_non_std <- CCA_non_std %>% 
  mutate(std_extension_corrected_unit = case_when(Unit == "mm_year"     ~ Rate,
                                                  Unit == "mm d-1"      ~ Rate * year_to_day,
                                                  Unit == "mm/day"      ~ Rate * year_to_day,
                                                  Unit == "µm/day"      ~ Rate * year_to_day   / mm_to_um,
                                                  Unit == "cm month-1"  ~ Rate * year_to_month * cm_to_mm,
                                                  Unit == "mm/4 months" ~ Rate * 3)) 

# recode NA values
CCA_non_std$Method_family[CCA_non_std$Method_family == "NA"] <- NA
CCA_non_std$Climate[CCA_non_std$Climate == "NA"] <- NA

## Fill the gaps
CCA_non_std$Climate[which(CCA_non_std$`Paper name` %in% c("McCoy and Ragazzola 2014", "Piazza et al 2022"))] = "Cool Temperate"
CCA_non_std$Method[which(CCA_non_std$`Paper name` == "O'Leary et al. 2017")] = "RGR"
CCA_non_std$Method_family[which(CCA_non_std$`Paper name` == "O'Leary et al. 2017")] = "BW or RGR"

# Add Genus informations
CCA_non_std$Genus   <- stringr::word(CCA_non_std$Species, 1)
CCA_non_std$Species <- paste(stringr::word(CCA_non_std$Species, 1), stringr::word(CCA_non_std$Species, 2), sep = "_")

## Some values are really high... Discuss about it w/ Steeve, Chris and Ben
table(CCA_non_std$Method_family, CCA_non_std$Genus)
nonstdstat = CCA_non_std %>% group_by(Method_family, Climate) %>% 
  summarise(growth_rate = mean(std_extension_corrected_unit), SD = sd(std_extension_corrected_unit))

CCA_non_std %>% ggplot() + geom_boxplot(aes(x = Method_family, y = std_extension_corrected_unit)) 

#### Incomplete studies ----
length(unique(CCA_CC_SC$`Paper name`))                                                  # 83 studies in total
length(unique(c(CCA_non_std$`Paper name`, CCA_biom$Paper.name, CCA_surf$`Paper name`))) # 50 studies used

data_with_NA  <- unique(CCA_CC_SC$`Paper name`[which(CCA_CC_SC$Method == "NA")])
data_complete <- unique(c(CCA_non_std$`Paper name`, CCA_biom$Paper.name, CCA_surf$`Paper name`))
length(data_with_NA[data_with_NA %notin% data_complete])                                # 16 incomplete studies

incomplete_data <- CCA_CC_SC %>% dplyr::filter(., `Paper name` %in% data_with_NA[data_with_NA %notin% data_complete],
                                               !grepl("%", Unit), Temperature != 31.95) %>% drop_na(., Rate)

# Vasquez-Elizondo and Enriquez 2016
incomplete_data$Error[which(incomplete_data$`Paper name` == "Vasquez-Elizondo and Enriquez 2016")] = c(0.01, 0.02, 0.06)
incomplete_data$Method[which(incomplete_data$`Paper name` == "Vasquez-Elizondo and Enriquez 2016")] = "Total alkalinity"
incomplete_data$Method_family[which(incomplete_data$`Paper name` == "Vasquez-Elizondo and Enriquez 2016")] = "TA anomaly"
incomplete_data$Unit[which(incomplete_data$`Paper name` == "Vasquez-Elizondo and Enriquez 2016")] = "µmol CaCO3/cm2/h"
# Studies from Russel et al (2009 and 2012) were non found
incomplete_data <- incomplete_data %>% dplyr::filter(., !grepl("Russel", `Paper name`))
# Studies from Ragazzola (2012, 2013)
incomplete_data$Unit[which(incomplete_data$`Paper name` %in% c("Ragazzola et al. 2012", "Ragazzola et al. 2013"))] = "mm_year"
incomplete_data$Method[which(incomplete_data$`Paper name` %in% c("Ragazzola et al. 2012", "Ragazzola et al. 2013"))] = "Alizarin"
incomplete_data$Method_family[which(incomplete_data$`Paper name` %in% c("Ragazzola et al. 2012", "Ragazzola et al. 2013"))] = "Staining"
# Study from Semesi et al 2009
incomplete_data$Unit[which(incomplete_data$`Paper name` == "Semesi et al. 2009")] = "nmol/gWW/min"
incomplete_data$Method[which(incomplete_data$`Paper name` == "Semesi et al. 2009")] = "TA"
incomplete_data$Method_family[which(incomplete_data$`Paper name` == "Semesi et al. 2009")] = "TA anomaly"
# Study from Short et al 2014
incomplete_data$Unit[which(incomplete_data$`Paper name` == "Short et al. 2014")] = "mg/cm-2/year"
incomplete_data$Method[which(incomplete_data$`Paper name` == "Short et al. 2014")] = "BW"
incomplete_data$Method_family[which(incomplete_data$`Paper name` == "Short et al. 2014")] = "BW or RGR"
incomplete_data$Species[which(incomplete_data$`Paper name` == "Short et al. 2014")] = "Hydrolithon sp."

# Std units
incomplete_data <- incomplete_data %>% 
  mutate(std_extension_corrected_unit = case_when(Unit == "mm_year" ~          Rate),
         std_surf_corrected_unit      = case_when(Unit == "mg/cm-2/year" ~     Rate,
                                                  Unit == "µmol CaCO3/cm2/h" ~ Rate * molar_mass_CaCO3 / mol_to_umol  * day_to_hour     * g_to_mg),
         std_biom_corrected_unit      = case_when(Unit == "nmol/gWW/min" ~     Rate                    / nmol_to_umol * hour_to_minute))
                                        
# Add Genus informations
incomplete_data$Genus   <- stringr::word(incomplete_data$Species, 1)
incomplete_data$Species <- paste(stringr::word(incomplete_data$Species, 1), stringr::word(incomplete_data$Species, 2), sep = "_")

# Add information to precedent dataset
incomplete_data_biom <- incomplete_data %>% drop_na(std_biom_corrected_unit) %>% 
  dplyr::select(., -c(std_surf_corrected_unit, std_extension_corrected_unit))
incomplete_data_surf <- incomplete_data %>% drop_na(std_surf_corrected_unit) %>% 
  dplyr::select(., -c(std_biom_corrected_unit, std_extension_corrected_unit))
incomplete_data_ext <- incomplete_data %>% drop_na(std_extension_corrected_unit) %>% 
  dplyr::select(., -c(std_surf_corrected_unit, std_biom_corrected_unit))

colnames(CCA_biom) <- colnames(incomplete_data_biom) ; CCA_biom_tot <- rbind(CCA_biom, incomplete_data_biom)
colnames(CCA_surf) <- colnames(incomplete_data_surf) ; CCA_surf_tot <- rbind(CCA_surf, incomplete_data_surf)
colnames(CCA_non_std) <- colnames(incomplete_data_ext) ; CCA_ext_tot <- rbind(CCA_non_std, incomplete_data_ext)

CCA_biom_tot <- CCA_biom_tot %>% dplyr::filter(std_biom_corrected_unit >= 0)

#### Build the global dataset ----
# Biomass
CCA_biom_tot <- CCA_biom_tot %>% dplyr::select(., c(`Paper name`, Genus, Climate, Temperature, Method_family, 
                                                    std_biom_corrected_unit,Error, Unit)) %>% 
  mutate(Standardization = rep("Biomass – umol_g_h", length(CCA_biom_tot$Method_family))) %>% 
  rename(., Rate_std = std_biom_corrected_unit)
# Surface
CCA_surf_tot <- CCA_surf_tot %>% dplyr::select(., c(`Paper name`, Genus, Climate, Temperature, Method_family, 
                                                    std_surf_corrected_unit,Error, Unit)) %>% 
  mutate(Standardization = rep("Surface – mg_cm2_day", length(CCA_surf_tot$Method_family))) %>% 
  rename(., Rate_std = std_surf_corrected_unit)
# Extension
CCA_ext_tot  <- CCA_ext_tot %>% dplyr::select(., c(`Paper name`, Genus, Climate, Temperature, Method_family, 
                                                   std_extension_corrected_unit,Error, Unit)) %>% 
  mutate(Standardization = rep("Extension – mm_year", length(CCA_ext_tot$Method_family))) %>% 
  rename(., Rate_std = std_extension_corrected_unit)

## Merge the 3 datasets
CCA_Growth <- rbind(CCA_biom_tot, CCA_surf_tot, CCA_ext_tot)
unique(CCA_Growth$`Paper name`) # 55 studies used
unique(CCA_Growth$Genus)        # 22 Genera

# Define the error as well
CCA_Growth <- CCA_Growth %>% mutate(., Error = as.numeric(Error)) %>% 
  mutate(std_error = case_when(Unit == "nmoles g-1 dry weight min-l" ~ Error / nmol_to_umol * hour_to_minute,
                               Unit == "mgCaCO3 g-1 d-1"             ~ Error / molar_mass_CaCO3 / g_to_mg / day_to_hour * mol_to_umol,
                               Unit == "mg g–1 d–1)"                 ~ Error / molar_mass_CaCO3 / g_to_mg / day_to_hour * mol_to_umol,
                               Unit == "umol_g_h"                    ~ Error,
                               Unit == "umol/g(dwt)/min"             ~ Error * hour_to_minute,
                               Unit == "mg CaCO3 mg-1 day-1."        ~ Error / molar_mass_CaCO3 / day_to_hour * mol_to_umol,
                               Unit == "g per g (day-1)"             ~ Error / molar_mass_CaCO3 / day_to_hour * mol_to_umol,
                               Unit == "µmolC/g day"                 ~ Error / day_to_hour,
                               Unit == "µmol CaCO3/g/day"            ~ Error / day_to_hour,
                               Unit == "mgCaCO3 g-1 DW h-1"          ~ Error / molar_mass_CaCO3 / g_to_mg * mol_to_umol,
                               Unit == "g CaCO3/cm3/year"            ~ Error * calcite_density / molar_mass_CaCO3 * mol_to_umol / 
                                                                               year_to_day / day_to_hour,
                               Unit == "nmol/gWW/min"                ~ Error / nmol_to_umol * hour_to_minute,
                               Unit == "mg_cm2_day"                  ~ Error,
                               Unit == "umol CaCO3 cm-2 h-1"         ~ Error * molar_mass_CaCO3 / mol_to_umol * day_to_hour * g_to_mg,
                               Unit == "mmol/cm**2/day"              ~ Error * molar_mass_CaCO3,
                               Unit == "µmol CaCO3/cm2/h"            ~ Error * molar_mass_CaCO3 / mol_to_umol * day_to_hour * g_to_mg,
                               Unit == "mg/cm-2/year"                ~ Error / year_to_day,
                               Unit == "µm/day"                      ~ Error * year_to_day / mm_to_um,
                               Unit == "mm/day"                      ~ Error * year_to_day,
                               Unit == "cm month-1"                  ~ Error * year_to_month * cm_to_mm,
                               Unit == "mm d-1"                      ~ Error * year_to_day,
                               Unit == "mm_year"                     ~ Error ,
                               Unit == "mm/4 months"                 ~ Error * 3)) %>% 
  dplyr::select(`Paper name`, Genus, Climate, Temperature, Method_family, Rate_std, std_error, Standardization)



