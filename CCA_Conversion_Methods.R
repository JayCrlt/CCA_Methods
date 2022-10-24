#### Set up environment ----
rm(list = ls()) ; options(cores = parallel::detectCores())
library(brms); library(Matrix); library(tidyverse) ; library(readxl) ; library(patchwork) ; library(tidybayes)
library(viridis) ; library(viridisLite) ; library(ggridges) ; library(hrbrthemes) ; library(posterior)

# Useful functions
'%notin%' <- function(x,y)!('%in%'(x,y))

# Datasets
CCA_CC_SC      <- read_excel("Data/CCA_CC_SC.xlsx")
LTER           <- read_excel("Data/MCR_LTER.xlsx") 
First_dataset  <- read_excel("Data/Kornder_MA.xlsx", sheet = "Sheet2") %>% 
                    distinct(Pubyear, Genus, Species, conX) %>% drop_na()
Second_dataset <- read_excel("Data/Kornder_MA.xlsx", sheet = "Sheet3") %>% 
                    distinct(Pubyear, Genus, Species, conX)
data_Niklas    <- rbind(First_dataset, Second_dataset) %>% 
                    mutate(., conX = conX/1000*365.25) 
Halimeda       <- read_excel("Data/Halimeda_Rates.xlsx")

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

# Clean Env.
draft_dataset = list(CCA_biom, CCA_non_std, CCA_surf) ; raw_data = CCA_CC_SC
rm(Biom_stat, CCA_biom, CCA_non_std, CCA_surf, incomplete_data, incomplete_data_biom, incomplete_data_ext, incomplete_data_surf,
   nonstdstat, Surfstat, CCA_biom_tot, CCA_ext_tot, CCA_surf_tot, CCA_CC_SC)

# Add Halimeda to the study
colnames(Halimeda) = colnames(CCA_Growth)
CCA_Growth <- rbind(CCA_Growth, Halimeda)

# Vizualize
CCA_Growth$Climate[which(CCA_Growth$Climate == "Cool Temperate")] = "Cool temperate"
CCA_Growth$Genus[which(CCA_Growth$`Paper name` %in% c("Ragazzola et al. 2012", "Ragazzola et al. 2013"))] = "Lithothamnion"
col_climate = c("#c23728", "#e1a692", "#a7d5ed", "#1984c5")
Data_viz <- CCA_Growth %>% group_split(Standardization) 

# Chris spotted errors
Data_viz[[1]]$Climate[which(Data_viz[[1]]$Genus == "Arthrocardia")] = "Cool temperate"

# Add Articulate CCAs
data_articulate = data.frame(Genus = c("Bosiella", "Corallina", "Calliarthron", "Amphiroa", "Arthrocardia", 
                                       "Ellisolandia", "Jania"),
                             Articulate = rep("Yes", 7))

# 2 studies far above than others...
GR_viz = vector("list", 8)
GR_viz[[1]] <- Data_viz[[1]] %>% dplyr::filter(., `Paper name` %notin% c("Graba-Landry et al. 2018", "Johnson et al 2019"),
                                               Genus != "Halimeda") %>%
  dplyr::filter(., Genus %notin% data_articulate$Genus) %>% 
  ggplot(aes(x = Genus, shape = Method_family, y = Rate_std, color = Climate)) + 
  geom_linerange(aes(ymin = Rate_std - std_error, ymax = Rate_std + std_error, color = Climate), 
                 position = position_jitter(seed = 123, width = 0.3)) + theme_bw() +
  geom_point(aes(fill = Climate, color = Climate), position = position_jitter(seed = 123, width = 0.3), size = 2) +
  geom_point(aes(fill = Climate), color = "black", position = position_jitter(seed = 123, width = 0.3), size = 2, show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_shape_manual(name = "Methods", values=c(21, 22, 23, 24)) +
  scale_fill_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_color_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_x_discrete(name = "") + 
  scale_y_continuous(name = expression("Calcification rate (µmol."*g^-1*".h"^-1*")"), limits = c(0,30))

GR_viz[[2]] <- Data_viz[[1]] %>% dplyr::filter(., `Paper name` %notin% c("Graba-Landry et al. 2018", "Johnson et al 2019"),
                                               Genus != "Halimeda") %>%
  dplyr::filter(., Genus %in% data_articulate$Genus) %>% 
  ggplot(aes(x = Genus, shape = Method_family, y = Rate_std, color = Climate)) + 
  geom_linerange(aes(ymin = Rate_std - std_error, ymax = Rate_std + std_error, color = Climate), 
                 position = position_jitter(seed = 123, width = 0.3)) + theme_bw() +
  geom_point(aes(fill = Climate, color = Climate), position = position_jitter(seed = 123, width = 0.3), size = 2) +
  geom_point(aes(fill = Climate), color = "black", position = position_jitter(seed = 123, width = 0.3), size = 2, show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_shape_manual(name = "Methods", values=c(21, 22, 23, 24)) +
  scale_fill_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_color_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_x_discrete(name = "") + scale_y_continuous(name = "") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

GR_viz[[3]] <- Data_viz[[1]] %>% dplyr::filter(., `Paper name` %notin% c("Graba-Landry et al. 2018", "Johnson et al 2019"),
                                               Genus == "Halimeda") %>%
  ggplot(aes(x = Genus, shape = Method_family, y = Rate_std, color = Climate)) + 
  geom_linerange(aes(ymin = Rate_std - std_error, ymax = Rate_std + std_error, color = Climate), 
                 position = position_jitter(seed = 123, width = 0.3)) + theme_bw() +
  geom_point(aes(fill = Climate, color = Climate), position = position_jitter(seed = 123, width = 0.3), size = 2) +
  geom_point(aes(fill = Climate), color = "black", position = position_jitter(seed = 123, width = 0.3), size = 2, show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_shape_manual(name = "Methods", values=c(21, 22, 23, 24)) +
  scale_fill_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_color_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_x_discrete(name = "") + scale_y_continuous(name = "") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

GR_viz[[4]] <- Data_viz[[2]] %>% dplyr::filter(., Genus %notin% c("Halimeda")) %>% 
  dplyr::filter(., Genus %notin% data_articulate$Genus) %>% 
  ggplot(aes(x = Genus, shape = Method_family, y = Rate_std, color = Climate)) + 
  geom_linerange(aes(ymin = Rate_std - std_error, ymax = Rate_std + std_error, color = Climate), 
                 position = position_jitter(seed = 123, width = 0.3)) + theme_bw() +
  geom_point(aes(fill = Climate, color = Climate), position = position_jitter(seed = 123, width = 0.3), size = 2) +
  geom_point(aes(fill = Climate), color = "black", position = position_jitter(seed = 123, width = 0.3), size = 2, show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_shape_manual(name = "Methods", values=c(21, 22, 23, 24)) +
  scale_fill_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_color_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_x_discrete(name = "") + scale_y_continuous(name = expression("Linear extension (mm."*yr^-1*")")) 

GR_viz[[5]] <- Data_viz[[2]] %>% dplyr::filter(., Genus %notin% c("Halimeda")) %>% 
  dplyr::filter(., Genus %in% data_articulate$Genus) %>% 
  ggplot(aes(x = Genus, shape = Method_family, y = Rate_std, color = Climate)) + 
  geom_linerange(aes(ymin = Rate_std - std_error, ymax = Rate_std + std_error, color = Climate), 
                 position = position_jitter(seed = 123, width = 0.3)) + theme_bw() +
  geom_point(aes(fill = Climate, color = Climate), position = position_jitter(seed = 123, width = 0.3), size = 2) +
  geom_point(aes(fill = Climate), color = "black", position = position_jitter(seed = 123, width = 0.3), size = 2, show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_shape_manual(name = "Methods", values=c(21, 22, 23, 24)) +
  scale_fill_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_color_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_x_discrete(name = "") + scale_y_continuous(name = "") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

GR_viz[[6]] <- Data_viz[[3]] %>% dplyr::filter(., Genus %notin% c("Halimeda")) %>% 
  dplyr::filter(., Genus %notin% data_articulate$Genus) %>% 
  ggplot(aes(x = Genus, shape = Method_family, y = Rate_std, color = Climate)) + 
  geom_linerange(aes(ymin = Rate_std - std_error, ymax = Rate_std + std_error, color = Climate), 
                 position = position_jitter(seed = 123, width = 0.3)) + theme_bw() +
  geom_point(aes(fill = Climate, color = Climate), position = position_jitter(seed = 123, width = 0.3), size = 2) +
  geom_point(aes(fill = Climate), color = "black", position = position_jitter(seed = 123, width = 0.3), size = 2, show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_shape_manual(name = "Methods", values=c(21, 22, 23, 24)) +
  scale_fill_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_color_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_x_discrete(name = "") + 
  scale_y_continuous(name = expression("Calcification rate (mg."*cm^-2*".day"^-1*")"), 
                     limits = c(0,4.1), breaks = seq(0,4,1))

GR_viz[[7]] <- Data_viz[[3]] %>% dplyr::filter(., Genus %notin% c("Halimeda")) %>% 
  dplyr::filter(., Genus %in% data_articulate$Genus) %>% 
  ggplot(aes(x = Genus, shape = Method_family, y = Rate_std, color = Climate)) + 
  geom_linerange(aes(ymin = Rate_std - std_error, ymax = Rate_std + std_error, color = Climate), 
                 position = position_jitter(seed = 123, width = 0.3)) + theme_bw() +
  geom_point(aes(fill = Climate, color = Climate), position = position_jitter(seed = 123, width = 0.3), size = 2) +
  geom_point(aes(fill = Climate), color = "black", position = position_jitter(seed = 123, width = 0.3), size = 2, show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_shape_manual(name = "Methods", values=c(21, 22, 23, 24)) +
  scale_fill_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_color_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_x_discrete(name = "") + 
  scale_y_continuous(name = "", 
                     limits = c(0,4.1), breaks = seq(0,4,1)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

GR_viz[[8]] <- Data_viz[[3]] %>% dplyr::filter(., Genus == "Halimeda") %>% 
  ggplot(aes(x = Genus, shape = Method_family, y = Rate_std, color = Climate)) + 
  geom_linerange(aes(ymin = Rate_std - std_error, ymax = Rate_std + std_error, color = Climate), 
                 position = position_jitter(seed = 123, width = 0.3)) + theme_bw() +
  geom_point(aes(fill = Climate, color = Climate), position = position_jitter(seed = 123, width = 0.3), size = 2) +
  geom_point(aes(fill = Climate), color = "black", position = position_jitter(seed = 123, width = 0.3), size = 2, show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_shape_manual(name = "Methods", values=c(21, 22, 23, 24)) +
  scale_fill_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_color_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) +
  scale_x_discrete(name = "") + 
  scale_y_continuous(name = "", 
                     limits = c(0,4.1), breaks = seq(0,4,1)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# So ? What's the best technique ?
# No enough values for CT scan
Data_viz_1 = Data_viz[[1]] %>% dplyr::filter(., Method_family != "X-ray CT Scan", 
                                             `Paper name` %notin% c("Graba-Landry et al. 2018", "Johnson et al 2019"),
                                             Genus %notin% c("Halimeda"))
Data_viz_2 = Data_viz[[2]] %>% dplyr::filter(., Genus %notin% c("Halimeda"))
Data_viz_3 = Data_viz[[3]] %>% dplyr::filter(., Genus %notin% c("Halimeda"))

fit_biom      <- brms::brm(Rate_std ~ Climate + Genus + (1 | Method_family), data = Data_viz_1, family = gaussian(), 
                           warmup = 3000, iter = 5000, control = list(max_treedepth = 10, adapt_delta = 0.99), core = 4)
fit_extension <- brms::brm(Rate_std ~ Climate + Genus + (1 | Method_family), data = Data_viz_2, family = gaussian(), 
                           warmup = 3000, iter = 5000, control = list(max_treedepth = 10, adapt_delta = 0.99), core = 4)
fit_surface   <- brms::brm(Rate_std ~ Climate + Genus + (1 | Method_family), data = Data_viz_3, family = gaussian(), 
                           warmup = 3000, iter = 5000, control = list(max_treedepth = 11, adapt_delta = 0.99), core = 4)

summary_extension = fit_extension %>% spread_draws(r_Method_family[condition,]) %>% summarise_draws()
summary_surface   = fit_surface %>% spread_draws(r_Method_family[condition,]) %>% summarise_draws()
summary_biomass   = fit_biom %>% spread_draws(r_Method_family[condition,]) %>% summarise_draws()

# So depending about what you're looking at
# If you're looking at growth, then, the staining technic looks more promising as the deviation is lower, with a similar amount of results
# If you're looking to calcification, it's BW <= TA <= X-ray CT-Scan < Isotops but TA always overestimate

summary           = rbind(summary_biomass, summary_extension, summary_surface)
summary$Methods   = chartr(".", " ", summary$condition)
vector_method     = vector("list", 8) ; for (i in 1:8) {vector_method[[i]] <- rnorm(100000, summary$mean[i], summary$sd[i])}
data_methods = data.frame(Methods = rep(c("BW or RGR", "Isotopes", "TA anomaly", "BW or RGR", "Staining", "BW or RGR", "TA anomaly", "X-ray CT Scan"), each = 100000),
                          Random_effect = abind::abind(vector_method), 
                          dataset = c(rep("Biomass", 300000), rep("Extension", 200000), rep("Surface", 300000))) 

data_methods = complete(data_methods, Methods, dataset, fill = list(Random_effect = NA)) %>% group_split(dataset)
summary_list = list() ; summary_list[[1]] = summary[1:3,] ; summary_list[[2]] = summary[4:5,] ; summary_list[[3]] = summary[6:8,]

colors = c("#ff595e", "#ffca3a", "#8ac926", "#1982c4", "#6a4c93")
color_list = list() ; color_list[[1]] = colors[c(1:2,4)] ; color_list[[2]] = colors[c(1,3)] ; color_list[[3]] = colors[c(1,4:5)]

Plot = list() ; for (i in 1:3) { Plot[[i]] <- data_methods[[i]] %>%
    ggplot(aes(y = Methods, x = Random_effect, fill = Methods, color = Methods)) + 
    geom_density_ridges(alpha=0.6, bandwidth=4, fill = NA, color = NA) + 
    theme_ipsum() +
    geom_segment(data = summary_list[[i]], aes(x = mean-sd, xend = mean+sd, y = Methods , yend = Methods, color = Methods), size = 2) +
    geom_point(data = summary_list[[i]], aes(x = mean, y = Methods, fill = Methods), size = 4, shape = 21, color = "black") +
    scale_x_continuous(limits = c(-5,10)) +
    scale_y_discrete(name = "") +
    scale_fill_manual(values = color_list[[i]]) + scale_color_manual(values = color_list[[i]]) + 
    theme(legend.position="none", panel.spacing = unit(0.1, "lines"), strip.text.x = element_text(size = 8)) + coord_flip() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) }

Figure_S1A <- GR_viz[[1]] + ggtitle("(a)") + theme(plot.title = element_text(face = "bold")) +
             GR_viz[[2]] + ggtitle("   ") + theme(plot.title = element_text(face = "bold")) +
             GR_viz[[3]] + ggtitle("   ") + theme(plot.title = element_text(face = "bold")) +
             plot_spacer() +
             GR_viz[[4]] + ggtitle("(c)") + theme(plot.title = element_text(face = "bold")) +
             GR_viz[[5]] + ggtitle("   ") + theme(plot.title = element_text(face = "bold")) +
             plot_spacer() +
             GR_viz[[6]] + ggtitle("(e)") + theme(plot.title = element_text(face = "bold")) +
             GR_viz[[7]] + ggtitle("   ") + theme(plot.title = element_text(face = "bold")) +
             GR_viz[[8]] + ggtitle("   ") + theme(plot.title = element_text(face = "bold")) +
  plot_layout(guides = "collect", nrow = 1, widths = c(9,4,1,2,9,3,3,9,1,1)) &
  scale_shape_manual(name = "Methods", limits = unique(CCA_Growth$Method_family), values = c(22, 21, 23, 24, 25)) &
  scale_color_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) &
  theme(plot.margin = unit(rep(0.1,4),"cm"))

Figure_S1B <- (Plot[[1]] + ggtitle("(b)") + theme(plot.title = element_text(face = "bold", size = 15)) +
              Plot[[2]] + ggtitle("(d)") + theme(plot.title = element_text(face = "bold", size = 15)) +
              Plot[[3]] + ggtitle("(f)") + theme(plot.title = element_text(face = "bold", size = 15))) 

Figure_S1 <- Figure_1A / Figure_1B

Figure_1A <-  GR_viz[[1]] + ggtitle("(a)") + theme(plot.title = element_text(face = "bold")) +
              GR_viz[[2]] + ggtitle("   ") + theme(plot.title = element_text(face = "bold")) +
              GR_viz[[3]] + ggtitle("   ") + theme(plot.title = element_text(face = "bold")) +
              plot_spacer() +
              GR_viz[[6]] + ggtitle("(c)") + theme(plot.title = element_text(face = "bold")) +
              GR_viz[[7]] + ggtitle("   ") + theme(plot.title = element_text(face = "bold")) +
              GR_viz[[8]] + ggtitle("   ") + theme(plot.title = element_text(face = "bold")) +
              plot_layout(guides = "collect", nrow = 1, widths = c(9,4,1,2,9,1,1)) &
              scale_shape_manual(name = "Methods", limits = c("Isotopes", "BW or RGR", "TA anomaly", "X-ray CT Scan"), 
                                 values = c(22, 21, 23, 24)) &
              scale_color_manual(values = col_climate, limits = c("Tropical", "Warm temperate", "Cool temperate", "Polar")) &
              theme(plot.margin = unit(rep(0.1,4),"cm"))

Plot = list() ; for (i in 1:3) { Plot[[i]] <- data_methods[[i]] %>% dplyr::filter(Methods != "Staining") %>% 
  ggplot(aes(y = Methods, x = Random_effect, fill = Methods, color = Methods)) + 
  geom_density_ridges(alpha=0.6, bandwidth=4, fill = NA, color = NA) + 
  theme_ipsum() +
  geom_segment(data = summary_list[[i]], aes(x = mean-sd, xend = mean+sd, y = Methods , yend = Methods, color = Methods), size = 2) +
  geom_point(data = summary_list[[i]], aes(x = mean, y = Methods, fill = Methods), size = 4, shape = 21, color = "black") +
  scale_x_continuous(limits = c(-5,5)) +
  scale_y_discrete(name = "") +
  scale_fill_manual(values = color_list[[i]]) + scale_color_manual(values = color_list[[i]]) + 
  theme(legend.position="none", panel.spacing = unit(0.1, "lines"), strip.text.x = element_text(size = 8)) + coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) }

Figure_1B <- (Plot[[1]] + ggtitle("(b)") + theme(plot.title = element_text(face = "bold", size = 15)) +
              Plot[[3]] + ggtitle("(d)") + theme(plot.title = element_text(face = "bold", size = 15))) 

### Making Figure 2 CCA vs Coral CaCO3 production rates ----
data_cca = Data_viz[[3]] %>% mutate(., Rate_std = Rate_std/1000*365.25, std_error = std_error/1000*365.25)

##### WORKING WITH QUANTILES
# CCA
Quantile_CCA_Data_all = data_cca %>% dplyr::filter(Climate %in% c("Warm temperate", "Tropical")) %>% 
  group_by(Genus) %>% summarise(.,         Q_0.01 = quantile(Rate_std, probs = 0.01),
                                           Q_0.25 = quantile(Rate_std, probs = 0.25),
                                           Q_0.50 = quantile(Rate_std, probs = 0.50),
                                           Q_0.75 = quantile(Rate_std, probs = 0.75),
                                           Q_0.99 = quantile(Rate_std, probs = 0.99))

Quantile_CCA_Data_avg = data.frame(Genus = "All CCA", 
                                   data_cca %>% dplyr::filter(Climate %in% c("Warm temperate", "Tropical")) %>% 
                                      dplyr::filter(Genus %notin% data_articulate$Genus) %>% summarise(., 
                                           Q_0.01 = quantile(Rate_std, probs = 0.01),
                                           Q_0.25 = quantile(Rate_std, probs = 0.25),
                                           Q_0.50 = quantile(Rate_std, probs = 0.50),
                                           Q_0.75 = quantile(Rate_std, probs = 0.75),
                                           Q_0.99 = quantile(Rate_std, probs = 0.99)))

# Corals
Quantile_Cor_Data_all = data_Niklas %>% group_by(Genus) %>% summarise(., 
                                           Q_0.01 = quantile(conX, probs = 0.01),
                                           Q_0.25 = quantile(conX, probs = 0.25),
                                           Q_0.50 = quantile(conX, probs = 0.50),
                                           Q_0.75 = quantile(conX, probs = 0.75),
                                           Q_0.99 = quantile(conX, probs = 0.99))

Quantile_Cor_Data_avg = data.frame(Genus = "All Corals", data_Niklas %>% summarise(., 
                                           Q_0.01 = quantile(conX, probs = 0.01),
                                           Q_0.25 = quantile(conX, probs = 0.25),
                                           Q_0.50 = quantile(conX, probs = 0.50),
                                           Q_0.75 = quantile(conX, probs = 0.75),
                                           Q_0.99 = quantile(conX, probs = 0.99)))

quantile_CCA_all = rbind(Quantile_CCA_Data_all, Quantile_CCA_Data_avg)
Quantile_Cor_all = rbind(Quantile_Cor_Data_all, Quantile_Cor_Data_avg)

summary_all <- rbind(Quantile_Cor_all %>% dplyr::filter(., Genus == "All Corals"),
                     quantile_CCA_all %>% dplyr::filter(., Genus == "All CCA"))

Figure_3A <- summary_all %>%
  ggplot(aes(x = Genus, y = Q_0.50, col = Genus), show.legend = F) + 
  geom_linerange(aes(ymin = Q_0.01, ymax = Q_0.99), size = 1, show.legend = F) +
  geom_linerange(aes(ymin = Q_0.25, ymax = Q_0.75), size = 3, show.legend = F) + theme_bw() +
  geom_point(size = 5, show.legend = F, aes(shape = Genus, fill = Genus), col = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(name = "") + 
  scale_shape_manual(name = "", limits = c("All Corals", "All CCA"), values = c(23, 21)) +
  scale_color_manual(values = c("violetred", "orange")) + scale_fill_manual(values = c("violetred", "orange")) +
  scale_y_continuous(name = expression("Calcification rate (g."*cm^-2*".yr"^-1*")"), 
                     limits = c(0,1), breaks = seq(0,1,0.1))

Figure_3Ba <- quantile_CCA_all %>% dplyr::filter(., Genus != "All CCA", Genus %notin% c("Amphiroa", "Halimeda")) %>% 
  ggplot(aes(x = Genus, y = Q_0.50)) + 
  geom_linerange(aes(ymin = Q_0.01, ymax = Q_0.99), size = 1, col = "violetred") +
  geom_linerange(aes(ymin = Q_0.25, ymax = Q_0.75), size = 3, col = "violetred") + theme_bw() +
  geom_point(size = 5, show.legend = F, shape = 21, fill = "violetred") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(name = "") + 
  scale_y_continuous(name = "", limits = c(0,1), breaks = seq(0,1,0.1)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

Figure_3Bb <- quantile_CCA_all %>% dplyr::filter(., Genus == "Amphiroa") %>% 
  ggplot(aes(x = Genus, y = Q_0.50)) + 
  geom_linerange(aes(ymin = Q_0.01, ymax = Q_0.99), size = 1, col = "violetred") +
  geom_linerange(aes(ymin = Q_0.25, ymax = Q_0.75), size = 3, col = "violetred") + theme_bw() +
  geom_point(size = 5, show.legend = F, shape = 21, fill = "violetred") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(name = "") + 
  scale_y_continuous(name = "", limits = c(0,1), breaks = seq(0,1,0.1)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

Figure_3Bc <- quantile_CCA_all %>% dplyr::filter(., Genus == "Halimeda") %>% 
  ggplot(aes(x = Genus, y = Q_0.50)) + 
  geom_linerange(aes(ymin = Q_0.01, ymax = Q_0.99), size = 1, col = "violetred") +
  geom_linerange(aes(ymin = Q_0.25, ymax = Q_0.75), size = 3, col = "violetred") + theme_bw() +
  geom_point(size = 5, show.legend = F, shape = 21, fill = "violetred") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(name = "") + 
  scale_y_continuous(name = "", limits = c(0,1), breaks = seq(0,1,0.1)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

Figure_3C <- Quantile_Cor_all %>% dplyr::filter(., Genus != "All Corals") %>% 
  ggplot(aes(x = Genus, y = Q_0.50)) + 
  geom_linerange(aes(ymin = Q_0.01, ymax = Q_0.99), size = 1, col = "orange") +
  geom_linerange(aes(ymin = Q_0.25, ymax = Q_0.75), size = 3, col = "orange") + theme_bw() +
  geom_point(size = 5, show.legend = F, shape = 23, fill = "orange") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(name = "") + 
  scale_y_continuous(name = "", limits = c(0,1), breaks = seq(0,1,0.1)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

Figure_3 <- Figure_3A + ggtitle("(a)") + theme(plot.title = element_text(face = "bold")) +
  plot_spacer() +
  Figure_3Ba + ggtitle("(b)") + theme(plot.title = element_text(face = "bold")) +
  Figure_3Bb + ggtitle("") + theme(plot.title = element_text(face = "bold")) +
  Figure_3Bc + ggtitle("") + theme(plot.title = element_text(face = "bold")) +
  plot_spacer() +
  Figure_3C + ggtitle("(c)") + theme(plot.title = element_text(face = "bold")) +
  plot_layout(widths = c(2, 1, 8, 1, 1, 1, 10))  & 
  theme(plot.tag = element_text(face = "bold"), text = element_text(size = 20), plot.margin = unit(rep(0.1,4),"cm"))

# Travis and Ben
# quantile(data_cca$Rate_std[-c(1:3)], probs = c(0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99))
quantile(data_cca$Rate_std, probs = c(0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99))
quantile(data_Niklas$conX , probs = c(0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99))

#### Figure Case Study Moorea ----
# Working for Figure 2 Mo'orea Case Study
MOZ_data = GR_viz[[3]]$data[c(7:19, 27:29),] %>% summarise(., GR = mean(Rate_std)*1000/365.25, sd = sd(Rate_std)*1000/365.25)
LTER_cca <- LTER %>% dplyr::filter(., Taxonomy_Substrate_Functional_Group == "Crustose Corallines")
corals_all = data_Niklas %>% group_by(Genus) %>% summarise(mean = mean(conX), sd = sd(conX))
corals_avg = data.frame(Genus = "All Corals", data_Niklas %>% summarise(mean = mean(conX), sd = sd(conX)))
corals_all = rbind(corals_all, corals_avg)
# LTER Dataset
# From Carlot et al. (2022)
SC = data.frame(Year = c(2005, 2008:2016), SC = c(3.07, 2.73, 2.46, 1.75, 1.58, 2.21, 1.99, 2.52, 2.30, 3.87))

## Datasets management
# CCA
LTER_cca <- LTER %>% dplyr::filter(., Taxonomy_Substrate_Functional_Group == "Crustose Corallines") %>% 
  dplyr::filter(., Year %in% SC$Year)
LTER_Summary_cca <- LTER_cca %>% group_by(Year, Site, Habitat, Transect, Quadrat) %>% 
  summarise(CCA_Cover = mean(Percent_Cover)) %>% group_by(Year, Site, Habitat, Transect) %>% 
  summarise(CCA_Cover = mean(CCA_Cover)) %>% group_by(Year, Site, Habitat) %>% 
  summarise(sd_CCA_Cover = sd(CCA_Cover), CCA_Cover = mean(CCA_Cover)) %>% group_by(Year, Site) %>% 
  summarise(CCA_Cover = mean(CCA_Cover), sd_CCA_Cover = mean(sd_CCA_Cover)) %>% 
  mutate(., CR = CCA_Cover/100 * 1.23)
LTER_CCA_avg = LTER_Summary_cca %>% group_by(Year) %>% 
  summarise(Cover = mean(CCA_Cover), sd_cover = sd(CCA_Cover), sd_CR = sd(CR), CR = mean(CR)) %>% 
  mutate(., Cover = as.numeric(Cover), sd_cover = as.numeric(sd_cover), CR = as.numeric(CR),
         sd_CR = as.numeric(sd_CR), Year = as.numeric(Year))

# Coral
LTER_coral <- LTER %>% dplyr::filter(., Taxonomy_Substrate_Functional_Group == "Coral") %>% 
  dplyr::filter(., Year %in% SC$Year)
LTER_Summary_coral <- LTER_coral %>% group_by(Year, Site, Habitat, Transect, Quadrat) %>% 
  summarise(CCA_Cover = mean(Percent_Cover)) %>% group_by(Year, Site, Habitat, Transect) %>% 
  summarise(CCA_Cover = mean(CCA_Cover)) %>% group_by(Year, Site, Habitat) %>% 
  summarise(sd_CCA_Cover = sd(CCA_Cover), CCA_Cover = mean(CCA_Cover)) %>% group_by(Year, Site) %>% 
  summarise(CCA_Cover = mean(CCA_Cover), sd_CCA_Cover = mean(sd_CCA_Cover)) %>% 
  mutate(., CR = CCA_Cover/100 * 3.09) %>% mutate(., Year = as.numeric(Year)) %>% 
  inner_join(., SC, by = "Year") %>% mutate(., CR_SC = CR * SC)
LTER_Coral_avg = LTER_Summary_coral %>% group_by(Year) %>% 
  summarise(sd_cover = sd(CCA_Cover), Cover = mean(CCA_Cover), sd_CR = sd(CR), CR = mean(CR), CR = mean(CR)) %>% 
  mutate(., Cover = as.numeric(Cover), sd_cover = as.numeric(sd_cover), CR = as.numeric(CR),
         sd_CR = as.numeric(sd_CR), Year = as.numeric(Year), CR_SC = CR * SC$SC, sd_CR_SC = sd_CR * SC$SC)

# Contribution
Contribution = data.frame(Year = LTER_Summary_cca$Year, 
                          Site = LTER_Summary_cca$Site,
                          Contribution = (LTER_Summary_cca$CR / (LTER_Summary_cca$CR + LTER_Summary_coral$CR_SC)) * 100)
Contribution_avg = Contribution %>% group_by(Year) %>% 
  summarise(Contribution_avg = mean(Contribution), sd = sd(Contribution)) %>% 
  mutate(., upr = Contribution_avg+sd, lwr = Contribution_avg-sd) %>% 
  mutate(., Year = as.numeric(Year), upr = as.numeric(upr), lwr = as.numeric(lwr))

## Figure Moorea Case study
Figure_2A <- ggplot(LTER_CCA_avg, aes(x = Year - 0.2, y = Cover)) + 
  geom_line(col = "pink", linetype = "dashed", size = .5) +
  geom_linerange(aes(ymin = Cover - sd_cover, ymax = Cover + sd_cover), col = "violetred") + theme_bw() +
  geom_line(data = LTER_Coral_avg, aes(x = Year + 0.2, y = Cover), col = "gold", linetype = "dashed", size = .5) +
  geom_linerange(data = LTER_Coral_avg, aes(x = Year + 0.2, ymin = Cover - sd_cover, ymax = Cover + sd_cover), col = "orange") +
  geom_point(size = 4, show.legend = F, shape = 21, fill = "violetred") +
  geom_point(data = LTER_Coral_avg, aes(x = Year + 0.2, y = Cover), size = 4, show.legend = F, shape = 23, fill = "orange") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(name = "", breaks = seq(2005, 2020,1)) +
  scale_y_continuous(name = "Cover (%)", limits = c(0,50), breaks = seq(0,50,10)) 

Figure_2B <- ggplot(LTER_CCA_avg, aes(x = Year - 0.2, y = CR)) + 
  geom_line(col = "pink", linetype = "dashed", size = .5) +
  geom_linerange(aes(ymin = CR - sd_CR, ymax = CR + sd_CR), col = "violetred") + theme_bw() +
  geom_line(data = LTER_Coral_avg, aes(x = Year + 0.2, y = CR_SC), col = "gold", linetype = "dashed", size = .5) +
  geom_linerange(data = LTER_Coral_avg, aes(x = Year + 0.2, ymin = CR_SC - sd_CR_SC, ymax = CR_SC + sd_CR_SC), col = "orange") +
  geom_point(size = 4, show.legend = F, shape = 21, fill = "violetred") +
  geom_point(data = LTER_Coral_avg, aes(x = Year + 0.2, y = CR_SC), size = 4, show.legend = F, shape = 23, fill = "orange") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(name = "", breaks = seq(2005, 2020,1)) + 
  scale_y_continuous(name = expression(atop(CaCO[3]~"production", paste("(kg."*m^-2*".yr"^-1*")"))), 
                     limits = c(0,6), breaks = seq(0,6,1)) 

Figure_2C <- ggplot(Contribution_avg, aes(x = as.numeric(Year), y = as.numeric(Contribution_avg))) + 
  geom_line(col = "pink", linetype = "dashed", size = .5) +
  geom_linerange(aes(ymin = as.numeric(Contribution_avg) - as.numeric(sd), 
                     ymax = as.numeric(Contribution_avg) + as.numeric(sd)), col = "violetred") + theme_bw() +
  geom_point(size = 4, show.legend = F, shape = 21, fill = "violetred") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(name = "", breaks = seq(2005,2020,1)) + 
  scale_y_continuous(name = expression(atop("CCA contribution to", paste("the"~CaCO[3]~"budget (%)"))), 
                     limits = c(0,100), breaks = seq(0,100,10)) 

Figure_2 <- Figure_2A + ggtitle("(a)") + theme(plot.title = element_text(face = "bold")) +
  Figure_2B + ggtitle("(b)") + theme(plot.title = element_text(face = "bold")) +
  Figure_2C + ggtitle("(c)") + theme(plot.title = element_text(face = "bold"))  & 
  theme(plot.tag = element_text(face = "bold"), text = element_text(size = 20), plot.margin = unit(rep(0.5,4),"cm")) &
  scale_shape_manual(name = "", values = c(21,23),                 limits = c("CCA", "Corals")) &
  scale_color_manual(name = "", values = c("violetred", "orange"), limits = c("CCA", "Corals")) 

ggsave(Figure_1A, filename = "Figures/Figure_5A.png", device = "png", height = 10, width = 25, units = "cm", dpi = 300)
ggsave(Figure_1B, filename = "Figures/Figure_5B.png", device = "png", height = 10, width = 25, units = "cm", dpi = 300)
ggsave(Figure_3, filename = "Figures/Figure_2.png", device = "png", height = 15, width = 35, units = "cm", dpi = 300)
ggsave(Figure_2, filename = "Figures/Figure_4.png", device = "png", height = 13, width = 35, units = "cm", dpi = 300)


