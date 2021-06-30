library(here)
library(tidyverse)

source(here("R/LICAR_functions.R"))

res <- get_rel_abundance(lipid_species = "SM 32:1", lipid_fragment_class = "SM", product_mz = "184", isotope_no = 3)
res <- get_rel_abundance(lipid_species = c("Cer d18:1/24:0", "Cer d18:1/24:1"), lipid_fragment_class = "Cer", product_mz = "184", isotope_no = 3)
res <- get_rel_abundance(lipid_species = c("Cer d18:1/24:0"), lipid_fragment_class = "Cer", product_mz = "184", isotope_no = 3)

