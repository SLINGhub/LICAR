library(here)
library(tidyverse)

source(here("R/LICAR_functions.R"))


d_input <- read_csv(here("data/test/PC26_0_SM30_1.csv"), trim_ws = TRUE,col_names = TRUE) |> 
  column_to_rownames("filename")

res <- isoCorrect(d_input, lipidClass = "SM", lipidGroup = "Head Group")

