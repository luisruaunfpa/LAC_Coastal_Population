# UNFPA LACRO - LAC Coastal Population Assessment #
# Luis de la Rua - July 2022 #

#Objective: merge all output tables from LAC_Coast_Pop_byCountry script and analyze results.

library(tidyverse)
library(readxl)
library(dplyr)
library(ggplot2)
library(plotly)
require(scales)

wd <- "G:/.shortcut-targets-by-id/1i7Dq5-HfwQYp5ApFZZAw3pZb8GM0S76S/UNFPA Luis de la Rua/Coastal_Population_LAC"
setwd(wd)

# bring all tables in result folder
files <- list.files("result", pattern = "*.csv$",full.names = T)

# combine all tables
merged <- files %>%  map_df(read_csv)



# describe
str(merged)

# make field name more consistent
# merged<-rename(merged,elev05=elev5)

write_csv(merged,"output/merged.csv")
