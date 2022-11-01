# UNFPA LACRO - LAC Coastal Population Assessment #
# Luis de la Rua - July 2022 #

#Objective: merge all output tables from LAC_Coast_Pop_byCountry script and analyze results.

library(tidyverse)
library(readxl)
library(dplyr)
library(ggplot2)
library(plotly)
require(scales)

wd <- "C:/Users/luisr/OneDrive - SPC/Documents/Personal/UNFPA/Coastal Population"
setwd(wd)
#combine all tables
merged<-dir("result",full.names = T) %>% map_df(read_csv)

#describe
str(merged)

# make field name more consistent
# merged<-rename(merged,elev05=elev5)

write_csv(merged,"output/merged.csv")

# #take a first look onto the dataset
# merged %>% group_by(country) %>% 
#   summarize(mean_per01 = mean(per_pop_b1),
#             mean_per05 = mean(per_pop_b5),
#             mean_per10 = mean(per_pop_b10))
# 
# #check for inconsistencies (pop05>pop10>pop20)
# merged %>% 
#   filter((pop_b1>pop_b5)|(pop_b5>pop_b10)|(pop_b1>pop_b10)|(per_pop_b1>1)|(per_pop_b5>1)|((per_pop_b10>1))) %>% arrange(country)
# 
# 
# ggplot(data=merged,mapping = aes(x=country,y=per_pop_b1))+geom_point(alpha = 0.6, color = "blue")
# ggplot(data=merged,mapping = aes(x=country,y=per_pop_b1))+geom_point(alpha = 0.1, aes(color = country))
# 
# 
# #ggplot(data=merged,mapping = aes(x=dem,y=elev05_per))+geom_point(alpha = 0.6, aes(color= country))
# 
# elev05_per_graph<-ggplot(data=merged,mapping = aes(x=country,y=elev05_per))+geom_point(alpha = 0.6, aes(color= dem))
# elev10_per_graph<-ggplot(data=merged,mapping = aes(x=country,y=elev10_per))+geom_point(alpha = 0.6, aes(color= dem))
# elev20_per_graph<-ggplot(data=merged,mapping = aes(x=country,y=elev20_per))+geom_point(alpha = 0.6, aes(color= dem))
# plot(elev20_per_graph)
# plot(elev10_per_graph)
# 
# #Population by country and threshold, and reshape it for the plot
# pop_alos<-merged %>% 
#   filter(dem=="ALOS") %>%
#   mutate(elev20_per=elev20_per-elev10_per) %>%  # substract elev20-elev10 for ggplot stack bars
#   mutate(elev20_per=elev20_per*100,elev10_per=elev10_per*100) %>% 
#   select(country, dem, elev10_per,elev20_per) %>% 
#   gather(key=elev,value = per,elev10_per,elev20_per)
# 
# alosgraph_1020_stk_bar<- 
#   ggplot(data=pop_alos,mapping = aes(x=country,y=per, width = 0.9,fill=elev)) + 
#   geom_bar(position = position_stack(reverse=T) ,stat="identity",color="black")
# 
# #Legend and colors
# alosgraph_1020_stk_bar<-alosgraph_1020_stk_bar +
#   scale_fill_manual(name="Elevation thresholds",
#                     values = c("elev10_per"="#a6bddb",
#                               "elev20_per"="#2b8cbe"),
#                     labels = c("elev10_per" = "LECZ-10",
#                                "elev20_per"="LECZ-20"))  
# #title
# alosgraph_1020_stk_bar<-alosgraph_1020_stk_bar + 
#   ggtitle("Proportion of population located in Low Elevation Coastal Zones")
# # Setting x and y labels
# alosgraph_1020_stk_bar<-alosgraph_1020_stk_bar + 
#   xlab("Countries") + ylab("Percentage of population (%)")
# alosgraph_1020_stk_bar
# 
# 
# #interactive need to rework this to improve titles
# stk_bar_interact<-ggplotly(alosgraph_1020_stk_bar,tooltip = "text")
# 
# #barchart showing populations
# 
# pop_alos2<-merged %>% 
#   filter(dem=="ALOS") %>%
#   #mutate(elev20=elev20-elev10) %>% 
#   select(country,elev10,elev20, dem) %>% 
#   gather(key=elev,value = pop,elev10,elev20)
# 
# leczpop_plot <-
#   ggplot(data=pop_alos2,mapping = aes(x=country,y=pop, width = 0.9,fill=elev)) + 
#   geom_bar(position = position_dodge2(reverse=F) ,stat="identity",color="black")
# #Legend and colors
# leczpop_plot<-leczpop_plot+
#   scale_fill_manual(name="Elevation thresholds",
#                     values = c("elev10"="#a6bddb",
#                                "elev20"="#2b8cbe"),
#                     labels = c("elev10" = "LECZ-10",
#                                "elev20"="LECZ-20"))  
# #title
# leczpop_plot<-leczpop_plot + 
#   ggtitle("Population located in Low Elevation Coastal Zones")
# # Setting x and y labels
# leczpop_plot<-leczpop_plot + 
#   xlab("Countries") + ylab("Population (pers.)")
# leczpop_plot<-leczpop_plot + 
# scale_y_continuous(labels = scales::comma)
# 
# leczpop_plot 
# 
# 
# leczpop_plot 
