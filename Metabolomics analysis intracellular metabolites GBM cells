---
title: "IGC GBM cells"
author: "Pauline Mencke"
date: "10/1/2020"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


# Loading all the needed libraries and setting a theme for when knitting
```{r Library loading, include=FALSE}

library(data.table)
library(janitor)
library(tidyxl)
library(reshape2)
library(pillar)
library(anchors)
library(tidyverse)
library(ggplot2)
library(dunn.test)
library(ggpubr) 
library(dplyr)
library(rstatix)
theme_set(theme_classic())
knitr::opts_chunk$set(echo = TRUE)

```

# Loading the data
```{r Loading of Data}


setwd('C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/R analysis metabolomics/Intracellular/')

#load data table containing MIDs

#remove nans in Excel with ctrl H 

TableImport = readxl::read_excel('C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/IGC metabolomics/Astrocytes and GBM cells/Analysis/raw data/20201001_Z3HV3_ZZZBH_4GFXW_J3VTI_PME3_PME7_PME8_PME9_IGC_MIDs.xlsx', range = "A1:DF297")

#remove all NAs (this removes also the disturbing metabolite names in the mass column)
Table = na.omit(TableImport)



```


# Prepare the tables to aqcuire the FinalTable
```{r Preparing tables}
#calculate the mean of each technical triplicate for each M of each metabolite (meaning all rows but only always for 3 columns)

#melt the table so that the data will be tidy (every experimental condition now has one value in a column)
Table_V0 <- melt(Table, id.vars=c("Metabolite_name", "Mass")) 
Table_V1 <- Table_V0%>% 
  rename(Full_Condition_Name = variable)

#separate metabolite name into metabolite name only 
Table_V2 <- Table_V1 %>% 
  separate(Metabolite_name, into = c("metabolite", "shit"), sep = "_")

#silence shit from metabolite name 
Table_V2$shit = NULL

#separate mass name into mass name only 
Table_V2 <- Table_V2 %>% 
  separate(Mass, into = c("crap", "Mass"), sep = " ")

#silence crap from mass 
Table_V2$crap = NULL



#separate full condiiton name so that the info on technical replicates and biological replicates can be extracted 
Table_V3 <- Table_V2 %>% 
  separate(Full_Condition_Name, into = c("V1", "V2", "V3", "V4", "V5", "V6", "V7"), sep = "_") %>% 
  #glue the name with the info needed   

  mutate(Celltype = glue::glue("{V3}"))
      


Table_GBM <- Table_V3[grep("^GBM", Table_V3$Celltype), ] %>% 
  mutate(Name = glue::glue("{V4}_{V5}_{V6}")) %>% 
  mutate(Condition = glue::glue("{V6}")) %>% 
  mutate(Cellline = glue::glue("{V4}")) %>% 
  mutate(Technical_Replicate = case_when(V7 == "1"~ 1,
                                         V7 == "2"~ 2,
                                         V7 == "3"~ 3)) %>%  
  dplyr::select(metabolite, Mass, Name, Cellline, Condition, Technical_Replicate, value)

write_csv(Table_GBM, "C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/IGC metabolomics/Astrocytes and GBM cells/Analysis/output files/Table_GBM.csv")

Table_GBM_means <- Table_GBM %>% 
  ungroup %>% 
  group_by(metabolite, Mass, Name) %>% 
  summarise("Mean_Technical_Replicates" = mean(value),
         "Stdev_Technical_Replicates" = sd(value))  
  #as we want to have only one value for each mean of the technical replicates, we want to shrink the table by using the distinct function 
  # distinct(metabolite, Mass, Name, Technical_Replicate, .keep_all = TRUE)



Table_LN229 <- Table_GBM_means[grep("^LN229", Table_GBM_means$Name), ]
Table_U251 <- Table_GBM_means[grep("^U251", Table_GBM_means$Name), ]
Table_U87 <- Table_GBM_means[grep("^U87", Table_GBM_means$Name), ]




write_csv(Table_LN229, "C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/IGC metabolomics/Astrocytes and GBM cells/Results/Table_LN229.csv")
write_csv(Table_U251, "C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/IGC metabolomics/Astrocytes and GBM cells/Results/Table_U251.csv")
write_csv(Table_U87, "C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/IGC metabolomics/Astrocytes and GBM cells/Results/Table_U87.csv")



```


# Automated image generator
```{r automated plots including stats, fig.height= 8, fig.width= 15}

#Define where you want to have your plots
DataDirectory <- "C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/IGC metabolomics/Astrocytes and GBM cells/Results/Barplots_GBM_cells/U251/"

# Define the graph function
# Check what you want to show, and change the y value accordingly
make_graph <- function(Table_U251, x, DataDirectory){
  #define the metabolite name and set a name for the file
  metabolite_name <- x  
  name <- paste(DataDirectory, metabolite_name, ".png", sep="")
 
  #open image creating part of the function, set resolution
  png(name, width = 720, height = 720, units = "px") 
  
 
  
  #Make the graph
  Graph <-
  ggplot(Data, aes(x = Mass, y = Mean_Technical_Replicates, col = Name))+
  geom_bar(aes(fill = Name, color = Name), colour = "black", stat="identity", position=position_dodge(0.8))+  
  geom_errorbar(aes(ymin = Ymin, ymax = Ymax, fill = Name), position=position_dodge(width=0.8), colour="black")+
  geom_point(position=position_dodge(width=0.8))+
  ggtitle(metabolite_name)+
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
  #Add statistics
  # stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.001, bracket.nudge.y = 0, dodge = 0.8, hide.ns = TRUE)
  
  #by printing, the graph is saved at the datadirectory   
  print(Graph)
  
  #the png module is closed
  dev.off()
  
}

# #Check what you want to show, mean_technical or mean_biological, and adjust the stat.test, weight, Ymin and Ymax accordingly
# get.stats<- function(Data, x){
#   metabolite_name <- x
#   
#   stat.test <- Data %>% 
#   group_by(Mass) %>%
#   dunn_test(Mean_Technical_Replicates ~ Name, p.adjust.method = "BH")
# 
#   #Here the xy position is added to the stat.test table. When we only want to show certain comparisons, we can fill these in here. However, since you have many different types of names, this can be a bit difficult to automate... can for sure be done on a subset of plots though.
#   stat.test <- stat.test %>%
#   add_xy_position(x = "Mass", fun = "max",  step.increase = 0.03, dodge = 0.8)
#   
# }

#Run the for loop and it should work!
#You should change the Ymin and Ymax if you change what you want to show (now it is Mean_Biological_Replicates and its stdev)
for(metabolite_name in unique(Table_U251$metabolite)){
  Data <- Table_U251 %>% 
    filter(metabolite == metabolite_name) %>% 
    mutate(Ymin = (Mean_Technical_Replicates-Stdev_Technical_Replicates),
         Ymax = (Mean_Technical_Replicates+Stdev_Technical_Replicates)) %>% 
    group_by(Mass)
  
  
  make_graph(Data, metabolite_name, DataDirectory)
}


```


