---
title: "Metabolomics analysis intracellular metabolites"
author: "Pauline Mencke"
date: "19/8/2022"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


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
library(dlookr)
library(png)
library(grid)
library(scales)
library(magick)
library(ggpattern)
theme_set(theme_classic())
knitr::opts_chunk$set(echo = TRUE)

```

```{r Loading of Data}
#load data table containing MIDs
setwd('C:/Users/pauline.mencke/Documents/Pauline/Lab/1. Main paper thesis/Astrocytes/metabolite tracing/LCMS/Gln LCMS 10052022/')
#remove nans in Excel with ctrl H 
TableImport = readxl::read_excel('C:/Users/pauline.mencke/Documents/Pauline/Lab/1. Main paper thesis/Astrocytes/metabolite tracing/LCMS/Gln LCMS 10052022/Raw data/LC_MS_Astrocytes_Gln_tracing.xlsx')



  

```

```{r Preparing tables}
#melt the table so that the data will be tidy (every experimental condition now has one value in a column)
Table_all <- TableImport %>%  
  dplyr::select(metabolite, sample, isotopologue, isotopologue_fraction) %>% 
  separate(sample, into = c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8"), sep = "_") %>% 
  mutate(Full_Condition_Name = glue::glue("{V6}")) %>% 
  #glue the name with the info needed   
  #create new column indicating the technical replicate and general name (without number at the end)
  mutate(Technical_Replicate = V8) %>%
  mutate(Biological_Replicate = V7) %>% 
  mutate(Condition = glue::glue("{V5}")) %>%
  rename(Name = Full_Condition_Name) %>%  
  rename(Mass = isotopologue) %>%  
  rename(value = isotopologue_fraction) %>%
  #make an easy overview table
  dplyr::select(metabolite, Mass, Name, Condition, Technical_Replicate, Biological_Replicate, value)




Table_areas = TableImport %>%  
  dplyr::select(metabolite, sample, isotopologue, corrected_area) %>% 
  separate(sample, into = c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8"), sep = "_") %>% 
  mutate(Full_Condition_Name = glue::glue("{V4}_{V6}_{V5}")) %>% 
  #glue the name with the info needed   
  #create new column indicating the technical replicate and general name (without number at the end)
  mutate(Technical_Replicate = V8) %>%
  mutate(Biological_Replicate = V7) %>%  
  rename(Name = Full_Condition_Name) %>%  
  rename(Mass = isotopologue) %>%  
  #make an easy overview table
  dplyr::select(metabolite, Mass, Name, Technical_Replicate, Biological_Replicate, corrected_area)
                

```

```{r Statistics - Calculating means}

Table_Mean_Technical_Replicates <- Table_all %>%
  group_by(metabolite, Mass, Name, Biological_Replicate) %>%
  mutate("Mean_Technical_Replicates" = mean(value),
         "Stdev_Technical_Replicates" = sd(value)) %>%
    #keep only one value for each mean of the technical replicates and shrink the table by using the distinct function
  distinct(metabolite, Mass, Name, .keep_all = TRUE)

Table_Mean_Biological_Replicates <- Table_Mean_Technical_Replicates %>%
  group_by(metabolite, Mass, Name) %>%
  mutate("Mean_Biological_Replicates" = mean(Mean_Technical_Replicates),
         "Stdev_Biological_Replicates" = sd(Mean_Technical_Replicates)) %>%  
  #keep only one value for each mean of the technical replicates and shrink the table by using the distinct function
  distinct(metabolite, Mass, Name, .keep_all = TRUE)



Table_Mean_Technical_Replicates_area <- Table_areas %>%
  group_by(metabolite, Mass, Name, Biological_Replicate) %>%
  mutate("Mean_Technical_Replicates" = mean(corrected_area),
         "Stdev_Technical_Replicates" = sd(corrected_area)) %>%
    #keep only one value for each mean of the technical replicates and shrink the table by using the distinct function
  distinct(metabolite, Mass, Name, .keep_all = TRUE)

Table_Mean_Biological_Replicates_area <- Table_Mean_Technical_Replicates_area %>%
  group_by(metabolite, Mass, Name) %>%
  mutate("Mean_Biological_Replicates" = mean(Mean_Technical_Replicates),
         "Stdev_Biological_Replicates" = sd(Mean_Technical_Replicates)) %>%  
  #keep only one value for each mean of the technical replicates and shrink the table by using the distinct function
  distinct(metabolite, Mass, Name, .keep_all = TRUE)

write_csv(Table_Mean_Biological_Replicates_area, "C:/Users/pauline.mencke/Documents/Pauline/Lab/1. Main paper thesis/Astrocytes/metabolite tracing/LCMS/Gln LCMS 10052022/Results/Table_Mean_Biological_Replicates_area.csv")


Table_Sums_areas <- Table_Mean_Biological_Replicates_area %>%
  group_by(metabolite, Name) %>%
  mutate("Sum_area" = sum(Mean_Biological_Replicates)) %>%  
  #keep only one value for each mean of the technical replicates and shrink the table by using the distinct function
  distinct(metabolite, Name, .keep_all = TRUE)


write_csv(Table_Sums_areas, "C:/Users/pauline.mencke/Documents/Pauline/Lab/1. Main paper thesis/Astrocytes/metabolite tracing/LCMS/Gln LCMS 10052022/Results/Table_Sums_areas.csv")




```

```{r Calculating Ratios}



#M+4 citric acid / M+2 citric acid
M4citricacid <- Table_Mean_Technical_Replicates %>%  
  filter(metabolite == "Citric acid") %>%  
  filter(Mass == "4")
M2citricacid <- Table_Mean_Technical_Replicates %>%  
  filter(metabolite == "Citric acid") %>%  
  filter(Mass == "2")
M4citricacidM2citricacid <- cbind(M4citricacid, M2citricacid) %>%  
    mutate(M4citricacidM2citricacidRatio = Mean_Technical_Replicates...8/Mean_Technical_Replicates...17) %>%  
  select(Name...3, Technical_Replicate...5, Biological_Replicate...6, M4citricacidM2citricacidRatio) %>%  
  rename(Name = Name...3) %>% 
  rename(Technical_Replicate = Technical_Replicate...5) %>% 
  rename(Biological_Replicate = Biological_Replicate...6) 

  
write_csv(M4citricacidM2citricacid, "C:/Users/pauline.mencke/Documents/Pauline/Lab/1. Main paper thesis/Astrocytes/metabolite tracing/LCMS/Gln LCMS 10052022/Results/M4citricacidM2citricacidRatio.csv")



#M+4 succinic acid / M+5 glutamic acid
M4succinicacid <- Table_Mean_Technical_Replicates %>%  
  filter(metabolite == "Succinic acid") %>%  
  filter(Mass == "4")
M5glutamicacid <- Table_Mean_Technical_Replicates %>%  
  filter(metabolite == "L-Glutamic acid") %>%  
  filter(Mass == "5")
M4succinicacidM5glutamicacid <- cbind(M4succinicacid, M5glutamicacid) %>%  
    mutate(M4succinicacidM5glutamicacidRatio = Mean_Technical_Replicates...8/Mean_Technical_Replicates...17) %>%  
  select(Name...3, Technical_Replicate...5, Biological_Replicate...6, M4succinicacidM5glutamicacidRatio) %>%  
  rename(Name = Name...3) %>% 
  rename(Technical_Replicate = Technical_Replicate...5) %>% 
  rename(Biological_Replicate = Biological_Replicate...6) 

  
write_csv(M4succinicacidM5glutamicacid, "C:/Users/pauline.mencke/Documents/Pauline/Lab/1. Main paper thesis/Astrocytes/metabolite tracing/LCMS/Gln LCMS 10052022/Results/M4succinicacidM5glutamicacidRatio.csv")





#M+5 glutamic acid / M+4 succinic acid
M5glutamicacid <- Table_Mean_Technical_Replicates %>%  
  filter(metabolite == "L-Glutamic acid") %>%  
  filter(Mass == "5")
M4succinicacid <- Table_Mean_Technical_Replicates %>%  
  filter(metabolite == "Succinic acid") %>%  
  filter(Mass == "4")
M5glutamicacidM4succinicacid <- cbind(M5glutamicacid, M4succinicacid) %>%  
    mutate(M5glutamicacidM4succinicacidRatio = Mean_Technical_Replicates...8/Mean_Technical_Replicates...17) %>%  
  select(Name...3, Technical_Replicate...5, Biological_Replicate...6, M5glutamicacidM4succinicacidRatio) %>%  
  rename(Name = Name...3) %>% 
  rename(Technical_Replicate = Technical_Replicate...5) %>% 
  rename(Biological_Replicate = Biological_Replicate...6) 

  
write_csv(M5glutamicacidM4succinicacid, "C:/Users/pauline.mencke/Documents/Pauline/Lab/1. Main paper thesis/Astrocytes/metabolite tracing/LCMS/Gln LCMS 10052022/Results/M5glutamicacidM4succinicacidRatio.csv")



```

```{r Exporting of data}

write_csv(Table_Mean_Technical_Replicates, "C:/Users/pauline.mencke/Documents/Pauline/Lab/1. Main paper thesis/Astrocytes/metabolite tracing/LCMS/Gln LCMS 10052022/ResultsLC_MS_Astrocytes_Gln_tracing_results_mean_technical_reps.csv")
write_csv(Table_Mean_Biological_Replicates, "C:/Users/pauline.mencke/Documents/Pauline/Lab/1. Main paper thesis/Astrocytes/metabolite tracing/LCMS/Gln LCMS 10052022/Results/LC_MS_Astrocytes_Gln_tracing_results_mean_biological_reps.csv")

Table_Mean_Technical_Replicates <- Table_Mean_Technical_Replicates %>%  
  filter (Mass == 0)

write_csv(Table_Mean_Technical_Replicates, "C:/Users/pauline.mencke/Documents/Pauline/Lab/1. Main paper thesis/Astrocytes/metabolite tracing/LCMS/Gln LCMS 10052022/Results/LC_MS_Astrocytes_Gln_tracing_results_mean_technical_reps_M0.csv")

FinalTable <- Table_Mean_Biological_Replicates #%>%  
  #filter (Mass == 0)

FinalTable$Name <- factor(FinalTable$Name, levels = c("C4", "C4mut", "DelPGC", "DelP"))

 
FinalTable <- na.omit(FinalTable)

```


```{r Automated image generator - Plotting of data} {r automated plots including stats, fig.height= 8, fig.width= 15}


#define where you want to have your plots
DataDirectory <- "C:/Users/pauline.mencke/Documents/Pauline/Lab/1. Main paper thesis/Astrocytes/metabolite tracing/LCMS/Gln LCMS 10052022/Results/MID_plots/M_all_new_colors/"
#define the graph function
 make_graph <- function(Data, x, DataDirectory){#, stat.test){
   #define the metabolite name and set a name for the file
   metabolite_name <- x
   name <- paste(DataDirectory, metabolite_name, ".png", sep="")
   #open image creating part of the function, set resolution
   png(name, width = 720, height = 720, units = "px")
   #make the graph
   
   Graph <-
 ggplot(Data, aes(x = Mass, y = Mean_Biological_Replicates, col = Name))+
  geom_bar(aes(fill = Name), stat="identity", colour = "black", position=position_dodge(0.8))+
  geom_errorbar(aes(ymin = Ymin, ymax = Ymax), width=.2,
                 position=position_dodge(.8))+
  #geom_errorbar(aes(ymin = Ymin, ymax = Ymax), position=position_dodge(width = 0.9), colour="black")+
  #geom_errorbar(aes(ymin = Ymin, ymax = Ymax, fill = Name), stat = "identity", width=0.8, position=position_dodge2(0.8), colour = 'black')+
  #geom_point(position=position_dodge(width=0.8))+
  ggtitle(metabolite_name)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=40))+
  theme(axis.text.x = element_text(size=20))+
  theme(axis.text.y = element_text(size=20))+
  theme(legend.text=element_text(size=20))+
  theme(axis.title.x = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20))+
  theme(legend.title = element_blank())+
  theme(legend.text = element_text(size=10))+
  # theme(legend.position = "none")+
     
    scale_x_continuous(breaks=0:length(FinalTable$Mass))+


  #set colours
  scale_fill_manual(values=c("chartreuse3", "blue2", "darkolivegreen2", "darkslategray2")) #+
  # # #scale_fill_manual(values=c("orange", "blue", "cyan", "red", "green", "yellow"))+
  # #                            #, "red", "green", "yellow", "red", "green", "yellow", "red", "green", "yellow")) +


  # #add statistics
  # #stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.001, bracket.nudge.y = 0, dodge = 0.8, hide.ns = TRUE)
  # #by printing, the graph is saved at the data directory
  print(Graph)
  #the png module is closed
  dev.off()
}

#run the for loop
for(metabolite_name in unique(FinalTable$metabolite)){

    Data <- FinalTable %>%
    filter(metabolite == metabolite_name) %>%
    mutate(Ymin = (Mean_Biological_Replicates-Stdev_Biological_Replicates),
    Ymax = (Mean_Biological_Replicates+Stdev_Biological_Replicates)) #%>%
    #
    #group_by(Mass)

    # pwc[[metabolite_name]] %>%
    # add_xy_position(x = "Name", fun = "max",  step.increase = 0.05, dodge = 0.8) -> stat.test

  make_graph(Data, metabolite_name, DataDirectory)#, stat.test = stat.test)
}



# setwd("C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/IGC metabolomics/Astrocytes and GBM cells/Results/Analysis_IGC_GBM_Astros/Barplots_Astrocytes/n2_3_4/")
# pdf('Summary_all_metabolites_astro_all_ns.pdf')
# lapply(ll <- list.files(patt='.*[.]png'),function(x){
#   img <- as.raster(readPNG(x))
#   grid.newpage()
#   grid.raster(img, interpolate = FALSE)
#
# })
# dev.off()



``` 





