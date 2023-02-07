---
title: "Metabolomics analysis intracellular metabolites"
author: "Pauline Mencke"
date: "1/8/2021"
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
setwd('C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/R analysis metabolomics/Intracellular/')
#remove nans in Excel with ctrl H 
TableImportN1 = readxl::read_excel('C:/Users/pauline.mencke/Documents/Pauline/Lab/1. Main paper thesis/Astrocytes/metabolite tracing/GCMS/IGC Astrocytes paper/GlcC13_N1/raw data/20210211_NHZDR_PMe13_IGC_MIDs.xlsx', range = "A1:AF386")
#remove all NAs (this removes also the disturbing metabolite names in the mass column)
TableN1 = na.omit(TableImportN1)
TableImportN2 = readxl::read_excel('C:/Users/pauline.mencke/Documents/Pauline/Lab/1. Main paper thesis/Astrocytes/metabolite tracing/GCMS/IGC Astrocytes paper/GlcC13_N2/raw data/2021024_0BMYB_PMe14_IGC_MIDs.xlsx', range = "A1:AL400")
#remove all NAs (this removes also the disturbing metabolite names in the mass column)
TableN2 = na.omit(TableImportN2)
TableImportN3 = readxl::read_excel('C:/Users/pauline.mencke/Documents/Pauline/Lab/1. Main paper thesis/Astrocytes/metabolite tracing/GCMS/IGC Astrocytes paper/GlcC13_N3/raw data/20210304_7AE94_PMe15_IGC_MIDs.xlsx', range = "A1:AF368")
#remove all NAs (this removes also the disturbing metabolite names in the mass column)
TableN3 = na.omit(TableImportN3)
TableImportN4_5 = readxl::read_excel('C:/Users/pauline.mencke/Documents/Pauline/Lab/1. Main paper thesis/Astrocytes/metabolite tracing/GCMS/IGC Astrocytes paper/GlcC13_N4_N5/raw data/20210503_TBSLB_24KSC_PMe16_PMe17_IGC_MID_formatted.xlsx', range = "A1:BP392")
#remove all NAs (this removes also the disturbing metabolite names in the mass column)
TableN4_5 = na.omit(TableImportN4_5)

```

```{r Preparing tables}
#melt the table so that the data will be tidy (every experimental condition now has one value in a column)
Table_V0N1 <- melt(TableN1, id.vars=c("Metabolite_name", "Mass")) %>%  
  rename(Full_Condition_Name = variable) %>% 
#separate metabolite name into metabolite name only 
  separate(Metabolite_name, into = c("metabolite", "other"), sep = "_") 
#silence other  
Table_V0N1$other = NULL
#separate mass name into mass name only 
Table_V1N1 <- Table_V0N1 %>% 
  separate(Mass, into = c("Number", "Mass"), sep = " ")
#silence Number 
Table_V1N1$Number = NULL
#separate full condiiton name so that the info on technical replicates and biological replicates can be extracted 
Table_V2_N1 <- Table_V1N1 %>% 
  separate(Full_Condition_Name, into = c("V1", "V2", "V3", "V4", "V5", "V6"), sep = "_") %>% 
  mutate(Full_Condition_Name = glue::glue("{V3}_{V4}_{V5}")) %>% 
  #glue the name with the info needed   
  #create new column indicating the technical replicate and general name (without number at the end)
  mutate(Technical_Replicate = V6) %>%
  rename(Name = Full_Condition_Name) %>%
  mutate(condition = glue::glue("{V5}")) %>% 
  mutate(Biological_Replicate = 1) %>%  
#make an easy overview table
  dplyr::select(metabolite, Mass, Name, condition, Technical_Replicate, Biological_Replicate, value)



#melt the table so that the data will be tidy (every experimental condition now has one value in a column)
Table_V0N2 <- melt(TableN2, id.vars=c("Metabolite_name", "Mass")) %>%  
  rename(Full_Condition_Name = variable) %>% 
#separate metabolite name into metabolite name only 
  separate(Metabolite_name, into = c("metabolite", "other"), sep = "_") 
#silence other  
Table_V0N2$other = NULL
#separate mass name into mass name only 
Table_V1N2 <- Table_V0N2 %>% 
  separate(Mass, into = c("Number", "Mass"), sep = " ")
#silence Number 
Table_V1N2$Number = NULL
#separate full condiiton name so that the info on technical replicates and biological replicates can be extracted 
Table_V2_N2 <- Table_V1N2 %>% 
  separate(Full_Condition_Name, into = c("V1", "V2", "V3", "V4", "V5", "V6"), sep = "_") %>% 
  mutate(Full_Condition_Name = glue::glue("{V3}_{V4}_{V5}")) %>% 
  #glue the name with the info needed   
  #create new column indicating the technical replicate and general name (without number at the end)
  mutate(Technical_Replicate = V6) %>%
  rename(Name = Full_Condition_Name) %>%
  mutate(condition = glue::glue("{V5}")) %>% 
  mutate(Biological_Replicate = 2) %>% 
#make an easy overview table
  dplyr::select(metabolite, Mass, Name, condition, Technical_Replicate, Biological_Replicate, value)




#melt the table so that the data will be tidy (every experimental condition now has one value in a column)
Table_V0N3 <- melt(TableN3, id.vars=c("Metabolite_name", "Mass")) %>%  
  rename(Full_Condition_Name = variable) %>% 
#separate metabolite name into metabolite name only 
  separate(Metabolite_name, into = c("metabolite", "other"), sep = "_") 
#silence other  
Table_V0N3$other = NULL
#separate mass name into mass name only 
Table_V1N3 <- Table_V0N3 %>% 
  separate(Mass, into = c("Number", "Mass"), sep = " ")
#silence Number 
Table_V1N3$Number = NULL
#separate full condiiton name so that the info on technical replicates and biological replicates can be extracted 
Table_V2_N3 <- Table_V1N3 %>% 
  separate(Full_Condition_Name, into = c("V1", "V2", "V3", "V4", "V5", "V6"), sep = "_") %>% 
  mutate(Full_Condition_Name = glue::glue("{V3}_{V4}_{V5}")) %>% 
  #glue the name with the info needed   
  #create new column indicating the technical replicate and general name (without number at the end)
  mutate(Technical_Replicate = V6) %>%
  rename(Name = Full_Condition_Name) %>%
  mutate(condition = glue::glue("{V5}")) %>% 
  mutate(Biological_Replicate = 3) %>% 
#make an easy overview table
  dplyr::select(metabolite, Mass, Name, condition, Technical_Replicate, Biological_Replicate, value)




#melt the table so that the data will be tidy (every experimental condition now has one value in a column)
Table_V0N4_5 <- melt(TableN4_5, id.vars=c("Metabolite_name", "Mass")) %>%  
  rename(Full_Condition_Name = variable) %>% 
#separate metabolite name into metabolite name only 
  separate(Metabolite_name, into = c("metabolite", "other"), sep = "_") 
Table_V0N4_5 <- na.omit(Table_V0N4_5)
#silence other  
Table_V0N4_5$other = NULL
#separate mass name into mass name only 
Table_V1N4_5 <- Table_V0N4_5 %>% 
  separate(Mass, into = c("Number", "Mass"), sep = " ")
#silence Number 
Table_V1N4_5$Number = NULL
#separate full condiiton name so that the info on technical replicates and biological replicates can be extracted 
Table_V2_N4_5 <- Table_V1N4_5 %>% 
  separate(Full_Condition_Name, into = c("V1", "V2", "V3", "V4", "V5", "V6"), sep = "_") %>% 
  mutate(Full_Condition_Name = glue::glue("{V3}_{V4}_{V5}")) %>% 
  #glue the name with the info needed   
  #create new column indicating the technical replicate and general name (without number at the end)
  mutate(Technical_Replicate = V6) %>%
  rename(Name = Full_Condition_Name) %>%
  mutate(condition = glue::glue("{V5}")) %>% 
  mutate(Biological_Replicate = case_when(V1 == "PMe16"~ 4,
                                          V1 == "PMe17"~ 5)) %>% 
#make an easy overview table
  dplyr::select(metabolite, Mass, Name, condition, Technical_Replicate, Biological_Replicate, value)





Table_all <- rbind(Table_V2_N1, Table_V2_N2, Table_V2_N3, Table_V2_N4_5) #%>%  
  #filter(Mass == "(M0)") %>%  
  #filter(condition == "normoxia")  
  



```

```{r Statistics - Removing outliers}

#calculate the mean of the 9 values for each metabolite (technical and biological triplicates)
Table_Means <- Table_all %>%
  group_by(metabolite, Mass, Name) %>%
  mutate("Mean" = mean(value),
         "Stdev" = sd(value)) %>%
#estimate the coefficient of variation (CV=standard deviation / mean). As a rule of thumb, a CV >= 1 indicates a relatively high variation.
  mutate(CV=Stdev/Mean) %>%
  #mutate(deviation = ifelse(CV >= 1, TRUE, FALSE))
  mutate(deviation = ifelse(Stdev == (0.05*Mean), TRUE, FALSE))
#create a table with no deviated data
Table_no_outliers <- Table_Means %>%
  group_by(metabolite, Mass, Name) %>%
  filter(deviation == FALSE)
#create a table with deviated data
Table_outliers <- Table_Means%>%
  #keep only the rows with CV >= 1
  filter(deviation == TRUE)
#remove outliers in deviated data
#find Q1, Q3, and interquartile range for values in column A
Q1 <- quantile(Table_outliers$value, .25)
Q3 <- quantile(Table_outliers$value, .75)
IQR <- IQR(Table_outliers$value)
#only keep rows in dataframe that have values within 1.5*IQR of Q1 and Q3
Table_outliers_removed <- subset(Table_outliers, Table_outliers$value > (Q1 - 1.5*IQR) & Table_outliers$value < (Q3 + 1.5*IQR))
#bind data to obtain final df that contains only non deviated data
Data_for_statistics <- rbind(Table_no_outliers) %>%  #, Table_outliers_removed) %>%
  select(metabolite, Mass, Name, condition, Technical_Replicate, Biological_Replicate, value)



```


```{r Ratios}


# M+2 citrate / M+3 pyruvate
# M+3 aspartate / M+3 pyruvate
# M+4 citrate / M+2 citrate
# M2 Proline / M2 Glutamate 
# M3 Pyruvate / M3 3PG 
# M3 3PG / M2 Citrate 
# M3 Serine / M3 3PG
# M3 3PG / M6 Glc geht nicht weil Glc nicht IGC gemessen wird 



Table_for_Ratios <- Data_for_statistics %>%
  filter(Name %in% c("Astrocytes_wt_normoxia", "Astrocytes_A8mut_normoxia","Astrocytes_DelPGC13_normoxia", "Astrocytes_DelP_normoxia", "Astrocytes_wt_hypoxia", "Astrocytes_A8mut_hypoxia","Astrocytes_DelPGC13_hypoxia", "Astrocytes_DelP_hypoxia")) %>%  
  filter(condition == "normoxia") 


#M+2 aKG / M+2 Glu
M2aKG <- Table_for_Ratios %>%  
  filter(metabolite == "a-Ketoglutarate") %>%  
  filter(Mass == "(M2)")
M2Glu <- Table_for_Ratios %>%  
  filter(metabolite == "Glutamate") %>%  
  filter(Mass == "(M2)")
M2aKGM2Glu <- cbind(M2aKG, M2Glu) %>%  
  filter( value...7 > 0) %>% 
  filter( value...14 > 0) %>% 
  mutate(M2aKGM2GluRatio = value...7/value...14) %>%  
  select(Name...3, condition...4, Technical_Replicate...5, Biological_Replicate...6, M2aKGM2GluRatio) %>%  
  rename(Name = Name...3) %>% 
  rename(condition = condition...4) %>%
  rename(Technical_Replicate = Technical_Replicate...5) %>% 
  rename(Biological_Replicate = Biological_Replicate...6) 
# %>%  
#   group_by(Name) %>%  
#   mutate(MeanM2CitrateM3PyruvateRatio = mean(M2CitrateM3PyruvateRatio))
  
write_csv(M2aKGM2Glu, "C:/Users/pauline.mencke/Documents/Pauline/Lab/1. Main paper thesis/Astrocytes/metabolite tracing/GCMS/IGC Astrocytes paper/Ratios/M2aKGM2Glu.csv")




#M+2 citrate / M+3 pyruvate
M2Citrate <- Table_for_Ratios %>%  
  filter(metabolite == "Citrate") %>%  
  filter(Mass == "(M2)")
M3Pyruvate <- Table_for_Ratios %>%  
  filter(metabolite == "Pyruvate") %>%  
  filter(Mass == "(M3)")
M2CitrateM3Pyruvate <- cbind(M2Citrate, M3Pyruvate) %>%  
  filter( value...7 > 0) %>% 
  filter( value...14 > 0) %>% 
  mutate(M2CitrateM3PyruvateRatio = value...7/value...14) %>%  
  select(Name...3, condition...4, Technical_Replicate...5, Biological_Replicate...6, M2CitrateM3PyruvateRatio) %>%  
  rename(Name = Name...3) %>% 
  rename(condition = condition...4) %>%
  rename(Technical_Replicate = Technical_Replicate...5) %>% 
  rename(Biological_Replicate = Biological_Replicate...6) 
# %>%  
#   group_by(Name) %>%  
#   mutate(MeanM2CitrateM3PyruvateRatio = mean(M2CitrateM3PyruvateRatio))
  
write_csv(M2CitrateM3Pyruvate, "C:/Users/pauline.mencke/Documents/Pauline/Lab/1. Main paper thesis/Astrocytes/metabolite tracing/GCMS/IGC Astrocytes paper/Ratios/M2CitrateM3Pyruvate.csv")



#M+3 aspartate / M+3 pyruvate
M3Aspartate <- Table_for_Ratios %>%  
  filter(metabolite == "Aspartate") %>%  
  filter(Mass == "(M3)")
M3Pyruvate <- Table_for_Ratios %>%  
  filter(metabolite == "Pyruvate") %>%  
  filter(Mass == "(M3)")
M3AspartateM3Pyruvate <- cbind(M3Aspartate, M3Pyruvate) %>% 
  filter( value...7 > 0) %>% 
  filter( value...14 > 0) %>% 
  mutate(M3AspartateM3PyruvateRatio = value...7/value...14) %>%  
  select(Name...3, condition...4, Technical_Replicate...5, Biological_Replicate...6, M3AspartateM3PyruvateRatio) %>%  
  rename(Name = Name...3) %>% 
  rename(condition = condition...4) %>%
  rename(Technical_Replicate = Technical_Replicate...5) %>% 
  rename(Biological_Replicate = Biological_Replicate...6) 
  
write_csv(M3AspartateM3Pyruvate, "C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/IGC metabolomics/Astrocytes paper/Ratios/M3AspartateM3Pyruvate.csv")




#M+4 citrate / M+2 citrate
M4Citrate <- Table_for_Ratios %>%  
  filter(metabolite == "Citrate") %>%  
  filter(Mass == "(M4)")
M2Citrate <- Table_for_Ratios %>%  
  filter(metabolite == "Citrate") %>%  
  filter(Mass == "(M2)")
M4CitrateM2Citrate <- cbind(M4Citrate, M2Citrate) %>%
  filter( value...7 > 0) %>% 
  filter( value...14 > 0) %>% 
  mutate(M4CitrateM2CitrateRatio = value...7/value...14) %>%  
  select(Name...3, condition...4, Technical_Replicate...5, Biological_Replicate...6, M4CitrateM2CitrateRatio) %>%  
  rename(Name = Name...3) %>% 
  rename(condition = condition...4) %>%
  rename(Technical_Replicate = Technical_Replicate...5) %>% 
  rename(Biological_Replicate = Biological_Replicate...6) 
  
write_csv(M4CitrateM2Citrate, "C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/IGC metabolomics/Astrocytes paper/Ratios/M4CitrateM2Citrate.csv")






#M2 Proline/ M2 Glutamate
M2Proline <- Table_for_Ratios %>%  
  filter(metabolite == "Proline") %>%  
  filter(Mass == "(M2)")
M2Glutamate <- Table_for_Ratios %>%  
  filter(metabolite == "Glutamate") %>%  
  filter(Mass == "(M2)")
M2ProlineM2Glutamate <- cbind(M2Proline, M2Glutamate) %>%  
  filter( value...7 > 0) %>% 
  filter( value...14 > 0) %>% 
  mutate(M2ProlineM2GlutamateRatio = value...7/value...14) %>%  
  select(Name...3, condition...4, Technical_Replicate...5, Biological_Replicate...6, M2ProlineM2GlutamateRatio) %>%  
  rename(Name = Name...3) %>% 
  rename(condition = condition...4) %>%
  rename(Technical_Replicate = Technical_Replicate...5) %>% 
  rename(Biological_Replicate = Biological_Replicate...6) 
  
write_csv(M2ProlineM2Glutamate, "C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/IGC metabolomics/Astrocytes paper/Ratios/M2ProlineM2Glutamate.csv")





#M3 Pyruvate / M3 3PG
M3Pyruvate <- Table_for_Ratios %>%  
  filter(metabolite == "Pyruvate") %>%  
  filter(Mass == "(M3)")
M33PG <- Table_for_Ratios %>%  
  filter(metabolite == "3PG") %>%  
  filter(Mass == "(M3)")
M3PyruvateM33PG <- cbind(M3Pyruvate, M33PG) %>%  
  filter( value...7 > 0) %>% 
  filter( value...14 > 0) %>% 
  mutate(M3PyruvateM33PGRatio = value...7/value...14) %>%  
  select(Name...3, condition...4, Technical_Replicate...5, Biological_Replicate...6, M3PyruvateM33PGRatio) %>%  
  rename(Name = Name...3) %>% 
  rename(condition = condition...4) %>%
  rename(Technical_Replicate = Technical_Replicate...5) %>% 
  rename(Biological_Replicate = Biological_Replicate...6) 
  
write_csv(M3PyruvateM33PG, "C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/IGC metabolomics/Astrocytes paper/Ratios/M3PyruvateM33PG.csv")




#M3 3PG/ M2 Citrate 
M33PG <- Table_for_Ratios %>%  
  filter(metabolite == "3PG") %>%  
  filter(Mass == "(M3)")
M2Citrate <- Table_for_Ratios %>%  
  filter(metabolite == "Citrate") %>%  
  filter(Mass == "(M2)")
M33PGM2Citrate <- cbind(M33PG, M2Citrate) %>% 
  filter( value...7 > 0) %>% 
  filter( value...14 > 0) %>% 
  mutate(M33PGM2CitrateRatio = value...7/value...14) %>%  
  select(Name...3, condition...4, Technical_Replicate...5, Biological_Replicate...6, M33PGM2CitrateRatio) %>%  
  rename(Name = Name...3) %>% 
  rename(condition = condition...4) %>%
  rename(Technical_Replicate = Technical_Replicate...5) %>% 
  rename(Biological_Replicate = Biological_Replicate...6) 
  
write_csv(M33PGM2Citrate, "C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/IGC metabolomics/Astrocytes paper/Ratios/M33PGM2Citrate.csv")





#M3 Serine/ M3 3PG
M3Serine <- Table_for_Ratios %>%  
  filter(metabolite == "Serine") %>%  
  filter(Mass == "(M3)")
M33PG <- Table_for_Ratios %>%  
  filter(metabolite == "3PG") %>%  
  filter(Mass == "(M3)")
M3SerineM33PG <- cbind(M3Serine, M33PG) %>%  
  filter( value...7 > 0) %>% 
  filter( value...14 > 0) %>% 
  mutate(M3SerineM33PGRatio = value...7/value...14) %>%  
  select(Name...3, condition...4, Technical_Replicate...5, Biological_Replicate...6, M3SerineM33PGRatio) %>%  
  rename(Name = Name...3) %>% 
  rename(condition = condition...4) %>%
  rename(Technical_Replicate = Technical_Replicate...5) %>% 
  rename(Biological_Replicate = Biological_Replicate...6) 
  
write_csv(M3SerineM33PG, "C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/IGC metabolomics/Astrocytes paper/Ratios/M3SerineM33PG.csv")






```

```{r Statistics - Normality test and ANOVA multiple comparisons}
#calculate the mean of the technical replicates
Table_Mean_Technical_Replicates <- Data_for_statistics %>%
  group_by(metabolite, Mass, Name, Biological_Replicate) %>%
  mutate("Mean_Technical_Replicates" = mean(value),
         "Stdev_Technical_Replicates" = sd(value)) %>%  
  mutate(deviation = ifelse(Stdev_Technical_Replicates == (0.05*Mean_Technical_Replicates), TRUE, FALSE)) %>%   
  filter(deviation == FALSE) %>%  
   
  mutate(deviation2 = ifelse(value == (0.5+Mean_Technical_Replicates), TRUE, FALSE)) %>%  
  filter(deviation2 == FALSE) %>%  
  mutate(deviation3 = ifelse(value == (0.5-Mean_Technical_Replicates), TRUE, FALSE)) %>%  
  filter(deviation3 == FALSE) #%>% 
  #filter(Mass == "(M0)") #%>%  
  #filter(condition == "normoxia")
  #keep only one value for each mean of the technical replicates and shrink the table by using the distinct function
  #distinct(metabolite, Mass, Name, condition, Biological_Replicate,.keep_all = TRUE)

write_csv(Table_Mean_Technical_Replicates, "C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/IGC metabolomics/Astrocytes paper/Astrocytes_intracellular_results_mean_technical_rep_n1_5_all_values_normoxia_hypoxia.csv")


#Table_Mean_Technical_Replicates <- Table_Mean_Technical_Replicates %>%
  #filter(Mass == "(M0)") %>%
  #filter(condition == "normoxia") %>%
  #filter(Biological_Replicate == "1") %>%  
  #filter(Biological_Replicate %in% c("4", "5")) %>%  
  #filter(Name %in% c("Astrocytes_wt_normoxia", "Astrocytes_A8mut_normoxia","Astrocytes_DelPGC13_normoxia", "Astrocytes_DelP_normoxia"))
  # filter(Name %in% c("Astrocytes_wt_normoxia", "Astrocytes_A8mut_normoxia","Astrocytes_DelPGC13_normoxia", "Astrocytes_DelP_normoxia", "Astrocytes_wt_hypoxia", "Astrocytes_A8mut_hypoxia","Astrocytes_DelPGC13_hypoxia", "Astrocytes_DelP_hypoxia"))


#write_csv(Table_Mean_Technical_Replicates, "C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/IGC metabolomics/Astrocytes paper/Astrocytes_intracellular_results_mean_technical_rep_n1_5_all_values.csv")


Table_Mean_Technical_Replicates <- Table_Mean_Technical_Replicates %>%
  distinct(metabolite, Mass, Name, condition, Biological_Replicate,.keep_all = TRUE) %>%  
mutate(filterName = glue::glue("{metabolite}_{Mass}_{Name}")) #%>%  
  #filter(Mass == "(M0)") %>%
  #filter(condition == "normoxia")
#write_csv(Table_Mean_Technical_Replicates, "C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/IGC metabolomics/Astrocytes paper/Astrocytes_intracellular_results_mean_technical_rep_n1_5_new_outlier_test.csv")
#(Table_Mean_Technical_Replicates, "C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/IGC metabolomics/Astrocytes paper/Astrocytes_intracellular_results_mean_technical_rep_n1_5.csv")
  




# 
# 
# 
# # Transform the data into long format
# Table_with_meanslong <- Table_Mean_Technical_Replicates %>%  
#   select(Name, Mean_Technical_Replicates) %>%  
#   mutate(Name = glue::glue("{metabolite}_{Name}"))
# 
# Table_with_meanslong <- Table_with_meanslong[4:5]
# Table_with_meanslong <- as.data.frame(Table_with_meanslong) %>%
# pivot_longer(-Name, names_to = "variables", values_to = "value")
# Table_with_meanslong %>% sample_n(1075)
# 
# Table_with_meanslong <- Table_with_meanslong %>%  
#   separate(Name, into = c("V1", "V2", "V3", "V4"), sep = "_") %>%
#   mutate(metabolite = glue::glue("{V1}")) %>% 
#   mutate(CellLine = glue::glue("{V3}")) %>% 
#   mutate(Condition = glue::glue("{V4}")) %>% 
#   select(metabolite, CellLine, value) %>%  
#   ungroup() %>%  
#   group_by(metabolite, CellLine) %>%  
#   mutate(Mean = mean(value)) %>%  
#   distinct(metabolite, CellLine, Mean)  
# 
# stat.test <- Table_with_meanslong %>%
#   group_by(metabolite) %>%
#   t_test(Mean ~ CellLine) %>%
#   adjust_pvalue(method = "BH") %>%
#   add_significance()
# #stat.test
# stat.test <- Table_with_meanslong %>%
#  group_by(metabolite) %>%
#  emmeans_test(value ~ CellLine)
# 
# 
# #perform normality test prior to anova 
# stat.test.normality <- Table_with_meanslong %>%
#   group_by(CellLine, Condition, variables) %>%
# normality(value) %>%
# #keep only the values where p_value is greater than 0.01
#   filter(p_value > 0.01)
















# 
# #perform normality test of data before pushing the data into ANOVA or dunn test
# Data_normality_test <- Table_Mean_Technical_Replicates %>%
#   group_by(metabolite, Mass, Name) %>%
#   normality(Mean_Technical_Replicates) %>%
# #   #keep only the values where p_value is greater than 0.01
#   filter(p_value > 0.01) %>%  
# mutate(filterName = glue::glue("{metabolite}_{Mass}_{Name}"))
# #obtain the values for the normally distributed data
# Normally_distributed_data <- Table_Mean_Technical_Replicates %>%
#   filter(filterName %in% Data_normality_test$filterName) %>%
#   ungroup() %>%
#   filter(Mass == "(M0)") %>%
#   filter(condition == "normoxia")

#order the sample names in the desired order before doing the statistics
#as otherwise the statistics will not fit the order of the graph
# Normally_distributed_data$Name <- factor(Normally_distributed_data$Name, levels = c("Astrocytes_wt_normoxia", "Astrocytes_A8mut_normoxia","Astrocytes_DelPGC13_normoxia", "Astrocytes_DelP_normoxia"))
# Normally_distributed_data <-  na.omit(Normally_distributed_data)


# Normally_distributed_data$Name <- factor(Normally_distributed_data$Name, levels = c("Astrocytes_wt_normoxia", "Astrocytes_A8mut_normoxia","Astrocytes_DelPGC13_normoxia", "Astrocytes_DelP_normoxia", "Astrocytes_wt_hypoxia", "Astrocytes_A8mut_hypoxia","Astrocytes_DelPGC13_hypoxia", "Astrocytes_DelP_hypoxia"))


                                         
                                         # 
#                                          # levels = c("Astrocytes_DelP_normoxia", "Astrocytes_DelPGC13_normoxia", "Astrocytes_DelPGC8_normoxia", "Astrocytes_A8mut_normoxia", "Astrocytes_wt_normoxia", "Astrocytes_DJ1OE_normoxia"))
# #obtain a list containing a tibble for each metabolite with the corresponding data
# 
# #obtain a list containing a tibble for each metabolite with the corresponding data
# lDf <- split(Normally_distributed_data, Normally_distributed_data$metabolite)
# #ANOVA comparing the groups
# pwc <- map(lDf, ~.x  %>%
#              #group_by(Mass) %>%
#              emmeans_test(Mean_Technical_Replicates ~ Name, p.adjust.method = "bonferroni"))
# #pwc

  
```

```{r Exporting of data}
#calculate the mean of the biological replicates to plot the mean of the biological replicates 
Table_Mean_Biological_Replicates <- Table_Mean_Technical_Replicates %>%
  group_by(metabolite, Mass, Name) %>%
  mutate("Mean_Biological_Replicates" = mean(Mean_Technical_Replicates),
         "Stdev_Biological_Replicates" = sd(Mean_Technical_Replicates), 
         "SEM_Biological_Replicates" = sd(Mean_Technical_Replicates)/sqrt(length(Mean_Technical_Replicates))) %>%
  #mutate(deviation2 = ifelse(Stdev_Biological_Replicates >= 0.001, TRUE, FALSE)) %>%
  #filter(deviation2 == FALSE) %>%
  #keep only one value for each mean of the technical replicates and shrink the table by using the distinct function
  distinct(metabolite, Mass, Name, condition, .keep_all = TRUE)

#write_csv(Table_Mean_Biological_Replicates, "C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/IGC metabolomics/Astrocytes paper/Astrocytes_intracellular_results_mean_biological_rep_n1_5_normoxia_hypoxia.csv")

FinalTable <- Table_Mean_Biological_Replicates %>%  
  filter(condition == "normoxia")
FinalTable$Name <- factor(FinalTable$Name, levels = c("Astrocytes_wt_normoxia", "Astrocytes_A8mut_normoxia","Astrocytes_DelPGC13_normoxia", "Astrocytes_DelP_normoxia"))
FinalTable <- na.omit(FinalTable)

#FinalTable[- grep("REVERSE", FinalTable$Name),]

#FinalTable$Name <- factor(FinalTable$Name, levels = c("Astrocytes_wt_normoxia", "Astrocytes_A8mut_normoxia","Astrocytes_DelPGC13_normoxia", "Astrocytes_DelP_normoxia", "Astrocytes_wt_hypoxia", "Astrocytes_A8mut_hypoxia","Astrocytes_DelPGC13_hypoxia", "Astrocytes_DelP_hypoxia"))
# 
# 
FinalTable <- FinalTable %>%  
mutate(deviation2 = ifelse(Mean_Biological_Replicates <=0, TRUE, FALSE)) %>%
filter(deviation2 == FALSE) %>%  
  na.omit(FinalTable)
#write_csv(FinalTable, "C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/IGC metabolomics/Astrocytes paper/Astrocytes_intracellular_results_mean_biological_rep_n1_5_new_outlier_test.csv")


# sem <- function(x) sd(x)/sqrt(length(x))
   
# FinalTable <- FinalTable
# FinalTable$Name <- factor(FinalTable$Name, levels = c("Astrocytes_DelP_normoxia", "Astrocytes_DelPGC13_normoxia", "Astrocytes_DelPGC8_normoxia", "Astrocytes_A8mut_normoxia", "Astrocytes_wt_normoxia", "Astrocytes_DJ1OE_normoxia"))
                                                      
                                                      
                                                      #"Astrocytes_DJ1OE_normoxia"))
                                                      
                                                      #, "Astrocytes_DelP_hypoxia", "Astrocytes_DelPGC13_hypoxia", "Astrocytes_A8mut_hypoxia", "Astrocytes_wt_hypoxia", "Astrocytes_DJ1OE_hypoxia"))


```



```{r Automated image generator - Plotting of data} {r automated plots including stats, fig.height= 8, fig.width= 15}


#define where you want to have your plots
DataDirectory <- "C:/Users/pauline.mencke/Documents/Pauline/Lab/1. Main paper thesis/Astrocytes/metabolite tracing/GCMS/IGC Astrocytes paper/all_ns_MIDs/"
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
  #theme(legend.text=element_text(size=20))+
  theme(axis.title.x = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  # theme(axis.text.x=element_blank())+
  #labs(x = "M0")+
  #theme(legend.title=element_blank())+
  theme(legend.position = "none")+


  #set colours
  scale_fill_manual(values=c("chartreuse3", "blue2", "darkolivegreen2", "darkslategray2"))
  #scale_fill_manual(values=c("green", "red", "blue", "orange", "green", "red", "blue", "orange")) #+
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
    mutate(Ymin = (Mean_Biological_Replicates-SEM_Biological_Replicates),
    Ymax = (Mean_Biological_Replicates+SEM_Biological_Replicates)) #%>%
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


```{r Automated image generator - Plotting of data} {r automated plots including stats, fig.height= 8, fig.width= 15}


#define where you want to have your plots
DataDirectory <- "C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/IGC metabolomics/Astrocytes paper/all_ns_MIDs_hypoxia_normoxia/"
#define the graph function
 make_graph <- function(Data, x, DataDirectory){#, stat.test){
   #define the metabolite name and set a name for the file
   metabolite_name <- x
   name <- paste(DataDirectory, metabolite_name, ".png", sep="")
   #open image creating part of the function, set resolution
   png(name, width = 720, height = 720, units = "px")
   #make the graph
   Graph <-

     ggplot(Data, aes(x = Name, y = Mean_Biological_Replicates, fill = Name, pattern = condition)) +
     geom_col_pattern(position = position_dodge(width = 0.8), stat = "identity",
                      #pattern = c("none", "none", "none" , "none", "stripe", "stripe", "stripe", "stripe"),
                      #width = 0.75,
                      color = "black",
                      pattern_fill = "black",
                      pattern_angle = 45,
                      pattern_density = 0.1,
                      pattern_spacing = 0.025,
                      pattern_key_scale_factor = 0.6) +
                      #mapping=aes(pattern=condition) +
      scale_pattern_manual(values = c(hypoxia = "stripe", normoxia = "none")) +
      guides(pattern = guide_legend(override.aes = list(fill = "white")),
          fill = guide_legend(override.aes = list(pattern = "none"))) +
     theme(axis.text.x = element_blank()) +
     labs(x = "M0") +
     ggtitle(metabolite_name) +
     #theme_classic() +
     theme(plot.title = element_text(hjust = 0.5)) +

   #set colours
   scale_fill_manual(values=c("green", "red", "blue", "orange", "green", "red", "blue", "orange")) +
   # # #scale_fill_manual(values=c("orange", "blue", "cyan", "red", "green", "yellow"))+
   # #                            #, "red", "green", "yellow", "red", "green", "yellow", "red", "green", "yellow")) +


   geom_errorbar(aes(ymin = Ymin, ymax = Ymax), position=position_dodge2(width = 0.9), colour="black") +
   geom_errorbar(mapping=aes(ymin = Ymin, ymax = Ymax, fill = Name), position=position_dodge(width=0.8), colour="black")
   #geom_point(position=position_dodge(width=0.8)) +
   #add statistics
   #stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.001, bracket.nudge.y = 0, dodge = 0.8, hide.ns = TRUE)
   #by printing, the graph is saved at the data directory
   print(Graph)
   #the png module is closed
   dev.off()
 }
#
# #run the for loop
 for(metabolite_name in unique(FinalTable$metabolite)){
#
     Data <- FinalTable %>%
     filter(metabolite == metabolite_name) %>%
     mutate(Ymin = (Mean_Biological_Replicates-SEM_Biological_Replicates),
     Ymax = (Mean_Biological_Replicates+SEM_Biological_Replicates)) #%>%
#     #
#     #group_by(Mass)
#
#     # pwc[[metabolite_name]] %>%
#     # add_xy_position(x = "Name", fun = "max",  step.increase = 0.05, dodge = 0.8) -> stat.test
#
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
#


``` 





