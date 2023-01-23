---
title: "Extracellular metabolomics Astrocytes Glc"
author: "Pauline Mencke"
date: "4/22/2021"
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

```

```{r Loading of Data}
setwd('C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/R analysis metabolomics/Extracellular/')

#load data table 
TableImport = readxl::read_excel('C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/Extracellular metabolomics/Astrocytes/raw data Astrocytes paper GlcC13_N1_5/20210422_LP564_PMe18_MGC_results_formatted.xlsx', range = "A1:FD17") 


#melt the table so that each value has a unique identifier (each row is unique instead of heaving the cell line name as column names)
TableImport <- melt(TableImport, id.vars=c("Metabolite")) 
Table_V1 <- TableImport %>% 
  rename(Full_Condition_Name = variable)

#Adjust the cell line and condition name so that it contains only the info needed 
Table_V2 <- Table_V1 %>% 
  separate(Full_Condition_Name, into = c("V1", "V2", "V3", "V4", "V5", "V6", "V7"), sep = "_") %>% 
  mutate(Full_Condition_Name = glue::glue("{V3}_{V4}_{V6}")) %>%  
  mutate(condition = glue::glue("{V6}"))


Table_Astrocytes <- Table_V2[grep("^Astrocytes", Table_V2$Full_Condition_Name), ] 
Table_Astrocytes <- Table_Astrocytes %>%  
  mutate(Biological_Replicate = case_when(V5 == "1"~ 1,
                                         V5 == "2"~ 2,
                                         V5 == "3"~ 3,
                                         V5 == "4"~ 4,
                                         V5 == "5"~ 5)) %>%  
  mutate(Technical_Replicate = case_when(V7 == "1"~ 1,
                                         V7 == "2"~ 2,
                                         V7 == "3"~ 3)) %>% 
  rename(Name = Full_Condition_Name) %>%  
  select(Metabolite, Name, condition, Biological_Replicate, Technical_Replicate, value)


Table_Astrocytes_Blank <- Table_V2[grep("^FreshMedium", Table_V2$V3), ] 
Table_Astrocytes_Blank <- Table_Astrocytes_Blank %>%  
  select(Metabolite, value) %>%  
  rename(Fresh_value = value)
Table_Astrocytes_Blank$Name = c("Astrocytes_Fresh")

Table_Astrocytes_Standard <- Table_V2[grep("^STD", Table_V2$Full_Condition_Name), ]   
Table_Astrocytes_Standard <- Table_Astrocytes_Standard %>%  
  mutate(Name = glue::glue("{V3}_{V4}")) %>% 
  mutate(Technical_Replicate = V5) %>%  
  select(Metabolite, Name, Technical_Replicate, value)




#before we start with the calculation for the different cell types, we need to plot the standard curves for the measured peak areas of each metabolite and the corresponding concentration and determine the linear regression equation and evaluate the R squared value. Then, we can use the equation from the respective standard curve of each metabolite to calculate the concentration of each metabolite in the medium using the peak areas for each sample 


#calculate the mean of the technical triplicates for the standards 
Table_Standards_means <- Table_Astrocytes_Standard %>% 
  group_by(Metabolite, Name) %>% 
  mutate("Mean" = mean(value),
         "Stdev" = sd(value)) %>% 
  distinct(Metabolite, Name, value, Mean, Stdev)

#create a column that determines whether there are outliers in the data
Table_Standards_outliers <- Table_Standards_means %>%
  mutate(CV=Stdev/Mean) %>%
  mutate(deviation = ifelse(CV >= 1, TRUE, FALSE))

#create a table with no deviated data
Table_Standards_no_outliers <- Table_Standards_outliers %>%
  group_by(Metabolite, Name) %>%
  filter(deviation == FALSE)

#create a table with deviated data
Table_Standards_outliers <- Table_Standards_outliers %>%
  #keep only the rows with CV >= 1
  filter(deviation == TRUE)
#remove outliers in deviated data
#find Q1, Q3, and interquartile range for values in column A
Q1 <- quantile(Table_Standards_outliers$value, .25)
Q3 <- quantile(Table_Standards_outliers$value, .75)
IQR <- IQR(Table_Standards_outliers$value)

#only keep rows in dataframe that have values within 1.5*IQR of Q1 and Q3
Table_Standards_outliers_removed <- subset(Table_Standards_outliers, Table_Standards_outliers$value > (Q1 - 1.5*IQR) & Table_Standards_outliers$value < (Q3 + 1.5*IQR))


#bind data to obtain final df that contains only non deviated data
Table_Standards_without_outliers <- rbind(Table_Standards_no_outliers, Table_Standards_outliers_removed) %>%
  distinct(Metabolite, Name, Mean)



Table_Standards_without_outliers$Name <- factor(Table_Standards_without_outliers$Name, levels = c("STD_10", "STD_20", "STD_40", "STD_60", "STD_80", "STD_100"))
                                     
#add info on protein level for normalization to protein level  
Table_Astrocytes <- Table_Astrocytes %>% 
  mutate(binding_Name = glue::glue("{Name}_{Technical_Replicate}_{Biological_Replicate}"))
Table_protein_quanti_Astrocytes = readxl::read_excel('C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/Protein quantification for normalization extracellular and IGC/Astrocytes/Astrocytes paper/final/Astrocytes protein quantification n1_5 31.03.2021.xlsx', range = "A1:B163")
Table_protein_quanti_Astrocytes <- Table_protein_quanti_Astrocytes %>% 
   mutate(binding_Name = Name)
p <- left_join(Table_Astrocytes, Table_protein_quanti_Astrocytes, by = "binding_Name")
Table_final_Astrocytes <- p %>%
  dplyr::select(Metabolite, Name.x, condition, Technical_Replicate, Biological_Replicate, value, protein_well) %>% 
  rename(Name = Name.x)

                                     
                                     

```


```{r linear regression and calculation of mol per metabolite - trying to do it for multiple metabolites }
#Formatting table of standards for calculations 
Table_Standards_mM = readxl::read_excel('C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/Extracellular metabolomics/Astrocytes/raw data Astrocytes paper GlcC13_N1_5/Standards_Astrocytes paper.xlsx', range = "A1:C115")
Table_Standards_peak_areas <- Table_Standards_without_outliers %>%
  mutate(binding_name = glue::glue("{Metabolite}_{Name}"))
Table_Standards_mM <- Table_Standards_mM %>%
   mutate(binding_name = glue::glue("{Metabolite}_{Name}"))
p <- left_join(Table_Standards_peak_areas, Table_Standards_mM, by = "binding_name")
data <- p %>%
  dplyr::select(Metabolite.x, Name.x, Mean, Standard_conc_in_mM) %>%
  rename(Name = Name.x) %>%
  rename(Metabolite = Metabolite.x)


#calculate linear regression 
lin_reg <- data.frame(NULL)
for(metabolite_name in unique(data$Metabolite)){
  m <- data %>% 
    filter(Metabolite == metabolite_name)  
  m <- summary(lm(m$Mean~m$Standard_conc_in_mM))    # run model
  lin_reg[metabolite_name, 1] <- names(data)[metabolite_name]           # print variable name
  lin_reg[metabolite_name, 2] <- m$coefficients[1,1]   # intercept
  lin_reg[metabolite_name, 3] <- m$coefficients[2,1]
  lin_reg[metabolite_name, 4] <- m$r.squared
                                    
}
names(lin_reg) <- c("y.variable", "intercept", "coef.x", "r.squared")
lin_reg <- lin_reg %>%  
  mutate(Metabolite = row.names(lin_reg))
lin_reg$y.variable = NULL


#"y = intercept + coef.x * x"
#we want to know x (mM on x axis) 
#x = (y+intercept)/coef.x Q
#y are the peak areas from our samples 

Table_for_all <- left_join(Table_final_Astrocytes, lin_reg, by = "Metabolite") 
FinalTable <- left_join(Table_for_all, Table_Astrocytes_Blank, by = "Metabolite") %>% 
  rename(Name = Name.x) 
FinalTable$Name.y = NULL 
FinalTable <- FinalTable %>%  
  mutate(conc_mMol = ((value+intercept)/coef.x )) %>% 
  mutate(conc_Mol = conc_mMol/1000 ) %>%  
  mutate(Blank_conc_mMol = ((Fresh_value+intercept)/coef.x )) %>% 
  mutate(Blank_conc_Mol = Blank_conc_mMol/1000 ) %>%  
#normalize to blank (fresh medium)
  mutate(Delta_conc_Mol_fresh_Blank = abs(conc_Mol - Blank_conc_Mol)) %>% 
  mutate(Delta_percentage_from_Blank = abs((Delta_conc_Mol_fresh_Blank/Blank_conc_Mol)*100)) %>%  
  mutate(Number_of_Moles = Delta_conc_Mol_fresh_Blank * 0.001) %>% 
#normalize to protein 
  mutate(Number_of_Moles_per_protein = Number_of_Moles/protein_well) %>% 
#normalize to incubation with tracer (per default 48hours), to calculate uptake and release rates in mol per hour 
  mutate(Moles_per_protein_per_hour = Number_of_Moles_per_protein/48 ) %>% 
#calculate fmol 
  mutate(fmol_protein_hour = Moles_per_protein_per_hour*1000000000000000) %>%  
  mutate(abs_fmol_protein_hour = abs(Moles_per_protein_per_hour*1000000000000000))

write_csv(FinalTable, "C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/Extracellular metabolomics/Astrocytes/output files Astrocytes paper/Astrocytes_paper_extracellular_results.csv")

FinalTableFiltered <- FinalTable %>%  
 filter(condition == "normoxia")
 #filter(Delta_percentage_from_Blank >= 10)

write_csv(FinalTableFiltered, "C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/Extracellular metabolomics/Astrocytes/output files Astrocytes paper/Astrocytes_paper_extracellular_results_filtered.csv")

FinalTableFiltered <- FinalTableFiltered %>%
  group_by(Metabolite, Name, Biological_Replicate) %>%  
  mutate("Mean_Technical_Replicates" = abs(mean(Delta_conc_Mol_fresh_Blank)),
         "Stdev_Technical_Replicates" = abs(sd(Delta_conc_Mol_fresh_Blank))) %>%
  
  #keep only one value for each mean of the technical replicates and shrink the table by using the distinct function
  select(Metabolite, Name, condition, Biological_Replicate, Mean_Technical_Replicates, Stdev_Technical_Replicates) %>%  
  #mutate_all(~replace(., is.na(.), 0)) %>%
  distinct(Metabolite, Name, condition, Biological_Replicate, .keep_all = TRUE) %>%  
  filter(condition == "normoxia")

FinalTableFiltered <- FinalTableFiltered %>%
  group_by(Metabolite, Name) %>%
  mutate("Mean_Biological_Replicates" = mean(Mean_Technical_Replicates),
         "Stdev_Biological_Replicates" = sd(Stdev_Technical_Replicates)) %>% 

#FinalTableFiltered <- na.omit(FinalTableFiltered) %>%  
  select(Metabolite, Name, condition, Biological_Replicate, Mean_Biological_Replicates, Stdev_Biological_Replicates) %>%  
  distinct(Metabolite, Name, condition, .keep_all = TRUE)

target = c("Astrocytes_wt_normoxia", "Astrocytes_A8mut_normoxia","Astrocytes_DelPGC13_normoxia", "Astrocytes_DelP_normoxia")
FinalTableFiltered <- filter(FinalTableFiltered, Name %in% target)

```


```{r linear regression graphs}

###########################################################################
#plotting of linear regression for multiple metabolites

for(metabolite_name in unique(data$Metabolite)){
  Data <- data %>%
    filter(Metabolite == metabolite_name)
  
 #define where you want to have your plots   
 mypath <- file.path("C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/Extracellular metabolomics/Astrocytes/plots Astrocytes paper/plots standards Astrocytes paper/",paste(metabolite_name, ".png", sep = ""))

 png(file=mypath)
    mytitle = paste(metabolite_name)
    plot(x = Data$Standard_conc_in_mM, y = Data$Mean, cex = 1.3, pch = 16, xlab = "conc in mM", ylab = "Mean peak areas", main = metabolite_name)
    #abline(lm(data$Mean ~ data$Standard_conc_in_mM))
    
     
 dev.off()
 
}




###########################################################################
```

```{r Statistics - Normality test and ANOVA multiple comparisons}
#calculate the mean of the technical replicates
Table_Mean_Technical_Replicates <- FinalTable %>%
  group_by(Metabolite, Name, Biological_Replicate) %>%  
  mutate("Mean_Technical_Replicates" = abs(mean(fmol_protein_hour)),
         "Stdev_Technical_Replicates" = abs(sd(fmol_protein_hour))) %>%
  #keep only one value for each mean of the technical replicates and shrink the table by using the distinct function
  select(Metabolite, Name, condition, Biological_Replicate, Mean_Technical_Replicates, Stdev_Technical_Replicates) %>%  
  #mutate_all(~replace(., is.na(.), 0)) %>%
  distinct(Metabolite, Name, condition, Biological_Replicate, .keep_all = TRUE) %>%  
na.omit(Table_Mean_Technical_Replicates)
write_csv(Table_Mean_Technical_Replicates, "C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/Extracellular metabolomics/Astrocytes/output files Astrocytes paper/Astrocytes_paper_extracellular_results_mean_technical_rep.csv")
#perform normality test of data before pushing the data into ANOVA or dunn test
Data_normality_test <- Table_Mean_Technical_Replicates %>%
  group_by(Metabolite, Name) %>%
  normality(Mean_Technical_Replicates) %>%
  #keep only the values where p_value is greater than 0.01
  filter(p_value > 0.01)
#obtain the values for the normally distributed data
Normally_distributed_data <- Table_Mean_Technical_Replicates %>%
  filter(Metabolite %in% Data_normality_test$Metabolite) %>%
  ungroup() %>%  
  filter(condition == "normoxia")
target = c("Astrocytes_wt_normoxia", "Astrocytes_A8mut_normoxia","Astrocytes_DelPGC13_normoxia", "Astrocytes_DelP_normoxia")
Normally_distributed_data <- filter(Normally_distributed_data, Name %in% target)
# #obtain a list containing a tibble for each metabolite with the corresponding data
lDf <- split(Normally_distributed_data, Normally_distributed_data$Metabolite)
#ANOVA comparing the groups
pwc <- map(lDf, ~.x  %>%
           emmeans_test(Mean_Technical_Replicates ~ Name, p.adjust.method = "bonferroni"))
#pwc

  
FinalTableP <- Normally_distributed_data %>%  
  filter(condition == "normoxia") %>%  
  filter(Metabolite == "Glucose")


stat.test <- FinalTableP  %>% 
  t_test(Mean_Technical_Replicates ~ Name, paired = TRUE) %>%
  add_significance()
stat.test



```


```{r Exporting of data}
#calculate the mean of the biological replicates to plot the mean of the biological replicates 
Table_Mean_Biological_Replicates <- Normally_distributed_data %>%
  ungroup() %>%
  group_by(Metabolite, Name) %>%
  mutate("Mean_Biological_Replicates" = mean(Mean_Technical_Replicates),
     "Stdev_Biological_Replicates" = sd(Mean_Technical_Replicates))

FinalTableP <- Normally_distributed_data %>%  
 filter(condition == "normoxia") %>%  
 filter(Metabolite == "Glucose") %>%  
 ungroup() %>%
 group_by(Metabolite, Name) %>%
 mutate("Mean_Biological_Replicates" = mean(Mean_Technical_Replicates),
     "Stdev_Biological_Replicates" = sd(Mean_Technical_Replicates))


  #    %>%
  # distinct(Metabolite, Name, .keep_all = TRUE))
write_csv(Table_Mean_Biological_Replicates, "C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/Extracellular metabolomics/Astrocytes/output files Astrocytes paper/Astrocytes_paper_extracellular_results_mean_biological_rep.csv")
FinalTable <- Table_Mean_Biological_Replicates
target = c("Astrocytes_wt_normoxia", "Astrocytes_A8mut_normoxia","Astrocytes_DelPGC13_normoxia", "Astrocytes_DelP_normoxia")
FinalTable <- filter(FinalTable, Name %in% target)
FinalTable$Name <- factor(FinalTable$Name, levels = c("Astrocytes_wt_normoxia", "Astrocytes_A8mut_normoxia","Astrocytes_DelPGC13_normoxia", "Astrocytes_DelP_normoxia"))


```


# Automated image generator
```{r automated plots including stats, fig.height= 8, fig.width= 15}
#define where you want to have your plots
DataDirectory <- "C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/Extracellular metabolomics/Astrocytes/plots Astrocytes paper/Delta_Medium_vs_Fresh/"
#define the graph function
make_graph <- function(Data, x, DataDirectory, stat.test){
  #define the metabolite name and set a name for the file
  metabolite_name <- x  
  name <- paste(DataDirectory, metabolite_name, ".png", sep="")
  #open image creating part of the function, set resolution
  png(name, width = 720, height = 720, units = "px") 
  #make the graph
  Graph <-
  
  
    
  ggplot(Data, aes(x = Name, y = Mean_Biological_Replicates, col = Name))+
  geom_bar(aes(fill = Name), stat="identity", colour = "black", position=position_dodge(0.8))+  
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
  theme(axis.text.x=element_blank())+
  
  #theme(legend.title=element_blank())+
  
    # 
    # 
    # ggplot(Data, aes(x = Mass, y = Mean_Biological_Replicates, fill = Name, pattern = condition)) + 
    # geom_col_pattern(position = position_dodge(width = 0.9), stat = "identity",
    #                  #pattern = c("none", "none", "none" , "stripe", "stripe", "stripe"), 
    #                  width = 0.75, 
    #                  color = "black", 
    #                  pattern_fill = "black",
    #                  pattern_angle = 45,
    #                  pattern_density = 0.1,
    #                  pattern_spacing = 0.025,
    #                  pattern_key_scale_factor = 0.6) +
    #                  #mapping=aes(pattern=condition)) +
    #  scale_pattern_manual(values = c(hypoxia = "stripe", normoxia = "none")) +
    #  guides(pattern = guide_legend(override.aes = list(fill = "white")),
    #      fill = guide_legend(override.aes = list(pattern = "none"))) +
    # 
   
  #set colours
  scale_fill_manual(values=c("green", "red", "blue", "orange"))+
                             #, "red", "green", "yellow", "red", "green", "yellow", "red", "green", "yellow")) +
           
  geom_errorbar(mapping=aes(ymin = Ymin, ymax = Ymax, fill = Name), position=position_dodge(width=0.8), colour="black")# +
  #geom_point(position=position_dodge(width=0.8))+
  # ggtitle(metabolite_name)+
  # theme_classic() +
  # theme(plot.title = element_text(hjust = 0.5))#+
    
    
  #add statistics
  #stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.001, bracket.nudge.y = 0, dodge = 0.8, hide.ns = TRUE)
  
  
  
  #by printing, the graph is saved at the data directory   
  print(Graph)
  #the png module is closed
  dev.off()
  
}


#run the for loop 
for(metabolite_name in unique(FinalTableFiltered$Metabolite)){
  Data <- FinalTableFiltered %>% 
    filter(Metabolite == metabolite_name) %>% 
    mutate(Ymin = (Mean_Biological_Replicates-Stdev_Biological_Replicates),
         Ymax = (Mean_Biological_Replicates+Stdev_Biological_Replicates)) 
    # pwc[[metabolite_name]] %>%
    # add_xy_position(x = "Name", fun = "max",  step.increase = 0.05, dodge = 0.8)->stat.test
    
  make_graph(Data, metabolite_name, DataDirectory, stat.test = stat.test)
}

# 
# 
# 
# 
# setwd("C:/Users/pauline.mencke/Documents/Pauline/Lab/Metabolite tracing/Extracellular metabolomics/plots Astrocytes/")
# pdf('Summary_all_metabolites_extracellular_astro_Glc.pdf')
# lapply(ll <- list.files(patt='.*[.]png'),function(x){
#   img <- as.raster(readPNG(x))
#   grid.newpage()
#   grid.raster(img, interpolate = FALSE)
# 
# })
# dev.off()
# 

```




