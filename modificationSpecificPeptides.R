

library(tidyverse)
library(dendextend)
library(cluster)
library(RColorBrewer)
library(gplots)
library(janitor)
library(rstatix)
library(ggpubr)
library(dplyr)
library(data.table)
library(devtools)
library(ggbiplot)
library(ggplot2)

setwd("C:/Users/A.Kotronoulas/Dropbox (Volition)/Aris/!Science/!Projects/CNC proteomics/EpiQMAx IP nov 2019/Epimax 2019 multi")

#### Prepare dataset

# read csv
mod_pep<-read_csv("modificationSpecificPeptides.csv")
metadata_mod_pep<-read_csv("metadata_mod_pep.csv")

# select columns of intensities
mod_pep<-mod_pep[, c(3:15)]

# In order to find the average between the 2 technical replicates, I convert all zeros to NAs.
mod_pep[mod_pep == 0] <- NA

# Find average
mod_pep$Sample_1<-(mod_pep$`Intensity 1_1`+mod_pep$`Intensity 1_2`)/2
mod_pep$Sample_2<-(mod_pep$`Intensity 2_1`+mod_pep$`Intensity 2_2`)/2
mod_pep$Sample_3<-(mod_pep$`Intensity 3_1`+mod_pep$`Intensity 3_2`)/2
mod_pep$Sample_4<-(mod_pep$`Intensity 4_1`+mod_pep$`Intensity 4_2`)/2
mod_pep$Sample_5<-(mod_pep$`Intensity 5_1`+mod_pep$`Intensity 5_2`)/2
mod_pep$Sample_6<-(mod_pep$`Intensity 6_1`+mod_pep$`Intensity 6_2`)/2

# Select only the columns of average for each peptide
mod_pep<-mod_pep[,c(1, 14:19)]


#### From now on I will analyse the data as follows for each group (CRC or control):
#       1. Peptides that have NAs for all samples will be turned to zero
#       2. Peptides with one NA (up to 33%) will be imputed
#       3. Peptides with more than one NA (more than 33%) but not all NA will be discarded   


# For CRC samples:
# select CRC samples 
mod_pep_CRC<-mod_pep[,c(1:4)]
# select peptides with all NA and replace with 0
cnt_na_CRC <- apply(mod_pep_CRC, 1, function(z) sum(is.na(z)))
mod_pep_CRC_All_na<-mod_pep_CRC[cnt_na_CRC > 2,]
mod_pep_CRC_All_na[is.na(mod_pep_CRC_All_na)] <- 0
# select peptides with one NA and imp ute
mod_pep_CRC_one_na<-mod_pep_CRC[cnt_na_CRC < 2,]
mod_pep_CRC_one_na$Sample_1[is.na(mod_pep_CRC_one_na$Sample_1)] <- 
  rowMeans(mod_pep_CRC_one_na[2:4], na.rm = TRUE)[is.na(mod_pep_CRC_one_na$Sample_1)] 
mod_pep_CRC_one_na$Sample_2[is.na(mod_pep_CRC_one_na$Sample_2)] <- 
  rowMeans(mod_pep_CRC_one_na[2:4], na.rm = TRUE)[is.na(mod_pep_CRC_one_na$Sample_2)] 
mod_pep_CRC_one_na$Sample_3[is.na(mod_pep_CRC_one_na$Sample_3)] <- 
  rowMeans(mod_pep_CRC_one_na[2:4], na.rm = TRUE)[is.na(mod_pep_CRC_one_na$Sample_3)] 
# bind df
mod_pep_CRC_final<-rbind(mod_pep_CRC_All_na, mod_pep_CRC_one_na)
mod_pep_CRC_final[2:4] <- lapply(mod_pep_CRC_final[2:4], as.numeric)


# For CONTROL samples:
# select control samples
 mod_pep_control<-mod_pep[,c(1, 5:7)]
 # select peptides with all NA and replace with 0
 cnt_na_control <- apply(mod_pep_control, 1, function(z) sum(is.na(z)))
 mod_pep_control_All_na<-mod_pep_control[cnt_na_control > 2,]
 mod_pep_control_All_na[is.na(mod_pep_control_All_na)] <- 0
 # select peptides with one NA and impute
 mod_pep_control_one_na<-mod_pep_control[cnt_na_control < 2,]
 mod_pep_control_one_na$Sample_4[is.na(mod_pep_control_one_na$Sample_4)] <- 
   rowMeans(mod_pep_control_one_na[2:4], na.rm = TRUE)[is.na(mod_pep_control_one_na$Sample_4)] 
 mod_pep_control_one_na$Sample_5[is.na(mod_pep_control_one_na$Sample_5)] <- 
   rowMeans(mod_pep_control_one_na[2:4], na.rm = TRUE)[is.na(mod_pep_control_one_na$Sample_5)] 
 mod_pep_control_one_na$Sample_6[is.na(mod_pep_control_one_na$Sample_6)] <- 
   rowMeans(mod_pep_control_one_na[2:4], na.rm = TRUE)[is.na(mod_pep_control_one_na$Sample_6)] 
 # bind df
 mod_pep_control_final<-rbind(mod_pep_control_All_na, mod_pep_control_one_na)
 mod_pep_control_final[2:4] <- lapply(mod_pep_control_final[2:4], as.numeric)
 

### Join the 2 dfs from CRC and CONTROL 
 df_final<-inner_join(mod_pep_CRC_final,mod_pep_control_final, by="Unique Sequence" )

#### From what I observed in the final df, I can split the peptides into the following categories:
#          1. Peptides that are !0 in CRC and 0 in all Control samples (CRC_Unique_peptides) 
#          2. Peptides that are !0 in Control and 0 in all CRC samples (Control_Unique_peptides)
#          3. Peptides that are !0 in all CRC and Control (Common_peptides)
#          4. Peptides that are 0 in all CRC and Control (zero_peptides)
 
CRC_Unique_peptides<-df_final %>%  filter(Sample_4==0&Sample_5==0&Sample_1!=0)
Control_Unique_peptides<-df_final %>%  filter(Sample_1==0&Sample_2==0&Sample_4!=0)
Common_peptides<-df_final%>%  filter(Sample_1!=0&Sample_2!=0&Sample_4!=0)
zero_peptides<-df_final %>%  filter(Sample_4==0&Sample_5==0&Sample_1==0)


### From now on I will do stat analysis only for the common peptides

Common_peptides<-column_to_rownames(Common_peptides, var = "Unique Sequence")
Common_peptides_log<-log(Common_peptides)
Common_peptides_log<-rownames_to_column(Common_peptides_log, var = "Unique Sequence")


hist(Common_peptides_log$Sample_1, breaks=30, col=rgb(0.3,0,0,0.5), xlab="log(intensity)", 
     ylab="nm of peptides", main="distribution of peptide log(intensity) in each sample" )

# Second with add=T to plot on top
hist(Common_peptides_log$Sample_2, breaks=30, xlim=c(0,300), col=rgb(0.6,0,0,0.5), add=T)
# third with add=T to plot on top
hist(Common_peptides_log$Sample_3, breaks=30, xlim=c(0,300), col=rgb(1,0,0,0.5), add=T)

hist(Common_peptides_log$Sample_4, breaks=30, xlim=c(0,300), col=rgb(0,0.3,0,0.5), add=T)
hist(Common_peptides_log$Sample_5, breaks=30, xlim=c(0,300), col=rgb(0,0.3,0,0.5), add=T)
hist(Common_peptides_log$Sample_3, breaks=30, xlim=c(0,300), col=rgb(0,0.3,0,0.5), add=T)



################# PCA #####################

# transpose

t_Common_peptides_log <- as.data.frame(t(Common_peptides_log))
t_Common_peptides_log <-t_Common_peptides_log %>%  row_to_names(row_number = 1)
t_Common_peptides_log<-rownames_to_column(t_Common_peptides_log, var = "sample")


# Assign Groups
t_Common_peptides_log<-t_Common_peptides_log %>% mutate(Group = 
                                      case_when(t_Common_peptides_log$sample == c("Sample_1","Sample_2", 
                                                                        "Sample_3") ~ 'CRC',
                                                t_Common_peptides_log$sample == c("Sample_4","Sample_5", 
                                                                        "Sample_6") ~ 'Control' )) %>%  group_by(Group)

t_com_pep_log<-t_Common_peptides_log[, c(1, 187, 2:186)]
t_com_pep_log[3:187] <- sapply(t_com_pep_log[3:187], as.numeric)
t_com_pep_log<-column_to_rownames(t_com_pep_log, var = "sample")


# then
t_com_pep_log.pca <- prcomp(t_com_pep_log[,-1], center = TRUE,scale. = TRUE)

fviz_eig(t_com_pep_log.pca) #Visualize eigenvalues (screen plot).

t_com_pep_log.scores<- t_com_pep_log.pca$x#this will give you the principal components and sd
summary(t_com_pep_log.pca)

t_com_pep_log.loadings<-t_com_pep_log.pca$rotation #this will give you the loadings
print(t_com_pep_log.loadings)

t_com_pep_log.scores<- t_com_pep_log.pca$x #this will give you the scores
print(t_com_pep_log.scores)

ggbiplot(t_com_pep_log.pca)

ggbiplot(t_com_pep_log.pca,ellipse=TRUE,obs.scale = 1, var.scale = 1,var.axes=FALSE, 
         groups=t_com_pep_log$Group )+
  scale_colour_manual(name="Origin", values= c("forest green", "red3"))+
  theme_minimal()

