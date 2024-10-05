library(readxl)
library(readr)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(tidyr)
library(iNEXT)
library(VennDiagram)
library(RColorBrewer)
library(devtools)
library(tibble)
library(purrr)

#Permisson requried install_github not working #Note, deleted by Feng
#install_github("https://github.com/BULQI/LungVirome/tree/Fumi/cha")


setwd("~/Datasets/Archive_Code/Partition_Analysis/")
source("Sequencing_Function_Codes.R")

wd <- "C:/Users/faihara/Desktop/faihara/Desktop/Lab/Sequencing/iNEXT/PartitionFIles"

wd <- "~/Datasets/Archive_Code/Partition_Analysis"
setwd(wd)

#Read all Partition assignments files
Raw_Data <- list()

for (x in 1:length(list.files(pattern = "PartitionAssignments"))) {
#Read PartitionAssignmentAnalysis results from Cloanalyst
    Raw_Data[[x]] <- read_delim(
      list.files(pattern = paste0("PartitionAssignmentsAnalysis_P",x
                                  )
                 ),
                              delim = "\t",
                              escape_double = FALSE,
                              trim_ws = TRUE
      )
    
    #Assign names of sample to proper list
  names(Raw_Data)[[x]] <- as.character(paste0("P", x))
}

#Collect total number of clones for each QID
QID_Counts <- list()

for (x in 1:length(Raw_Data)) {
  
  QID_Counts[[x]] <- Raw_Data[[x]] %>%
    group_by(QID) %>%
    dplyr::count(QID) %>%
    dplyr::rename(Counts = n)
  
  names(QID_Counts)[[x]] <- as.character(paste0("P", x))
}

#Collect total number of clones in a QID by sample----
Sample_Counts <- list()
for (x in 1:length(Raw_Data)) {
  
  Sample_Counts[[x]] <- Raw_Data[[x]] %>%
    group_by(SampleID) %>%
    count(QID) %>%
    rename(Counts = n)
  
  names(Sample_Counts)[[x]] <- as.character(paste0("P", x))
}

#Cha function from QID_Counts----
#QID_Counts wide format to prepare for Cha, for single sample
QID_Counts_Wide <- pivot_wider(Sample_Counts$P1,
                               values_from = Counts,
                               names_from = SampleID,
                               values_fill = 0
                               )

#QID_Counts wide format for Cha, All sample
QID_Counts_Wide_All <- list()

for (x in 1:length(QID_Counts)) {
  
  QID_Counts_Wide_All[[x]] <- pivot_wider(Sample_Counts[[x]],
                                          values_from = Counts,
                                          names_from = SampleID,
                                          values_fill = 0)
  
  names(QID_Counts_Wide_All)[[x]] <- as.character(paste0("P", x))
}

colnames(c(QID_Counts_Wide_All[[1]],
         QID_Counts_Wide_All[[2]]))

#Create comparison list for loop
Comparison.table <- as.list()

#Attempt to recapituate HillRatio function
  combinations <- as.matrix(c())
  
#HillRatio: Calculate Hill ratios from iNEXT Hill numbers----
HR <- list()

#HillRatio functional code lost
#for (x in 1:length(QID_Counts_Wide)) {
#  matrix <- as.matrix(QID_Counts_Wide_All[[x]])
#  HR[[x]] <- HillRatio(matrix, Call_iNEXT = F)
  
#  names(HR)[[x]] <- as.character(paste0("P", x))
#}

#Manual entry to generate hill ratios
matrix <- as.matrix(QID_Counts_Wide_All[[x]])
HR[[x]] <- HillAnalysis(matrix[,2],matrix[,4], Call_iNEXT = T)

names(HR)[[x]] <- as.character(paste0("P", x))

##Condense to one df
HR_DF <- data.frame(do.call(rbind, HR))
colnames(HR_DF) <- c("Ratio1", "Ratio2", "Ratio3")

##Add Patient number column
HR_DF <- HR_DF %>%
  mutate(Patient_Num = row_number())

#Graph----
Plot_DF <- pivot_longer(HR_DF,
                        cols = starts_with("Ratio"),
                        names_to = "Diversity_Ratio")
Plot_DF$Patient_Num <- as.factor(Plot_DF$Patient_Num)


ggplot(Plot_DF, aes(x = Diversity_Ratio,
                    y = value,
                    color = Patient_Num)
       ) +
  
  geom_point(size = 4)

#Save graph
destination <- "C:/Users/faihara/Desktop/faihara/Desktop/Lab/Writings/Figures/Paired_Sample/iNEXT"

ggsave(paste0("Diversity_Ratio_iNEXT.png"),
       width = 5,
       height = 6.5,
       path = destination)
