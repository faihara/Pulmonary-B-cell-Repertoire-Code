library(readxl)
library(readr)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(tidyr)
library(iNEXT)

#Read partition data
setwd("~/Datasets/Archive_Code/Partition_Analysis/")
source("Sequencing_Function_Codes.R")

Raw_Data1 <- read_delim("2021-01-15_PartitionAssignmentsAnalysis_P1.txt",
           delim = "\t",
           escape_double = FALSE,
           trim_ws = TRUE
           )

#Downsample by smallest sample (size = 0 will downsample by smallest)
Raw_Down_P1 <- Downsample(Raw_Data1, size = 1000)

#Check downsample
Raw_Down_P1 %>%
  group_by(SampleID) %>%
  count(SampleID)

#Subset general information-------------------------
#Collect unique QIDs
QID_Unique <- Raw_Data1 %>%
  distinct(QID)

#Collect unique sample IDs
SampleID_Unique <- Raw_Data1 %>%
  distinct(SampleID)

#Number of q-clones in a sample
SampleID_Counts1 <- Raw_Data1 %>%
  group_by(SampleID) %>%
  dplyr::count(SampleID)

#Collect total number of clones in a QID
QID_Counts <- Raw_Data1 %>%
  group_by(QID) %>%
  dplyr::count(QID) %>%
  rename(Counts = n)

#Collect total number of clones in a QID by sample
Sample_Counts <- Raw_Data1 %>%
  group_by(SampleID) %>%
  dplyr::count(QID) %>%
  rename(Counts = n)


#Analysis------------------------------------
#Find how often a clone appears across samples (Degeneracy)
QID_Degeneracy <- Raw_Data1 %>%
  group_by(SampleID) %>%
  distinct(QID) %>%
  ungroup() %>%
  count(QID) %>%
  rename(Degeneracy = n)

#Collect total number of clones in a QID
QID_Counts <- Raw_Data1 %>%
  group_by(QID) %>%
  count(QID) %>%
  rename(Counts = n)

QID_Degeneracy <- left_join(QID_Degeneracy,
                            QID_Counts,
                            by = "QID") #add QIDs and tissue information

##Save data
write.table(QID_Degeneracy,
            "QID_Degeneracy.txt",
            sep = "\t",
            row.names = FALSE)

QID_Degeneracy_Count <- count(QID_Degeneracy,
                              Degeneracy) #Count by degeneracy



##Subset by tissue and degeneracy-----------------------------------------------

##Dictate keywords for selecting tissue (Required for DatabySample function)
SubsetKey <- c('B', "L")

##Calculate QID counts by tissue type
QID_Counts_Tissue <- DatabySample(Raw_Data1)
##Convert list to tibble
QID_Counts_LB <- tibble(do.call(rbind, QID_Counts_Tissue))

##Rename B to Blood
QID_Counts_LB$Sample <- str_replace(QID_Counts_LB$Sample,
                                    pattern = "B",
                                    replacement = "Blood")

##Rename L to Lung
QID_Counts_LB$Sample <- str_replace(QID_Counts_LB$Sample,
                                    pattern = "L",
                                    replacement = "Lung")

##Count Degeneracy
QID_Counts_LB_Degeneracy <- QID_Counts_LB %>%
  group_by(Sample) %>%
  distinct(QID) %>%
  ungroup() %>%
  count(QID) %>%
  rename(Degeneracy = n) %>%
  left_join(QID_Counts_LB, by = "QID")

##Subset shared QIDs
QID_Counts_LB_Shared <- QID_Counts_LB_Degeneracy[which(QID_Counts_LB_Degeneracy$Degeneracy == 2), ]

##Subset tissue unique
QID_Counts_LB_Unique <- QID_Counts_LB_Degeneracy[which(QID_Counts_LB_Degeneracy$Degeneracy == 1), ]

##Save data
write.table(QID_Counts_LB_Shared,
           "QID_Counts_LB_Shared.txt",
           sep = "\t",
           row.names = FALSE)

write.table(QID_Counts_LB_Unique,
            "QID_Counts_LB_Unique.txt",
            sep = "\t",
            row.names = FALSE)

#Plot Degeneracy
x <- c("Blood" = QID_Counts_LB_Degeneracy[which(QID_Counts_LB_Degeneracy$Sample == "Blood"),
                                          1],
       "Lung" = QID_Counts_LB_Degeneracy[which(QID_Counts_LB_Degeneracy$Sample == "Lung"),
                                         1]
       )


ggVennDiagram(x,
              label_alpha = 0,
              category.names = c("Blood", "Lung"),
              set_color = "black",
              lty = 1
              ) +
  ggplot2::scale_fill_gradient2(low = "#F4FAFE",
                                high = "#4981BF") #+
  #geom_sf(color = "black")



##QID analysis------------------

#Find top 10 unique QID
Top_LB_Unique <- Top_QID(Raw_Data)

##Save data
write.table(Top_LB_Unique,
            "Top_LB_Unique.txt",
            sep = "\t",
            row.names = FALSE)



##Entropy-----------------------
Blood <- filter(Raw_Data1,
                grepl(("B"),
                      SampleID
                      )
                )

Lung <- filter(Raw_Data1,
               grepl(("L"),
                     SampleID
                     )
               )

QID_Counts <- Blood %>%
  group_by(QID) %>%
  count(QID) %>%
  rename(Counts = n)

tmp <-QID_Unique <- QID_Counts %>%
  distinct(QID)

tmp2 <-QID_Unique <- QID_Counts %>%
  distinct(QID)

pb <- txtProgressBar(min = 0, max = nrow(QID_Counts), style = 3)

for (x in 1:nrow(QID_Counts)) {
  
  xi <- QID_Counts[x,2]
  n <- nrow(Blood)
  
  
  p_hat <- (xi*(xi-1))/(n*(n-1))
  #p_hat <- xi/n
  tmp[x,2] <- log(p_hat)
  tmp2[x,2] <- xi*log(p_hat)
  
  setTxtProgressBar(pb, x)
}

close(pb)

print(sum(tmp2$Counts))

for (x in 1:nrow(QID_Counts)) {
  
  xi <- QID_Counts[x,2]
  n <- nrow(Blood)
  
  
  p_hat <- (xi*(xi-1))/(n*(n-1))
  #p_hat <- xi/n
  tmp[x,2] <- log(p_hat)
  tmp2[x,2] <- n*p_hat*log(p_hat)
}


#Calculate Jaccard between tissue or within tissue-------------------
##Dictate keywords for selecting tissue (Required for DatabySample function)
Subset <- c(unique(Raw_Data$SampleID), "L", "B")


Jaccard <- tibble()

for (x in 1:length(Subset)) {
  SubsetKey <- Subset[x]
  
  for (y in 1:length(Subset)) {
    SubsetKey <- c(Subset[x], Subset[y])
    ##Calculate QID counts by tissue type
    QID_Counts_Subset <- DatabySample(Raw_Data)
    
    ##Convert list to tibble
    QID_Counts_Subset <- tibble(do.call(rbind, QID_Counts_Subset))
    
    ##Count Degeneracy
    QID_Counts_Subset_Degeneracy <- QID_Counts_Subset %>%
      group_by(Sample) %>%
      distinct(QID) %>%
      ungroup() %>%
      count(QID) %>%
      rename(Degeneracy = n) %>%
      left_join(QID_Counts_Subset, by = "QID")
    
    ##Subset shared QIDs
    QID_Counts_Subset_Shared <- QID_Counts_Subset_Degeneracy[which(QID_Counts_Subset_Degeneracy$Degeneracy == 2), ]
    
    #Jaccard_Index <- length of shared QIDs / length of all QIDs
    Shared_QID <- length(unique(QID_Counts_Subset_Shared$QID))
    All_QID <- length(unique(QID_Counts_Subset_Degeneracy$QID))
    
    #Calculate Jaccard
    Jaccard[x,y] <- Shared_QID/All_QID
  }
}
colnames(Jaccard) <- c(Subset)
Jaccard <- cbind(Subset,Jaccard)

write.table(Jaccard,
            "Jaccard.txt",
            sep = "\t",
            row.names = FALSE)

#write.table(Jaccard,"Jaccard_table.txt",sep = "\t",row.names = FALSE)

#filter QIDs more then 10
tmp <- filter(Jaccard, Total_Count > 10)

Jaccard_Over <- filter(Jaccard, Jaccard >= 50)

write.table(Jaccard_Over,
            "Jaccard_Over.txt",
            sep = "\t",
            row.names = FALSE)
Jaccard_Under <- filter(Jaccard, Jaccard <50)

write.table(Jaccard_Under,
            "Jaccard_Under.txt",
            sep = "\t",
            row.names = FALSE)



##Simpson Index-------------------------------------------
#Subset Max QID counts
QID_Max <- Sample_Counts %>%
  group_by(SampleID) %>%
  filter(Counts == max(Counts))

#Subset Clone counts by tissue
Sample_Counts_List <-CountsbySample(Raw_Data)

#Calculate Diversity by SampleID
Homogeneity <- GetHomogeneity(Raw_Data)

Simpson_table <- tibble(do.call(rbind, Homogeneity))

##Simpson index is for homogeneity :: S=1 -> homogenetic. Calculate 1/Simpson for relative Diversity final product
Simpson_table <- mutate(Simpson_table, "1/simpson" = 1/Homogeneity)

write.table(Simpson_table,
            "Simpson_table.txt",
            sep = "\t",
            row.names = FALSE)


##for Downsampleing--------
###Find sample with smallest n
tmp <- Raw_Data %>% group_by(SampleID) %>% count()
min(tmp$n)

###Downsample Raw data to SampleID with smallest n
Raw_Data_Down <- Raw_Data %>% group_by(SampleID) %>% sample_n(min(tmp$n), replace = F)
#iNEXT--------------------------------

##Dictate keywords for selecting tissue (Required for DatabySample function)
SubsetKey <- c('B', "L")

##Calculate QID counts by tissue type
QID_Counts_Tissue <- DatabySample(Raw_Down_P1)

#iNEXT requirements: cannot use tibble, all data must be numerical (no QID column)
tibble(rbind(QID_Counts_Tissue[[1]])) -> Blood
Blood <- Blood[,-3]

tibble(rbind(QID_Counts_Tissue[[2]])) -> Lung
Lung <- Lung[,-3]

Blood_Lung <- full_join(x = Blood, y = Lung, by = "QID")
colnames(Blood_Lung) <- c("QID", "Blood", "Lung")

Blood_Lung <- replace_na(Blood_Lung,
                         list(Blood = 0,
                              Lung = 0)
                         )

Blood_Lung$QID <- as.character(Blood_Lung$QID)
Blood_Lung <- data.frame(Blood_Lung[,2:3])

iNEXT_Output1 <- iNEXT(Blood_Lung,
                       q= c(0, 1, 2),
                       datatype = "abundance")

write.table(iNEXT_Output1$AsyEst,
            "iNEXT_Output_Down.txt",
            sep = "\t",
            row.names = FALSE)

ggiNEXT(iNEXT_Output1, type = 1, facet.var = "order") +
  facet_wrap(~order, scales = "free")

##Manual ggINEXT----------
##Subset iNEXT estimates into a dataframe
Blood_iNEXT <- tibble(do.call(rbind, iNEXT_Output1$iNextEst[1]))
Blood_iNEXT[,10] <- "Blood"
colnames(Blood_iNEXT)[10] <- "Tissue"

Lung_iNEXT <- tibble(do.call(rbind, iNEXT_Output1$iNextEst[2]))
Lung_iNEXT[,10] <- "Lung"
colnames(Lung_iNEXT)[10] <- "Tissue"

##Combine Lung and Blood
tmp <- rbind(Blood_iNEXT,
             Lung_iNEXT)

ggplot(tmp[which(tmp$order == 2), ],
       aes(x = m,
           y = qD,
           group = order,
           color = Tissue)
       )+
  geom_jitter()



##Further analysis of iNEXT results-----

###Read data
iNEXT_Est <- read_delim("iNEXT_Output.txt",
                       delim = "\t",
                       escape_double = FALSE,
                       trim_ws = TRUE
                       )

###select diversity tests to subset  
Divtests<- c(unique(iNEXT_Est$Diversity))
###Empty df to record
tmp2 <- tibble()

##Calculate delta of estimates from mean based on each transformed index
for (x in 1:length(Divtests)) {
  tmp <- iNEXT_Est[grep(as.character(Divtests[x]), iNEXT_Est$Diversity),]
  
  if (x == 1) {
    tmp <- tmp %>%
      mutate("trans_est" = Estimator)
    
    tmp[,9] <- mean(tmp$Estimator)
    colnames(tmp)[9] <- "mean_est"
    
    tmp <- tmp %>%
      mutate("delta_est" = trans_est - mean_est)
    
    tmp2 <- rbind(tmp2, tmp)
  }
  else if (x == 2) {
    tmp <- tmp %>%
      mutate("trans_est" = log(Estimator))
    
    tmp[,9] <- log(mean(tmp$Estimator))
    colnames(tmp)[9] <- "mean_est"
    
    tmp <- tmp %>%
      mutate("delta_est" = trans_est - mean_est)
    
    tmp2 <- rbind(tmp2, tmp)
  }
  else if (x == 3){
    tmp <- tmp %>%
      mutate("trans_est" = 1/Estimator)
    
    tmp[,9] <- 1/mean(tmp$Estimator)
    colnames(tmp)[9] <- "mean_est"
    
    tmp <- tmp %>%
      mutate("delta_est" = trans_est - mean_est)
    
    tmp2 <- rbind(tmp2, tmp)
  }
  
}

#Change factors for graphing
tmp2$Diversity <- factor(tmp2$Diversity,
                         levels = c(Divtests))

#Plot delta
ggplot(tmp2, aes(x = Site, y = delta_est))+
  geom_point(size = 4)+
  facet_wrap(facets = "Diversity", scales = "free")+
  geom_hline(yintercept = 0)

ggplot(tmp2[1:2,], aes(x = Site, y = trans_est))+
  geom_point(size = 4)+
  geom_point(aes(x = Site, y = mean_est))

#Plot estimator
ggplot(tmp2, aes(x = Site, y = Estimator, shape = Site))+
  geom_point(size = 4)+
  facet_wrap(facets = "Diversity",
             scales = "free")+
  geom_errorbar(data = filter(tmp2,
                              Diversity == "Species richness"),
                aes(ymin = LCL,
                    ymax = UCL),
                width = .2,
                position = position_dodge(0.7)
                )+
  geom_errorbar(data = filter(tmp2,
                              Diversity == "Shannon diversity"),
                aes(ymin = LCL,
                    ymax = UCL),
                width = .2,
                position = position_dodge(0.7)
                )+
  geom_errorbar(data = filter(tmp2,
                              Diversity == "Simpson diversity"),
                aes(ymin = LCL,
                    ymax = UCL),
                width = .2,
                position = position_dodge(0.7)
                )


#Subset by delta----------------------------------
##Create empty files for recording
QID_Counts_LB_Shared_Close <- tibble()
QID_Counts_LB_Shared_Far <- tibble()

##Subset by clone size
for (x in 1:(nrow(QID_Counts_LB_Shared)/2)) {
  tmp <- rbind(QID_Counts_LB_Shared[(2*x-1),],
               QID_Counts_LB_Shared[2*x,])
  if ((tmp[1,3] - tmp[2,3])^2 < 20) {
    QID_Counts_LB_Shared_Close <- rbind(QID_Counts_LB_Shared_Close,tmp)
  }
  else{
    QID_Counts_LB_Shared_Far <- rbind(QID_Counts_LB_Shared_Far,tmp)
  }
}

##Save data
write.table(QID_Counts_LB_Shared_Close,
            "QID_Counts_LB_Shared_Close.txt",
            sep = "\t",
            row.names = FALSE)

write.table(QID_Counts_LB_Shared_Far,
            "QID_Counts_LB_Shared_Far.txt",
            sep = "\t",
            row.names = FALSE)
