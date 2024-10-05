#Subset counts of each QID by SampleID. Written by Dr. Thomas B. Kepler
CountsbySample <- function(Data, SampleID = SampleID) {
  require(stringr)
  require(dplyr)
  
  #Collect unique sample IDs
  SampleID_Unique <- Data %>%
    distinct({{SampleID}})
  
  #Collect total number of clones in a QID by sample
  Sample_Counts <- Data %>%
    group_by({{SampleID}}) %>%
    count(QID)
  
  #Create empty list to save data
  Sample_Counts_List <- list()
  
  #Separate Sample_Counts by Sample ID
  for (x in 1:dim(SampleID_Unique)[1]) {
    
    df <- Sample_Counts %>%
      filter({{SampleID}} == SampleID_Unique[x,1]) #Subset by SampleID
    
    names(df) <- SampleID_Unique[x,1] #Give each df name
    
    colnames(df) <- c("SampleID", "QID", "Count" ) #apply column headers
    
    Sample_Counts_List[[x]] <- df #Save in list
  }
  
  return(Sample_Counts_List)
}

#Downsample by tissue. Written by Dr. Fumiaki Aihara
Downsample <- function(Data, Group = SampleID, size = 0) {
  
  ###Find sample with smallest n
  tmp <- Data %>%
    group_by({{Group}}) %>%
    count()
  
  #
  if (size > 0) {
    Raw_Data_Down <- Data %>%
      group_by({{Group}}) %>%
      sample_n({{size}},
               replace = F)
  }
  
  else if (size == 0){
    ###Downsample Raw data to SampleID with smallest n
    Raw_Data_Down <- Data %>%
      group_by({{Group}}) %>%
      sample_n(min(tmp$n),
               replace = F)
  }
  
  return(Raw_Data_Down)
}

#Estimate diversity within a sampleID
GetHomogeneity <- function(Data, SampleID = SampleID, QID = QID) {
  require(stringr)
  require(dplyr)
  Homogeneity <- list() #create save file
  
  ##Collect unique sample IDs
  SampleID_Unique <- Data %>%
    distinct({{SampleID}})
  
  ##Collect total number of clones in a QID by sample
  Sample_Counts <- Data %>%
    group_by({{SampleID}}) %>%
    count({{QID}})
  
  ##Empty list to save counts of each QID by SampleID
  Sample_Counts_List <- list()
  
  ##Separate Sample_Counts by Sample ID
  for (x in 1:dim(SampleID_Unique)[1]) {
    
    df <- Sample_Counts %>%
      filter({{SampleID}} == SampleID_Unique[x,1]) #Subset by SampleID
    
    names(df) <- SampleID_Unique[x,1] #Give df name
    
    colnames(df) <- c("SampleID", "QID", "Count" ) #apply column headers
    
    Sample_Counts_List[[x]] <- df #Save in list
  }
  
  ##Estimate diversity
  ###Set progress bar
  pb <- txtProgressBar(min = 0,
                       max = nrow(SampleID_Unique),
                       style = 3
                       )
  
  for (y in 1:nrow(SampleID_Unique)) {
    
    tmp <- tibble(Sample_Counts_List[[y]]) #Subset given SampleID table from main list
    
    tmp2 <- tibble() #Create save file
    
    ###Calculation
    for (x in 1:nrow(tmp)) {
      n <- sum(tmp$Count)
      #Calculate (xi * (xi-1))/(n*(n-1)) for each QID where xi:: # of seq in a QID, n:: total # of Seq
      tmp2[x,1] <- ((tmp[x,3])*(tmp[x,3]-1))/(n*(n-1))
      
      #Sum results and save in Diversity list
      Homogeneity[[y]] <- cbind(tmp[1,1],
                              "Homogeneity" = sum(tmp2)
                              )
    }
    #update progress bar
    setTxtProgressBar(pb, y)
  }
  return(Homogeneity)
  close(pb)
}

#Subset Raw data by sample
DatabySample<- function(Data, QID = QID, SubsetKey = c("B", "L")) {
  
  require(stringr)
  require(dplyr)
  
  QID_Counts_Tissue <- list()
  
  #Separate Sample_Counts by Sample ID
  for (x in 1:length(SubsetKey)) {
    tmp <- filter(Data,
                    grepl(as.character(SubsetKey[x]),
                          SampleID)
                  )
    
    #Collect total number of clones in a QID
    df <- tmp%>%
      group_by({{QID}}) %>%
      dplyr::count({{QID}}) #%>%
      #rename(Counts = n)
    
    tmp2 <- semi_join(df, tmp, by = 'QID')
    
    colnames(df) <- c("QID", "Count") #apply column headers
    
    df$Sample <- as.character(SubsetKey[x])
    
    QID_Counts_Tissue[[x]] <- df #Save in list
    names(QID_Counts_Tissue)[[x]] <- as.character(SubsetKey[x])
  }
  return(QID_Counts_Tissue)
}

#Top 10 QID
Top_QID <- function(Raw_Data) {
  
  require(stringr)
  require(dplyr)
  
  BloodorLung <- c('B', "L") #Keywords to dictate tissue type, change as needed
  
  QID_Counts_Tissue <- DatabySample(Raw_Data) #Subset Raw data by sample
  
  QID_Counts_LB <- tibble(do.call(rbind,
                                  QID_Counts_Tissue
                                  )
                          ) #Convert subsetted data to tibble
  
  #Find unique QID
  ##Count Degeneracy
  QID_Counts_LB_Degeneracy <- QID_Counts_LB %>%
    group_by(Sample) %>%
    distinct(QID) %>%
    ungroup() %>%
    count(QID) %>%
    rename(Degeneracy = n) %>%
    left_join(QID_Counts_LB,
              by = "QID"
              )
  
  ##Subset tissue unique QID
  QID_Counts_LB_Unique <- QID_Counts_LB_Degeneracy[which(QID_Counts_LB_Degeneracy$Degeneracy == 1),
                                                   ]
  
  ##Subset top largest unique QID by sample
  Top_LB_Unique <- QID_Counts_LB_Unique %>%
    group_by(Sample) %>%
    top_n(10, Count) #Change threshold needed here
  
  ##Calculate proportion of each top cluster by all counts
  Top_LB_Unique <- Top_LB_Unique %>%
    mutate(Percent_Total = Count/dim(Raw_Data)[1])
  
  ##Create save file for for loop
  Save_List <- list()
  
  ##Calculate proportion of each top cluster by all counts within a tissue
  for (x in 1:2) {
    ##Subset raw data by sample
    Raw_Sample <- filter(Raw_Data,grepl(as.character(BloodorLung[x]),
                                        SampleID
                                        )
                         )
    
    ##Calculate proportion of top QID by sample
    tmp <- filter(Top_LB_Unique,
                  Sample == as.character(BloodorLung[x])
                  ) #Subset by sample
    
    Proportion <- tmp %>%
      mutate(Percent_Sample = Count/dim(Raw_Sample)[1]) #Calculate proportion
    
    Save_List[[x]] <- Proportion #Save in list
  }
  Top_LB_Unique <- tibble(do.call(rbind,
                                  Save_List)
                          ) #Convert list to tibble
  
  ##Rename B to Blood
  Top_LB_Unique$Sample <- str_replace(Top_LB_Unique$Sample,
                                      pattern = "B",
                                      replacement = "Blood")
  
  ##Rename L to Lung
  Top_LB_Unique$Sample <- str_replace(Top_LB_Unique$Sample,
                                      pattern = "L",
                                      replacement = "Lung")
  
  return(Top_LB_Unique)
}

#Jaccard index
jaccard_function <- function(Data) {
  require(dplyr)
  
  #Data preparation--------------------------------
  ##Calculate QID counts by tissue type
  QID_Counts_Tissue <- DatabySample(Data)
  
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
  
  ##Subset shared QIDs by Degeneracy
  QID_Counts_LB_Shared <- QID_Counts_LB_Degeneracy[which(QID_Counts_LB_Degeneracy$Degeneracy == 2), ]
  
  if (nrow(QID_Counts_LB_Shared)==0) {
    stop("No shared sequences")
  } 
  
  else {
    #Calculate Jaccard Index by CDR3----------------------------------
    ##Identify Unique QIDs in Shared
    QID_Unique_Shared <- QID_Counts_LB_Shared %>%
      distinct(QID)
    
    ##Subset shared QIDs from Raw_Data to Raw_Shared
    Raw_Shared <- tibble()
    
    for (x in 1:nrow(QID_Unique_Shared)) {
      tmp <- Data[which(Data$QID == as.numeric(QID_Unique_Shared[x,1]
      )
      ),
      ]
      
      Raw_Shared <- rbind(Raw_Shared,tmp)
    }
    
    ##Create save files
    tmp2 <- list()
    Output <- QID_Unique_Shared
    
    ##Calculate Jaccard Index
    pb <- txtProgressBar(min = 0,
                         max = nrow(QID_Unique_Shared),
                         style = 3
                         )
    
    for (x in 1:nrow(QID_Unique_Shared)) {
      #Subset by QID
      tmp <- filter(Raw_Shared,
                    QID == as.numeric(QID_Unique_Shared[x,1]
                    )
      )
      #Subset by tissue type
      setA <- filter(tmp,
                     grepl(("B"),
                           SampleID
                     )
      )
      
      setB <- filter(tmp,
                     grepl(("L"),
                           SampleID
                     )
      )
      
      #Calculate intersect and union
      I <- intersect(setA$CDR3, setB$CDR3)
      
      U <- union(setA$CDR3, setB$CDR3)
      
      #Calculate Jaccard
      tmp2 <-  length(I)/length(U)*100
      
      #Save info
      Output[x, 2] <- tmp2
      Output[x, 3] <- length(I)
      Output[x, 4] <- length(U)
      Output[x, 5] <- nrow(setA)
      Output[x, 6] <- nrow(setB)
      
      #update progress bar
      setTxtProgressBar(pb, x)
    }
    
    colnames(Output) <- c("QID",
                          "Jaccard",
                          "Intersect",
                          "Union",
                          "Seq_Count_Blood",
                          "Seq_Count_Lung"
    )
    
    Output <- mutate(Output,
                     Total_Count = Seq_Count_Blood + Seq_Count_Lung
    )
    return(Output)
    close(pb)
  }
}

#iNEXT down sample data process for lapply. Written by Dr. Fumiaki Aihara
Process_Div <- function(Data, target_down = 1000, SubsetKey = c('B', "L"), Patient_Numb = c(1:2), QID = QID) {
  #Downsample by smallest sample (size = 0 will downsample by smallest)----
  print("Downsampling")
  Raw_Down <- list()
  for (x in 1:2) {
    Raw_Down[[x]] <- Downsample(Raw_Data[[x]], size = target_down)
    names(Raw_Down)[[x]] <- as.character(paste0("P", x,"_Down"))
  }
  
  #iNEXT--------------------------------
  print("Preparing Data")
  ##Save file for final output
  iNEXT_Output <- list()
  
  ##Calculate QID counts by tissue type
  QID_Counts_Tissue <- list()
  
  for (x in 1:length(Raw_Down)) {
    QID_Counts_Tissue[[x]] <- DatabySample(Raw_Down[[x]])
  }
  
  #iNEXT requirements: cannot use tibble, all data must be numerical (no QID column)
  ##Subset QID counts data by tissue and patient number
  for (x in 1:length(Patient_Numb)) {
    assign(paste0("DF_B_P", Patient_Numb[x]), tibble(rbind(QID_Counts_Tissue[[x]][[1]][,-3])))
  }
  for (x in 1:length(Patient_Numb)) {
    assign(paste0("DF_L_P", Patient_Numb[x]), tibble(rbind(QID_Counts_Tissue[[x]][[2]][,-3])))
  }
  
  ##Join data within patient number
  QID_Count_Join <- list()
  for (x in 1:length({{SubsetKey}})) {
    for (y in 1:length(Patient_Numb)) {
      tmp <- mget(paste0("DF_", SubsetKey[x], "_P", Patient_Numb[y]))
      QID_Count_Join <- c(QID_Count_Join,tmp)
    }
  }
  ##Join data by patient number and reduce to one dataframe
  tmp <- reduce(QID_Count_Join, full_join, by = "QID")
  colnames(tmp)[2:ncol(tmp)] <- names(QID_Count_Join)
  
  #Replace NAs with 0's
  tmp[, 2:ncol(tmp)][is.na(tmp[, 2:ncol(tmp)])] <- 0
  
  #Remove QID column
  tmp2 <- as.data.frame(tmp[,-1])
  
  #Run iNEXT within patients
  results <- list()
  pb <- txtProgressBar(min = 0, max = length(Patient_Numb), style = 3)
  
  #iNEXT
  for (x in 1:length(Patient_Numb)) {
    print(paste0("Starting iNEXT P", x))
    setTxtProgressBar(pb, x)
    tmp3 <- select(tmp2, matches(paste0("P", x)))
    results[[length(results)+1]] <- list(iNEXT(tmp3,
                                               q = c(0, 1, 2),
                                               datatype = "abundance"))
    Sys.sleep(0.3)
  }
  close(pb)
  
  Est_list <- list()
  for (x in 1:length(results)) {
    Est_list[[x]] <- results[[x]][[1]]$AsyEst
  }
  return(Est_list)
}

#Hill number non-iNEXT variant. Written by Dr. Thomas B. Kepler
HillNumbers <- function(x)
{
  n <- length(x)
  total <- sum(x)
  p <- x/total
  h0 <- 0
  h1 <- 0
  h2 <- 0
  for (q in p)
  {
    if (q > 0)
    {
      h0 <- h0 + 1
      h1 <- h1 - q * log(q)
      h2 <- h2 + q^2
    }
  }
  return(c(h0,exp(h1),1/h2))
}

#Extract estimated values from iNEXT AsyEst result. Written by Dr. Thomas B. Kepler
extractEst <- function(x){
  AsyEst1 <- data.frame(x$AsyEst)
  colnames(AsyEst1) <- c(colnames(x$AsyEst))
  AsyEst1 <- tibble::rownames_to_column(AsyEst1, "Diversity")
  
  #wide format
  AsyEst_Wide <- AsyEst1 %>%
    select(Diversity, Estimator) %>%
    pivot_wider(#AsyEst1,
                #id_cols = Diversity,
                names_from = Diversity,
                values_from = Estimator)
  return(AsyEst_Wide)
}

# Divergence Ratios. Mathmatical formulas developed and written by Dr. Thomas B. Kepler. Implimentation and arguments written by Dr. Fumiaki Aihara
HillAnalysis <- function(x1, x2, Call_iNEXT=FALSE)
{
  n <- length(x1)
  if (length(x2) != n)
  {
    return("Lengths must be the same")
  }
  
  if (Call_iNEXT == TRUE){
    #Run iNEXT
    ##print(paste0("iNEXT sample ", x1))
    NEXT1 <- iNEXT(x1, q = c(0, 1, 2), datatype = "abundance")
    
    ##print(paste0("iNEXT sample ", x2))
    NEXT2 <- iNEXT(x2, q = c(0, 1, 2), datatype = "abundance")
    
    #Extract Estimates
    H1 <- extractEst(NEXT1)
    H2 <- extractEst(NEXT2)
    
    #Calculations
    phi1 <- sum(x1)
    phi2 <- sum(x2)
    tot <- phi1 + phi2
    phi1 <- phi1/tot
    phi2 <- phi2/tot
    
    y <- x1 + x2
    HP <- HillNumbers(y)
    rho0 <- 2 * HP[1]/(H1[1] + H2[1])
    rho1 <- HP[2] / ( H1[2]^phi1 * H2[2]^phi2)
    rho2 <- ( phi1^2 / H1[3] + phi2^2 / H2[3] ) * HP[3]/(phi1^2 + phi2^2)
    
    return (c(rho0,rho1,rho2))
  }
  else if (Call_iNEXT == FALSE){
  H1 <- HillNumbers(x1)
  H2 <- HillNumbers(x2)
  phi1 <- sum(x1)
  phi2 <- sum(x2)
  tot <- phi1 + phi2
  phi1 <- phi1/tot
  phi2 <- phi2/tot
  
  y <- x1 + x2
  HP <- HillNumbers(y)
  rho0 <- 2 * HP[1]/(H1[1] + H2[1])
  rho1 <- HP[2] / ( H1[2]^phi1 * H2[2]^phi2)
  rho2 <- ( phi1^2 / H1[3] + phi2^2 / H2[3] ) * HP[3]/(phi1^2 + phi2^2)
  
  return (c(rho0,rho1,rho2))
  }
}
