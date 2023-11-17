library(ggplot2)
library(tidyverse)
library(dplyr)
library(xlsx)
library(reshape)
library(corrplot)
library(janitor)
palette = colorRampPalette(c("blue", "white", "red")) (20)
#function to round off all digits in a data frame
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}
# Read tsv files using read.delim() method
df<-read.delim( "/Users/ramsey.najm/Desktop/GitHub/Projects/RN_CON/SC6-19_Metadata/meta.data copy.tsv",sep="\t")
is.data.frame(df)
#order the rows by ascending order of network size. In the  end unnecessary, but easier to red
dforder <- df%>%arrange(UVI.assignment.f4)
head(dforder)#Identify the unique genotypes in the dataset
dforder$Genotype
unique(dforder$Days.Category)
x <- unique(dforder$Genotype)
x
#generate individual data frames for the cell types containing only the tracing condition at day 150 and 
#remove NA values
#######################################################
dfH9 <- filter(dforder,  Genotype == 'H9', Condition == 'Tracing', Days.Category %in% c('150','D150'))
dfH9 <- dfH9 %>% drop_na(Network.Size.F4)
#######################################################
dfH9120 <- filter(dforder,  Genotype == 'H9', Condition == 'Tracing', Days.Category %in% c('120'))
dfH9120 <- dfH9120 %>% drop_na(Network.Size.F4)
#######################################################
df177<- filter(dforder,  Genotype == '177', Condition == 'Tracing', Days.Category %in% c('150','D150'))
df177 <- df177 %>% drop_na(Network.Size.F4)
#######################################################
df176<- filter(dforder,  Genotype == '176', Condition == 'Tracing', Days.Category %in% c('150','D150'))
df176 <- df176 %>% drop_na(Network.Size.F4)
#Subset disease conditions by their genotype, disese context (mutant (M) v. control (C), 
#and tracing and NA values)
dfp24_M <- filter(dforder, Genotype == 'TSC.p24', Condition == 'Tracing', Days.Category %in% c('150','D150'), grepl('TSC.p24.M',Tracing.Cell.Lines))
dfp24_M <- dfp24_M %>% drop_na(Network.Size.F4)
########################################################
dfp24_C<- filter(dforder, Genotype == 'TSC.p24', Condition == 'Tracing', Days.Category %in% c('150','D150'),
                 grepl('TSC.p24.C',Tracing.Cell.Lines))
dfp24_C <- dfp24_C %>% drop_na(Network.Size.F4)
########################################################
dfp5_M<- filter(dforder, Genotype == 'TSC.p5', Condition == 'Tracing', Days.Category %in% c('150','D150'),
                grepl('TSC.p5.M',Tracing.Cell.Lines))
dfp5_M <- dfp5_M %>% drop_na(Network.Size.F4)
########################################################
dfp5_C <- filter(dforder, Genotype == 'TSC.p5', Condition == 'Tracing', Days.Category %in% c('150','D150'),
                 grepl('TSC.p5.C',Tracing.Cell.Lines))
dfp5_C <- dfp5_C %>% drop_na(Network.Size.F4)
########################################################
dfARID_M<- filter(dforder, Genotype == 'ARID.p2', Condition == 'Tracing', Days.Category %in% c('150','D150'),
                  grepl('ARID1B.p2.M',Tracing.Cell.Lines))
dfARID_M <- dfARID_M %>% drop_na(Network.Size.F4)
########################################################
dfARID_C <- filter(dforder, Genotype == 'ARID.p2', Condition == 'Tracing', Days.Category %in% c('150','D150'),
                   grepl('ARID1B.p2.C',Tracing.Cell.Lines))
dfARID_C <- dfARID_C %>% drop_na(Network.Size.F4)
########################################################################################################################################################################
#Now we hve individual data sets I want to generate the percentage that all network sizes make up in those respective datasets
df.list <- list( 'dfH9_netfreq' = dfH9, 'dfH9120_netfreq' = dfH9120,'df177_netfreq' = df177, 'df176_netfreq' = df176, 'dfp24_M_netfreq' = dfp24_M,
                 'dfp24_C_netfreq' = dfp24_C,'dfp5_M_netfreq' = dfp5_M,'dfp5_C_netfreq' = dfp5_C,'dfARID_M_netfreq' = dfARID_M,
                 'dfARID_C_netfreq' = dfARID_C)

df.list.filt <- lapply(df.list, function(x){
  y <- x %>%
    group_by(Network.Size.F4) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  y <- as.data.frame(y)
  y$nNetwork <- y$n / y$Network.Size.F4
  #y <- filter(y, Network.Size.F4 >'1')
  y$freqperc <- 100*y$freq 
  y <- round_df(y, digits = 2)
  
})
#Make the names of df.list.filt = names of df.list
names(df.list.filt)<-names(df.list)
#assign a variable to the names of df.list.filt
data_names <- names(df.list.filt)
#pull out all individual data frames from the final list and save to the global environment
list2env(df.list.filt,envir=.GlobalEnv)
#write CSVs to a set directory
out <- "/Users/ramsey.najm/Desktop/GitHub/Projects/RN_CON/SC6-19_Metadata/Plots/NetworkFreq_All/"
for(i in 1:length(data_names)) {                              # Head of for-loop
  print(i)
  write.csv(get(data_names[i]), 
            file =  paste0(out, data_names[i],"_all.csv"), 
            row.names = FALSE) # Write CSV files to folder in directory specified by out
}
#######################################################################################################################################################################################################################################
#Network composition per cluster for all networks
df.list <- list( 'dfH9_netnumb' = dfH9, 'dfH9120_netnumb' = dfH9120,'df177_netnumb' = df177, 'df176_netnumb' = df176, 'dfp24_M_netnumb' = dfp24_M,
                 'dfp24_C_netnumb' = dfp24_C,'dfp5_M_netnumb' = dfp5_M,'dfp5_C_netnumb' = dfp5_C,'dfARID_M_netnumb' = dfARID_M,
                 'dfARID_C_netnumb' = dfARID_C)
#Generate function and apply it to all data frames in that list
df.list.filt <- lapply(df.list, function(x){
  tab <- table(x$UVI.assignment.f4, x$RNA_snn_res.0.1.ordered)
  y <- as.data.frame.matrix(tab)
  #filter for networks greater than 1
  #df <- y %>% filter(rowSums(.)>1)
})
#Make the names of df.list.filt = names of df.list
names(df.list.filt)<-names(df.list)
#assign a variable to the names of df.list.filt
data_names <- names(df.list.filt)
#pull out all individual data frames from the final list and save to the global environment
list2env(df.list.filt,envir=.GlobalEnv)
# #Add the row variance between columns to every row
# df.list <- list( 'dfH9_netnumb' = dfH9_netnumb, 'dfH9120_netnumb' = dfH9120_netnumb,'df177_netnumb' = df177_netnumb, 'df176_netnumb' = df176_netnumb, 'dfp24_M_netnumb' = dfp24_M_netnumb,
#                  'dfp24_C_netnumb' = dfp24_C_netnumb,'dfp5_M_netnumb' = dfp5_M_netnumb,'dfp5_C_netnumb' = dfp5_C_netnumb,'dfARID_M_netnumb' = dfARID_M_netnumb,
#                  'dfARID_C_netnumb' = dfARID_C_netnumb)
# 
# 
# df.list.filt <- lapply(df.list, function(x){
#   x$RowVariance <- apply(x, 1, var)
#   x
# })
# names(df.list.filt)<-names(df.list)
# #assign a variable to the names of df.list.filt
# data_names <- names(df.list.filt)
# #pull out all individual data frames from the final list and save to the global environment
# list2env(df.list.filt,envir=.GlobalEnv)
#write CSVs to a set directory
out <- "/Users/ramsey.najm/Desktop/GitHub/Projects/RN_CON/SC6-19_Metadata/Plots/NetworkNumb_All/"
for(i in 1:length(data_names)) {                              # Head of for-loop
  print(i)
  write.csv(get(data_names[i]), 
            file =  paste0(out, data_names[i],"_all.csv"), 
            row.names = FALSE) # Write CSV files to folder in directory specified by out
}
#######################################################################################################################################################################################################################################
#Heatmaps of network composition greater than 1 across cell lines
df.list <- list( 'dfH9_netnumb_mat' = dfH9, 'dfH9120_netnumb_mat' = dfH9120,'df177_netnumb_mat' = df177, 'df176_netnumb_mat' = df176, 'dfp24_M_netnumb_mat' = dfp24_M,
                 'dfp24_C_netnumb_mat' = dfp24_C,'dfp5_M_netnumb_mat' = dfp5_M,'dfp5_C_netnumb_mat' = dfp5_C,'dfARID_M_netnumb_mat' = dfARID_M,
                 'dfARID_C_netnumb_mat' = dfARID_C)
#Generate function and apply it to all data frames in that list
df.list.filt <- lapply(df.list, function(x){
  tab <- table(x$UVI.assignment.f4, x$RNA_snn_res.0.1.ordered)
  y <- as.data.frame.matrix(tab)
  #filter for networks greater than 1
  df <- y %>% filter(rowSums(.)>1)
  df <- as.matrix(df)
})
#Heatmap of network composition for every scSeq
names(df.list.filt)<-names(df.list)
data_names <- names(df.list.filt)
#pull out all individual data frames from the final list and save to the global environment
list2env(df.list.filt,envir=.GlobalEnv)
#Generate outgoing directory
out <- "/Users/ramsey.najm/Desktop/GitHub/Projects/RN_CON/SC6-19_Metadata/Plots/NetworkCompHM_over1/"
# make and save heat maps from the list of matrices generated in the above function
mapply(function(x, y){
  png(file=paste0(out, y, "_NetHM_over1.png"))
  heatmap(x, Rowv = NULL, Colv= NA, 
          xlab = 'Seurat Cluster' , ylab = 'Network',col = palette, main= paste0(y,'_NetHM > 1'))
  dev.off()
}, x = df.list.filt, y = data_names)

#######################################################################################################################################################################################################################################
#Heatmaps of network composition greater than 3 across cell lines
df.list <- list( 'dfH9_netnumb_mat' = dfH9, 'dfH9120_netnumb_mat' = dfH9120,'df177_netnumb_mat' = df177, 'df176_netnumb_mat' = df176, 'dfp24_M_netnumb_mat' = dfp24_M,
                 'dfp24_C_netnumb_mat' = dfp24_C,'dfp5_M_netnumb_mat' = dfp5_M,'dfp5_C_netnumb_mat' = dfp5_C,'dfARID_M_netnumb_mat' = dfARID_M,
                 'dfARID_C_netnumb_mat' = dfARID_C)
#Generate function and apply it to all data frames in that list
df.list.filt <- lapply(df.list, function(x){
  tab <- table(x$UVI.assignment.f4, x$RNA_snn_res.0.1.ordered)
  y <- as.data.frame.matrix(tab)
  #filter for networks greater than 3
  df <- y %>% filter(rowSums(.)>3)
  df <- as.matrix(df)
})
#Heatmap of network composition for every scSeq
names(df.list.filt)<-names(df.list)
data_names <- names(df.list.filt)
#pull out all individual data frames from the final list and save to the global environment
list2env(df.list.filt,envir=.GlobalEnv)
#Generate outgoing directory
out <- "/Users/ramsey.najm/Desktop/GitHub/Projects/RN_CON/SC6-19_Metadata/Plots/NetworkCompHM_over3/"
# make and save heat maps from the list of matrices generated in the above function
mapply(function(x, y){
  png(file=paste0(out, y, "_NetHM_over3.png"))
  heatmap(x, Rowv = NULL, Colv= NA, 
          xlab = 'Seurat Cluster' , ylab = 'Network',col = palette, main= paste0(y,'_NetHM > 3'))
  dev.off()
}, x = df.list.filt, y = data_names)

#######################################################################################################################################################################################################################################
#Network celltype composition by percentage
df.list <- list( 'dfH9_netperc' = dfH9, 'dfH9120_netperc' = dfH9120, 'df177_netperc' = df177, 'df176_netperc' = df176, 'dfp24_M_netperc' = dfp24_M,
                 'dfp24_C_netperc' = dfp24_C,'dfp5_M_netperc' = dfp5_M,'dfp5_C_netperc' = dfp5_C,'dfARID_M_netperc' = dfARID_M,
                 'dfARID_C_netperc' = dfARID_C)
#Generate function and apply it to all data frames in that list
df.list.filt <- lapply(df.list, function(x){
  tab <- table(x$UVI.assignment.f4, x$RNA_snn_res.0.1.ordered)
  y <- 100* tab/rowSums(tab)
  z <- as.data.frame.matrix(y)
})
#Make the names of the df that comes out of the function the same as the original list
names(df.list.filt)<-names(df.list)
#pull out all individual data frames from the final list and save to the global environment
list2env(df.list.filt,envir=.GlobalEnv)
#######################################################################################################################################################################################################################################
#Next let's do this for networks >3, what we did above was all networks so there are networks of size 1,2, or 3 
#that are not necessarily the most interesting

#Filter object for networks greater than 3
df.list <- list( 'dfH9_3' = dfH9, 'dfH9120_3' = dfH9120, 'df177_3' = df177, 
                 'df176_3' = df176, 'dfp24_M_3' = dfp24_M,
                 'dfp24_C_3' = dfp24_C,'dfp5_M_3' = dfp5_M,
                 'dfp5_C_3' = dfp5_C,'dfARID_M_3' = dfARID_M,
                 'dfARID_C_3' = dfARID_C)
#Generate function and apply it to all data frames in that list
df.list.filt <- lapply(df.list, function(x){
  x %>% filter( Network.Size.F4 > 3)
})
names(df.list.filt)<-names(df.list)
#pull out all individual data frames from the final list and save to the global environment
list2env(df.list.filt,envir=.GlobalEnv)
#######################################################################################################################################################################################################################################
#Make new object of network percentages but including only networks greater than 3
#######################################################################################################################################################################################################################################
df.list <- list( 'dfH9_netperc_3' = dfH9_3, 'dfH9120_netperc_3' = dfH9120_3, 'df177_netperc_3' = df177_3, 
                 'df176_netperc_3' = df176_3, 'dfp24_M_netperc_3' = dfp24_M_3,
                 'dfp24_C_netperc_3' = dfp24_C_3,'dfp5_M_netperc_3' = dfp5_M_3,'dfp5_C_netperc_3' = dfp5_C_3,
                 'dfARID_M_netperc_3' = dfARID_M_3,'dfARID_C_netperc_3' = dfARID_C_3)
#Generate function and apply it to all data frames in that list
df.list.filt <- lapply(df.list, function(x){
  tab <- table(x$UVI.assignment.f4, x$RNA_snn_res.0.1.ordered)
  y <- 100* tab/rowSums(tab)
  z <- as.data.frame.matrix(y)
})
#Make the names of the df that comes out of the function the same as the original list
names(df.list.filt)<-names(df.list)
#pull out all individual data frames from the final list and save to the global environment
list2env(df.list.filt,envir=.GlobalEnv)
#######################################################################################################################################################################################################################################
#Transposing network percentages so that we can make correlation heatmaps
df.list <- list( 'dfH9_netperc_3T' = dfH9_netperc_3, 'dfH9120_netperc_3T' = dfH9120_netperc_3, 
                 'df177_netperc_3T' = df177_netperc_3, 'df176_netperc_3T' = df176_netperc_3, 'dfp24_M_netperc_3T' = dfp24_M_netperc_3,
                 'dfp24_C_netperc_3T' = dfp24_C_netperc_3,'dfp5_M_netperc_3T' = dfp5_M_netperc_3,'dfp5_C_netperc_3T' = dfp5_C_netperc_3,
                 'dfARID_M_netperc_3T' = dfARID_M_netperc_3,'dfARID_C_netperc_3T' = dfARID_C_netperc_3)

df.list.filt <- lapply(df.list, function(x){
  tx <- t(x)
  as.data.frame(tx)
  cor(tx)
})
#match the names of the two lists (before and after performing the function)
names(df.list.filt)<-names(df.list)
#pull out all individual data frames from the final list and save to the global environment
list2env(df.list.filt,envir=.GlobalEnv)
#Heatmap of network composition for every scSeq
data_names <- names(df.list.filt)
#Generate outgoing directory
out <- "/Users/ramsey.najm/Desktop/GitHub/Projects/RN_CON/SC6-19_Metadata/Plots/NetComp_Correlation_over3/"
# make and save heat maps from the list of matrices generated in the above function
mapply(function(x, y){
  png(file=paste0(out, y, "_NetComp_over3.png"))
  heatmap(x, col = palette, symm = TRUE)
  dev.off()
}, x = df.list.filt, y = data_names)
#######################################################################################################################################################################################################################################
#######################################################################################################################################################################################################################################
#######################################################################################################################################################################################################################################
#Next let's do this for networks >1
#Filter object for networks greater than 1
df.list <- list( 'dfH9_1' = dfH9, 'dfH9120_1' = dfH9120, 'df177_1' = df177, 
                 'df176_1' = df176, 'dfp24_M_1' = dfp24_M,
                 'dfp24_C_1' = dfp24_C,'dfp5_M_1' = dfp5_M,
                 'dfp5_C_1' = dfp5_C,'dfARID_M_1' = dfARID_M,
                 'dfARID_C_1' = dfARID_C)
#Generate function and apply it to all data frames in that list
df.list.filt <- lapply(df.list, function(x){
  x %>% filter( Network.Size.F4 > 1)
})
names(df.list.filt)<-names(df.list)
#pull out all individual data frames from the final list and save to the global environment
list2env(df.list.filt,envir=.GlobalEnv)
#######################################################################################################################################################################################################################################
#Make new object of network percentages but including only networks greater than 3
#######################################################################################################################################################################################################################################
df.list <- list( 'dfH9_netperc_1' = dfH9_1, 'dfH9120_netperc_1' = dfH9120_1, 'df177_netperc_1' = df177_1, 
                 'df176_netperc_1' = df176_1, 'dfp24_M_netperc_1' = dfp24_M_1,
                 'dfp24_C_netperc_1' = dfp24_C_1,'dfp5_M_netperc_1' = dfp5_M_1,'dfp5_C_netperc_1' = dfp5_C_1,
                 'dfARID_M_netperc_1' = dfARID_M_1,'dfARID_C_netperc_1' = dfARID_C_1)
#Generate function and apply it to all data frames in that list
df.list.filt <- lapply(df.list, function(x){
  tab <- table(x$UVI.assignment.f4, x$RNA_snn_res.0.1.ordered)
  y <- 100* tab/rowSums(tab)
  z <- as.data.frame.matrix(y)
})
#Make the names of the df that comes out of the function the same as the original list
names(df.list.filt)<-names(df.list)
#pull out all individual data frames from the final list and save to the global environment
list2env(df.list.filt,envir=.GlobalEnv)
#######################################################################################################################################################################################################################################
#Transposing network percentages so that we can make correlation heatmaps
df.list <- list( 'dfH9_netperc_1T' = dfH9_netperc_1, 'dfH9120_netperc_1T' = dfH9120_netperc_1, 
                 'df177_netperc_1T' = df177_netperc_1, 'df176_netperc_1T' = df176_netperc_1, 'dfp24_M_netperc_1T' = dfp24_M_netperc_1,
                 'dfp24_C_netperc_1T' = dfp24_C_netperc_1,'dfp5_M_netperc_1T' = dfp5_M_netperc_1,'dfp5_C_netperc_1T' = dfp5_C_netperc_1,
                 'dfARID_M_netperc_1T' = dfARID_M_netperc_1,'dfARID_C_netperc_1T' = dfARID_C_netperc_1)

df.list.filt <- lapply(df.list, function(x){
  tx <- t(x)
  as.data.frame(tx)
  cor(tx)
})
#match the names of the two lists (before and after performing the function)
names(df.list.filt)<-names(df.list)
#pull out all individual data frames from the final list and save to the global environment
list2env(df.list.filt,envir=.GlobalEnv)
#Heatmap of network composition for every scSeq
data_names <- names(df.list.filt)
#Generate outgoing directory
out <- "/Users/ramsey.najm/Desktop/GitHub/Projects/RN_CON/SC6-19_Metadata/Plots/NetComp_Correlation_over1/"
# make and save heat maps from the list of matrices generated in the above function
mapply(function(x, y){
  png(file=paste0(out, y, "_NetComp_over1.png"))
  heatmap(x, col = palette, symm = TRUE)
  dev.off()
}, x = df.list.filt, y = data_names)
#######################################################################################################################################################################################################################################
################################################################################################################################################################################################
#Next, divide up each network size (i.e., network size 4) into individual scSeq experiments to make
#datapoints in the network size histrograms so we can do some statistics instead of just the average


dfH9_split <- dfH9 %>% group_split(dfH9$library)
dfH9120_split <- dfH9120 %>% group_split(dfH9120$library)
df177_split <- df177 %>% group_split(df177$library)
df176_split <- df176 %>% group_split(df176$library)
dfp24_M_split <- dfp24_M %>% group_split(dfp24_M$library)
dfp5_M_split <- dfp5_M %>% group_split(dfp5_M$library)
dfp24_C_split <- dfp24_C %>% group_split(dfp24_C$library)
dfp5_C_split <- dfp5_C %>% group_split(dfp5_C$library)
dfARID_M_split <- dfARID_M %>% group_split(dfARID_M$library)
dfARID_C_split <- dfARID_C %>% group_split(dfARID_C$library)
################################################################
dfH9_sc8.325 <- dfH9_split[[1]]
dfH9_sc8.325 <- as.data.frame(dfH9_sc8.325)

dfH9_sc10.949<- dfH9_split[[2]]
dfH9_sc10.949 <- as.data.frame(dfH9_sc10.949)

dfH9_sc11.179<- dfH9_split[[3]]
dfH9_sc11.179<- as.data.frame(dfH9_sc11.179)

dfH9_sc11.180<- dfH9_split[[4]]
dfH9_sc11.180 <- as.data.frame(dfH9_sc11.180)

dfH9_sc11.181<- dfH9_split[[5]]
dfH9_sc11.181 <- as.data.frame(dfH9_sc11.181)

dfH9_sc11.182<- dfH9_split[[6]]
dfH9_sc11.182 <- as.data.frame(dfH9_sc11.182)

dfH9_sc11.421<- dfH9_split[[7]]
dfH9_sc11.421 <- as.data.frame(dfH9_sc11.421)

dfH9_sc12.380<- dfH9_split[[8]]
dfH9_sc12.380 <- as.data.frame(dfH9_sc12.380)

dfH9_sc12.381<- dfH9_split[[9]]
dfH9_sc12.381 <- as.data.frame(dfH9_sc12.381)

dfH9_sc12.382<- dfH9_split[[10]]
dfH9_sc12.382 <- as.data.frame(dfH9_sc12.382)

dfH9_sc13.767<- dfH9_split[[11]]
dfH9_sc13.767 <- as.data.frame(dfH9_sc13.767)

dfH9_sc13.768<- dfH9_split[[12]]
dfH9_sc13.768 <- as.data.frame(dfH9_sc13.768)
###################################################
dfH9120_sc06.692 <- dfH9120_split[[1]]
dfH9120_sc06.692 <- as.data.frame(dfH9120_sc06.692)

dfH9120_sc09.239<- dfH9120_split[[2]]
dfH9120_sc09.239<- as.data.frame(dfH9120_sc09.239)
###################################################
df177_sc12.380 <- df177_split[[1]]
df177_sc12.380 <- as.data.frame(df177_sc12.380)

df177_sc12.381 <- df177_split[[2]]
df177_sc12.381 <- as.data.frame(df177_sc12.381)

df177_sc12.382 <- df177_split[[3]]
df177_sc12.382 <- as.data.frame(df177_sc12.382)

df177_sc15.033 <- df177_split[[4]]
df177_sc15.033 <- as.data.frame(df177_sc15.033)

df177_sc15.035 <- df177_split[[5]]
df177_sc15.035 <- as.data.frame(df177_sc15.035)

df177_sc17.041 <- df177_split[[6]]
df177_sc17.041 <- as.data.frame(df177_sc17.041)

df177_sc17.042 <- df177_split[[7]]
df177_sc17.042 <- as.data.frame(df177_sc17.042)
###################################################
df176_sc14.002 <- df176_split[[1]]
df176_sc14.002 <- as.data.frame(df176_sc14.002)

df176_sc14.003 <- df176_split[[2]]
df176_sc14.003 <- as.data.frame(df176_sc14.003)

df176_sc16.039 <- df176_split[[3]]
df176_sc16.039 <- as.data.frame(df176_sc16.039)

df176_sc16.040 <- df176_split[[4]]
df176_sc16.040 <- as.data.frame(df176_sc16.040)
###################################################
dfp24_M_sc15.033 <- dfp24_M_split[[1]]
dfp24_M_sc15.033 <- as.data.frame(dfp24_M_sc15.033)

dfp24_M_sc16.040 <- dfp24_M_split[[2]]
dfp24_M_sc16.040 <- as.data.frame(dfp24_M_sc16.040)

dfp24_M_sc17.041 <- dfp24_M_split[[3]]
dfp24_M_sc17.041 <- as.data.frame(dfp24_M_sc17.041)

dfp24_M_sc19.895 <- dfp24_M_split[[4]]
dfp24_M_sc19.895 <- as.data.frame(dfp24_M_sc19.895)
###################################################
dfp5_M_sc08.327 <- dfp5_M_split[[1]]
dfp5_M_sc08.327 <- as.data.frame(dfp5_M_sc08.327)

dfp5_M_sc15.033 <- dfp5_M_split[[2]]
dfp5_M_sc15.033 <- as.data.frame(dfp5_M_sc15.033)

dfp5_M_sc17.041 <- dfp5_M_split[[3]]
dfp5_M_sc17.041 <- as.data.frame(dfp5_M_sc17.041)
###################################################
dfp24_C_sc15.035 <- dfp24_C_split[[1]]
dfp24_C_sc15.035 <- as.data.frame(dfp24_C_sc15.035)

dfp24_C_sc17.042 <- dfp24_C_split[[2]]
dfp24_C_sc17.042 <- as.data.frame(dfp24_C_sc17.042)

dfp24_C_sc19.896 <- dfp24_C_split[[3]]
dfp24_C_sc19.896 <- as.data.frame(dfp24_C_sc19.896)
###################################################
dfp5_C_sc15.035 <- dfp5_C_split[[1]]
dfp5_C_sc15.035 <- as.data.frame(dfp5_C_sc15.035)

dfp5_C_sc16.040 <- dfp5_C_split[[2]]
dfp5_C_sc16.040 <- as.data.frame(dfp5_C_sc16.040)

dfp5_C_sc17.042 <- dfp5_C_split[[3]]
dfp5_C_sc17.042 <- as.data.frame(dfp5_C_sc17.042)
###################################################
dfARID_M_sc14.002 <- dfARID_M_split[[1]]
dfARID_M_sc14.002 <- as.data.frame(dfpARID_M_sc14.002)
###################################################
dfARID_C_sc14.003 <- dfARID_C_split[[1]]
dfARID_C_sc14.003 <- as.data.frame(dfARID_C_sc14.003)
#######################################################################################################################################################################################################################################
#I want to make csv's of the network frequencies for each experiment at an earlier filtering stage. Ideally this would be at F2, 
#but in this iteration of the metadata, I can only find F1. Thereofre I have:
#1) made a large list of data frames of each experiment from every cell line
#2)performed my network frequency function for networks that are at the F1 filtering stage
#3) written a for-loop that makes a csv of each and then puts it into the working directory

df.list <- list( 'dfH9_sc8.325_net_freq_F1' = dfH9_sc8.325, 'dfH9sc12.380_net_freq_F1' = dfH9_sc12.380,'dfH9_sc11.180_net_freq_F1'=dfH9_sc11.180,
                 'dfH9_sc11.182_net_freq_F1'=dfH9_sc11.182,"dfH9_sc11.181_net_freq_F1" = dfH9_sc11.181,'dfH9_sc10.949_net_freq_F1'=dfH9_sc10.949, 
                 'dfH9_sc12.382_net_freq_F1'=dfH9_sc12.382, 'dfH9_sc11.421_net_freq_F1'=dfH9_sc11.421,
                 'dfH9_sc13.767_net_freq_F1'=dfH9_sc13.767, 'dfH9_sc13.768_net_freq_F1'=dfH9_sc13.768, 'dfH9_sc11.179_net_freq_F1'=dfH9_sc11.179,
                 'dfH9_sc12.381_net_freq_F1'=dfH9_sc12.381,'dfH9120_sc06.692_net_freq_F1' = dfH9120_sc06.692, 'dfH9120_sc09.239_net_freq_F1' = dfH9120_sc09.239,
                 'df177_sc15.035_net_freq_F1' = df177_sc15.035, 'df177_sc12.380_net_freq_F1' = df177_sc12.380,
                 'df177_sc17.042_net_freq_F1' = df177_sc17.042,'df177_sc15.033_net_freq_F1' = df177_sc15.033,
                 'df177_sc17.041_net_freq_F1' = df177_sc17.041,'df177_sc12.382_net_freq_F1' = df177_sc12.382,
                 'df177_sc12.381_net_freq_F1' = df177_sc12.381,'df176_sc16.039_net_freq_F1' = df176_sc16.039, 'df176_sc16.040_net_freq_F1' = df176_sc16.040,
                 'df176_sc14.003_net_freq_F1' = df176_sc14.003,'df176_sc14.002_net_freq_F1' = df176_sc14.002,
                 'dfp24_M_sc19.895_net_freq_F1' = dfp24_M_sc19.895, 'dfp24_M_sc16.040_net_freq_F1' = dfp24_M_sc16.040,
                 'dfp24_M_sc17.041_net_freq_F1' = dfp24_M_sc17.041,'dfp24_M_sc15.033_net_freq_F1' = dfp24_M_sc15.033,
                 'dfp5_M_sc08.327_net_freq_F1' = dfp5_M_sc08.327, 'dfp5_M_sc15.033_net_freq_F1' = dfp5_M_sc15.033,
                 'dfp5_M_sc17.041_net_freq_F1' = dfp5_M_sc17.041,'dfp24_C_sc17.042_net_freq_F1' = dfp24_C_sc17.042, 'dfp24_C_sc19.896_net_freq_F1' = dfp24_C_sc19.896,
                 'dfp24_C_sc15.035_net_freq_F1' = dfp24_C_sc15.035,'dfp5_C_sc15.035_net_freq_F1' = dfp5_C_sc15.035, 'dfp5_C_sc17.042_net_freq_F1' = dfp5_C_sc17.042,
                 'dfp5_C_sc16.040_net_freq_F1' = dfp5_C_sc16.040,'dfARID_M_sc14.002_net_freq_F1' = dfARID_M_sc14.002,'dfARID_C_sc14.003_net_freq_F1' = dfARID_C_sc14.003)

df.list.filt <- lapply(df.list, function(x){
  y <- x %>%
    group_by(Network.Size.F1) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  y <- as.data.frame(y)
  y
  y$nNetwork <- y$n / y$Network.Size.F1
  #y <- filter(y, Network.Size.F4 >'1')
  y$freqperc <- 100*y$freq 
  y <- round_df(y, digits = 2)
  y
  
})
#Make the names of the df that comes out of the function the same as the original list
names(df.list.filt)<-names(df.list)
data_names <- names(df.list.filt)
#pull out all individual data frames from the final list and save to the global environment
list2env(df.list.filt,envir=.GlobalEnv)
#General For-loop info https://statisticsglobe.com/for-loop-in-r
#https://statisticsglobe.com/r-write-read-multiple-csv-files-for-loop
out <- "/Users/ramsey.najm/Desktop/GitHub/Projects/RN_CON/SC6-19_Metadata/Plots/NetworkFrequencies_byExp_F1/"
for(i in 1:length(data_names)) {                              # Head of for-loop
  print(i)
  write.csv(get(data_names[i]), 
            file =  paste0(out, data_names[i],".csv"), 
            row.names = FALSE) # Write CSV files to folder in directory specified by out
}
####################################################################################################################################################
#Generate Network Composition Heatmaps for each Exeriment
df.list <- list( 'dfH9_sc8.325' = dfH9_sc8.325, 'dfH9sc12.380' = dfH9_sc12.380,
                 'dfH9sc11.180'=dfH9_sc11.180, 'dfH9sc11.182'=dfH9_sc11.182,"dfH9sc11.181" = dfH9_sc11.181,
                 'dfH9sc10.949'=dfH9_sc10.949, 'dfH9sc12.382'=dfH9_sc12.382, 'dfH9sc11.421'=dfH9_sc11.421,
                 'dfH9sc13.767'=dfH9_sc13.767, 'dfH9sc13.768'=dfH9_sc13.768, 'dfH9sc11.179'=dfH9_sc11.179,
                 'dfH9sc12.381'=dfH9_sc12.381,'dfH9120_sc06.692' = dfH9120_sc06.692, 'dfH9120_sc09.239' = dfH9120_sc09.239,
                 'df177_sc15.035' = df177_sc15.035, 'df177_sc12.380' = df177_sc12.380,
                 'df177_sc17.042' = df177_sc17.042,'df177_sc15.033' = df177_sc15.033,
                 'df177_sc17.041' = df177_sc17.041,'df177_sc12.382' = df177_sc12.382,
                 'df177_sc12.381' = df177_sc12.381,'df176_sc16.039' = df176_sc16.039, 'df176_sc16.040' = df176_sc16.040,
                 'df176_sc14.003' = df176_sc14.003,'df176_sc14.002' = df176_sc14.002,
                 'dfp24_M_sc19.895' = dfp24_M_sc19.895, 'dfp24_M_sc16.040' = dfp24_M_sc16.040,
                 'dfp24_M_sc17.041' = dfp24_M_sc17.041,'dfp24_M_sc15.033' = dfp24_M_sc15.033,
                 'dfp5_M_sc08.327' = dfp5_M_sc08.327, 'dfp5_M_sc15.033' = dfp5_M_sc15.033,
                 'dfp5_M_sc17.041' = dfp5_M_sc17.041,'dfp24_C_sc17.042' = dfp24_C_sc17.042, 'dfp24_C_sc19.896' = dfp24_C_sc19.896,
                 'dfp24_C_sc15.035' = dfp24_C_sc15.035,'dfp5_C_sc15.035' = dfp5_C_sc15.035, 'dfp5_C_sc17.042' = dfp5_C_sc17.042,
                 'dfp5_C_sc16.040' = dfp5_C_sc16.040,'dfARID_M_sc14.002' = dfARID_M_sc14.002,'dfARID_C_sc14.003' = dfARID_C_sc14.003)

#Generate function and apply it to all data frames in that list
df.list.filt <- lapply(df.list, function(x){
  tab <- table(x$UVI.assignment.f4, x$RNA_snn_res.0.1.ordered)
  y <- as.data.frame.matrix(tab)
  #filter for networks greater than 1
  df <- y %>% filter(rowSums(.)>3)
  df <- as.matrix(df)
})
#Heatmap of network composition for every scSeq
names(df.list.filt)<-names(df.list)
data_names <- names(df.list.filt)
#Generate outgoing directory
out <- "/Users/ramsey.najm/Desktop/GitHub/Projects/RN_CON/SC6-19_Metadata/Plots/NetworkCompHM_by_Exp_over3/"
# make and save heat maps from the list of matrices generated in the above function
mapply(function(x, y){
  png(file=paste0(out, y, "_NetHM.png"))
  heatmap(x, Rowv = NULL, Colv= NA, 
          xlab = 'Seurat Cluster' , ylab = 'Network',col = palette, main= paste0(y,'_NetHM'))
  dev.off()
}, x = df.list.filt, y = data_names)
####################################################################################################################################################
