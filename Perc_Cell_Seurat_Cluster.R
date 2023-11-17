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
########################################################################################################################################################################
########################################################################################################################################################################
#Heatmaps of network composition greater than 1 across cell lines
df.list <- list( 'dfH9_colsum' = dfH9, 'dfH9120_colsum' = dfH9120,'df177_colsum' = df177, 'df176_colsum' = df176, 'dfp24_M_colsum' = dfp24_M,
                 'dfp24_C_colsum' = dfp24_C,'dfp5_M_colsum' = dfp5_M,'dfp5_C_colsum' = dfp5_C,'dfARID_M_colsum' = dfARID_M,
                 'dfARID_C_colsum' = dfARID_C)
#Generate function and apply it to all data frames in that list
df.list.filt <- lapply(df.list, function(x){
  tab <- table(x$UVI.assignment.f4, x$RNA_snn_res.0.1.ordered)
  y <- as.data.frame.matrix(tab)
  #filter for networks greater than X
  df <- y %>% filter(rowSums(.)>3)
  # Calculate the sum of each column
  column_sums <- colSums(df)
  # Calculate the total sum of all values in the data frame
  total_sum <- sum(unlist(df))
  # Calculate the percentage for each column
  column_percentages <- (column_sums / total_sum) * 100
  as.data.frame(column_percentages)
})

#make names of df.list.filt the same as the names of df.list
names(df.list.filt)<-names(df.list)
#assign a variable to the names of df.list.filt
data_names <- names(df.list.filt)
#pull out all individual data frames from the final list and save to the global environment
list2env(df.list.filt,envir=.GlobalEnv)


out <- "/Users/ramsey.najm/Desktop/GitHub/Projects/RN_CON/SC6-19_Metadata/Plots/Network_Cluster_perc/"
for(i in 1:length(data_names)) {                              # Head of for-loop
  print(i)
  write.csv(get(data_names[i]), 
            file =  paste0(out, data_names[i],"_perc.csv"), 
            row.names = TRUE) # Write CSV files to folder in directory specified by out
}



