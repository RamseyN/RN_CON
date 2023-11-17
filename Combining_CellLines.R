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
df.list <- list( 'dfH9_netnumb_mat' = dfH9, 'dfH9120_netnumb_mat' = dfH9120,'df177_netnumb_mat' = df177, 'df176_netnumb_mat' = df176, 'dfp24_M_netnumb_mat' = dfp24_M,
                 'dfp24_C_netnumb_mat' = dfp24_C,'dfp5_M_netnumb_mat' = dfp5_M,'dfp5_C_netnumb_mat' = dfp5_C,'dfARID_M_netnumb_mat' = dfARID_M,
                 'dfARID_C_netnumb_mat' = dfARID_C)
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
#pull out all individual data frames from the final list and save to the global environment
list2env(df.list.filt,envir=.GlobalEnv)
#Correct TSCp5_C to have the same number of columsn as TSCp24_C (missing column 6)
is.data.frame(dfp5_C_netnumb_mat)
is.matrix(dfp5_C_netnumb_mat)
#Make column 6 and fill it with 0's
dfp5_C_netnumb_mat <- cbind(dfp5_C_netnumb_mat, "6" = 0)
dfp5_C_netnumb_mat
#switch the orddr so column 6 is in the 7th column (including 0)
# Get the number of columns in the matrix
num_columns <- ncol(dfp5_C_netnumb_mat)
num_columns
# Check if the matrix has at least 11 columns
if (num_columns >= 11) {
  # Extract the 11th column
  column_11 <- dfp5_C_netnumb_mat[, 11]
  
  # Remove the 11th column from the matrix
  dfp5_C_netnumb_mat <- dfp5_C_netnumb_mat[, -11]
  
  # Insert the 11th column at the 7th position
  dfp5_C_netnumb_mat <- cbind(dfp5_C_netnumb_mat[, 1:6], column_11, dfp5_C_netnumb_mat[, 7:(num_columns - 1)])
  
  # Print the result
  print(dfp5_C_netnumb_mat)
} else {
  cat("Matrix must have at least 11 columns.")
}

# Rename the column
colnames(dfp5_C_netnumb_mat)[which(colnames(dfp5_C_netnumb_mat) == "column_11")] <- "6"

# Print the result
print(dfp5_C_netnumb_mat)

#make heatmap of combined TSC Cell Lines

df_list <- list('TSCp24C' = dfp24_C_netnumb_mat, 'TSCp5C' = dfp5_C_netnumb_mat)
lapply(df_list, colnames)

# Combine the dataframes row-wise
combined_df <- do.call(rbind, df_list)
heatmap(combined_df, Rowv = NULL, Colv= NA, 
        xlab = 'Seurat Cluster' , ylab = 'Network',col = palette, main = 'TSC CONTROL NETWORK >3')

###NOW DO IT FOR THE MUTANTS#############################MUCH EASIER#######################################

df_list <- list('TSCp24M' = dfp24_M_netnumb_mat, 'TSCp5M' = dfp5_M_netnumb_mat)


# Combine the dataframes row-wise
combined_df <- do.call(rbind, df_list)
heatmap(combined_df, Rowv = NULL, Colv= NA, 
        xlab = 'Seurat Cluster' , ylab = 'Network',col = palette, main= 'TSC MUTANT NETWORK >3')

