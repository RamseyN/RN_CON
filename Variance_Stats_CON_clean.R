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
#Ploting the variance of the network composition across cell lines (Number of Cells/Cluster)
df_list <- list( 'H9_netnumb' = dfH9, 'H9120_netnumb' = dfH9120,'177_netnumb' = df177, '176_netnumb' = df176, 'TSCp24M_netnumb' = dfp24_M,
                'TSCp24C_netnumb' = dfp24_C,'TSCp5M_netnumb' = dfp5_M,'TSCp5C_netnumb' = dfp5_C,'ARIDM_netnumb' = dfARID_M,
                'ARIDC_netnumb' = dfARID_C)
#Generate function and apply it to all data frames in that list
df_list <- lapply(df_list, function(x){
  tab <- table(x$UVI.assignment.f4, x$RNA_snn_res.0.1.ordered)
  y <- as.data.frame.matrix(tab)
  #filter for networks greater than 1
  df <- y %>% filter(rowSums(.)>3)
})
#Calculate Variance/Row
df_list <- lapply(df_list, function(x){
  x$RowVariance <- apply(x, 1, var)
  x
})

#Make the names of df_list.filt = names of df_list
names(df_list)<-names(df_list)
#pull out all individual data frames from the final list and save to the global environment
list2env(df_list,envir=.GlobalEnv)
# Identify unique column names
unique_column_names <- unique(unlist(lapply(df_list, names)))
# Ensure all dataframes have the same column structure
for (i in 1:length(df_list)) {
  missing_columns <- setdiff(unique_column_names, names(df_list[[i]]))
  extra_columns <- setdiff(names(df_list[[i]]), unique_column_names)
  
  # Add missing columns with NA values
  for (col in missing_columns) {
    df_list[[i]][, col] <- NA
  }
  
  # Remove extra columns not present in other dataframes
  if (length(extra_columns) > 0) {
    df_list[[i]] <- df_list[[i]][, unique_column_names]
  }
}

# Combine the dataframes row-wise
combined_df <- do.call(rbind, df_list)


# Check if there are any NA values in the entire dataframe
has_na <- any(is.na(combined_df))

# Output the results
if(has_na) {
  print("There are NA values in the dataframe.")
} else {
  print("There are no NA values in the dataframe.")
}

#Replace all NA values (This is likely because of cluster 10 is only represented in some cell lines 
#for the cell lines that don't have cells in cluster 10, I replace these with a 0 value)
# Assuming 'df' is your dataframe
combined_df <- combined_df %>% replace(is.na(.), 0)

# Check if there are any NA values in the entire dataframe
has_na <- any(is.na(combined_df))

# Output the results
if(has_na) {
  print("There are NA values in the dataframe.")
} else {
  print("There are no NA values in the dataframe.")
}


#Try this measure of row heterogeneity...
#This code snippet calculates the sum of squared differences for each row in the dataframe df. 
#It subtracts the mean of each row from each element in the row, squares these differences, and 
#then sums them up. Higher values of row_heterogeneity would indicate rows where one or more values 
#significantly deviate from the mean value of that row.
# Select columns excluding 12 and 13
#Initially I tried to remove the grop rows, that's before I figured out how to use selected columns

#combined_df <- combined_df[, !names(combined_df) %in% "Group"]
combined_df
# Example dataframe df
selected_columns <- combined_df[, c('0','1','2','3','4','5','6','7','8','9','10')]#Selecting for commonly populated columns
epsilon <- 0.0001  # Small value to avoid issues with 0 values

combined_df$RowHet <- apply(selected_columns, 1, function(row) sum((row - mean(row) + epsilon)^2))
combined_df
#Now try calculating the entropy from full columns
combined_df$RowEntropy <- apply(selected_columns, 1, function(row) entropy(row))
combined_df$RowEntropyML <- apply(selected_columns, 1, function(row) entropy(row, method='ML'))
combined_df$RowEntropyshrink <- apply(selected_columns, 1, function(row) entropy(row, method='shrink'))
combined_df$RowEntropyJeffreys <- apply(selected_columns, 1, function(row) entropy(row, method='Jeffreys'))

# Print the updated dataframe
print(combined_df)


# Extract the grouping variable from row names before the first underscore
grouping_variable <- sub("^(.*?)_.*", "\\1", rownames(combined_df))

# Add the grouping variable as a new column
combined_df <- combined_df %>%
  mutate(Group = grouping_variable)

# Define the custom order for the groups
custom_order <- c("H9", "H9120", "177", "176", "TSCp24M", "TSCp24C", "TSCp5M", "TSCp5C", "ARIDM", "ARIDC")

# Reorder the levels of the Group variable
combined_df$Group <- factor(combined_df$Group, levels = custom_order)

#Rearrange the GROUP column to be the first in the dataframe to make it more easily readible
# Reorder columns, making the current column 13 as column 1
combined_df <- combined_df %>% select(Group, everything())
# Print the updated dataframe
print(combined_df)


# Create a boxplot with grouping variable 'Group'
png(file =  "/Users/ramsey.najm/Desktop/GitHub/Projects/RN_CON/SC6-19_Metadata/Plots/Statistics/NetComp_Variance_CellLine_over3.png",
   width = 12, height = 6, units = "in", res = 300)
boxplot(RowVariance ~ Group, data = combined_df, outline = FALSE,
        main = "Network >3 Composition Variance",
        xlab = "Cell Line", ylab = "Variance")
dev.off()


####################################################################################################################################################
# Create the boxplot with individual data points
# p <- ggplot(combined_df, aes(x = Group, y = RowVariance)) +
#   geom_boxplot() +
#   geom_jitter(aes(x = Group, y = RowVariance), width = 0.2, alpha = 0.5) +  # Add jittered data points
#   labs(title = "Network >3 Composition Variance", y = "Variance in Network Composition (Cell-type Percentage") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# # Limit the number of data points shown to a maximum of 100 per group
# p <- p + coord_cartesian(ylim = c(min(combined_df$RowVariance), max(combined_df$RowVariance))) 
# 
# # Print the plot
# png(file =  "/Users/ramsey.najm/Desktop/GitHub/Projects/RN_CON/SC6-19_Metadata/Plots/Statistics/NetComp_Variance_CellLine_over3.png",
#     width = 12, height = 6, units = "in", res = 300)
# print(p)
# dev.off()

#This can be problematic with outliers, so just use normal boxplot for now

####################################################################################################################################################

boxplot(RowVariance ~ Group, data = combined_df, outline = FALSE,
        main = "Network >3 Composition Variance",
        xlab = "Cell Line", ylab = "Variance")

boxplot(RowHet ~ Group, data = combined_df, outline = FALSE,
        main = "Network >3 Composition Heterogeneity",
        xlab = "Cell Line", ylab = "Heterogeneity")

boxplot(RowEntropy ~ Group, data = combined_df, outline = FALSE,
        main = "Network >3 Composition Entropy",
        xlab = "Cell Line", ylab = "Entropy")





