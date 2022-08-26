#Clean Pipeline for Assessing Rabies Complexity and Bias from Amplicon Sequencing
##############################Source Code
#source('~/GitHub/Projects/RN_CON/Scripts/Load.packages.local.R')

require(DataInCode)    # https://github.com/vertesy/DataInCode
require(DatabaseLinke.R)    # https://github.com/vertesy/DatabaseLinke
require(CodeAndRoll2)    # https://github.com/vertesy/CodeAndRoll2
require(Stringendo)    # https://github.com/vertesy/Stringendo
require(ReadWriter)    # https://github.com/vertesy/ReadWriter
require(ggExpress)    # https://github.com/vertesy/ggExpress
require(MarkdownReports)    # https://github.com/vertesy/MarkdownReports
require(MarkdownHelpers)
require(Seurat.utils)    # https://github.com/vertesy/Seurat.utils
require(colorout)
require(Seurat)
require(future) # parallelization
require(doMC)
require(tictoc)
require(schex)
require(clipr)
require(ggplot2)
library(moments)
require(tidyverse)
source('~/GitHub/Projects/RN_CON/Scripts/Rocinante.R') #warning for AddGOGeneList.maual.R and IntersetWithExpressed.R
source('~/GitHub/Projects/RN_CON/Scripts/Functions.for.Connectomics.R')
source('~/GitHub/Projects/RN_CON/Scripts/ggAuxillaryFunctions.R') #essential for qPie to work
source('~/GitHub/Projects/RN_CON/Scripts/qpie.R') #Pie Plots
#Setup------------------------------------------------------------------------
setup_MarkdownReports(OutDir = '~/Dropbox (VBC)/Group Folder Knoblich/Users/Ramsey/Connectomics/Analysis/RV_Amplicon_QC/', scriptname = "RV_Quality_Control_RN.R")
OutDirOrig = OutDir
# Read In ------------------------
#Make a list of all the UVI_Counts we want to plot
ls.df <- list(
  'pGFP1' <- read.simple.tsv('/Volumes/groups/knoblich/users/burkard/Abel.Vertesy/R11846/164962/UVI_count.tsv'), #GFP-Plasmid UVI Rep 1
  'pGFP2' <- read.simple.tsv('/Volumes/groups/knoblich/users/burkard/Abel.Vertesy/R11851/164978/UVI_count.tsv'), #GFP-Plasmid UVI Rep 2
  'pCrim' <- read.simple.tsv('/Volumes/groups/knoblich/users/burkard/Abel.Vertesy/R12357_bulk/173978/UVI_count.tsv'), #Crimson-Plasmid UVI
  'pCrimGG14' <- read.simple.tsv('/Volumes/groups/knoblich/users/burkard/Abel.Vertesy/R13959_R14001/204098/UVI_count.tsv'), #M barcoded (GG14) Crimson-Plasmid UVI
  'vGFP' <- read.simple.tsv('/Volumes/groups/knoblich/users/burkard/Abel.Vertesy/R12357_bulk/173979/UVI_count.tsv'), #GFP-Virus UVI
  'vCrim' <- read.simple.tsv('/Volumes/groups/knoblich/users/burkard/Abel.Vertesy/R12357_bulk/173980/UVI_count.tsv'), #Crimson-Virus UVI
  'HEKGFP1' <- read.simple.tsv('/Volumes/groups/knoblich/users/burkard/Abel.Vertesy/R11846/164961/UVI_count.tsv'), #GFP-HEK UVI Rep 1
  'HEKGFP2' <- read.simple.tsv('/Volumes/groups/knoblich/users/burkard/Abel.Vertesy/R11851/164977/UVI_count.tsv'), #GFP-HEK UVI Rep 2
  'HEKCrim' <- read.simple.tsv('/Volumes/groups/knoblich/users/burkard/Abel.Vertesy/R12357_bulk/173981/UVI_count.tsv'), #Crimson-HEK UVI
  'HEKCrim_RV6' <- read.simple.tsv('/Volumes/groups/knoblich/users/burkard/Abel.Vertesy/R13062/187475/UVI_count.tsv') 
)

#Add Metadata
Samples <- c('pGFP1','pGFP2','pCrim','pCrimGG14','vGFP','vCrim', 'HEKGFP1', 'HEKGFP2','HEKCrim', 'HEKCrim_RV6')
SampleTypes <- c('Plasmid','Plasmid','Plasmid','Plasmid','Virus','Virus','HEK.Cells', 'HEK.Cells','HEK.Cells','HEK.Cells')
Fluorophore <- c('GFP-UVI','GFP-UVI','Crimson-UVI','Crimson-UVI', 'GFP-UVI', 'Crimson-UVI','GFP-UVI','GFP-UVI', 'Crimson-UVI','Crimson-UVI')
SampleSpec <- kpp(SampleTypes, Fluorophore)

Samples.Plasmid <- Samples[SampleTypes ==  'Plasmid']
Samples.RabV <- Samples[SampleTypes ==  'HEK.Cells']
Samples.RabV.Raw <- Samples[SampleTypes ==  'Virus']

#Start Statistics
(ls.n.vec <- lapply(ls.df, as.named.vector))
N.Barcodes.Detected <- unlapply(ls.n.vec, l)
Depth.VS.N.Barcodes.Detected <- tibble(Depth, N.Barcodes.Detected, Samples, SampleTypes)
Depth.VS.N.Barcodes.Detected$SampleTypes <- factor(Depth.VS.N.Barcodes.Detected$SampleTypes, 
                             levels = c('Plasmid','Virus','HEK.Cells'))
#Bar Plots
{
  create_set_Original_OutDir(); create_set_SubDir("Barplots")
{
p <- ggbarplot(Depth.VS.N.Barcodes.Detected, 'Samples','N.Barcodes.Detected'
          , fill = 'SampleTypes', color = 'SampleTypes', label = TRUE)+ 
    theme(text = element_text(size = 16))
  png('Depth_vs_N_Barcodes_Detected.png',width = 1100, height = 500)
  print(p)
  dev.off()
}
#Simple barplot laying depicting N.Barcodes Detected by Samples
##############################Barplots of bias###########################


"Skew"
Skew <- unlapply(ls.n.vec, skewness)
Skewdf <- tibble(Skew, Samples, SampleTypes)
Skewdf$SampleTypes <- factor(Skewdf$SampleTypes, 
                              levels = c('Plasmid','Virus','HEK.Cells'))
{
  p <- ggbarplot(Skewdf, 'Samples','Skew'
          , fill = 'SampleTypes', color = 'SampleTypes', label = TRUE
          , lab.nb.digits = 1, title = 'Skew')+ 
    theme(text = element_text(size = 16))
png('Skew.png', width = 1100, height = 500)
print(p)
dev.off()
}

##############################Measure of Inequality###########################
"Inequality - Gini" #How do I labele each column on X axis with sample names?
Gini <- unlapply(ls.n.vec, DescTools::Gini)
Ginidf <- tibble(Gini, Samples, SampleTypes)
Ginidf$SampleTypes <- factor(Ginidf$SampleTypes, 
                             levels = c('Plasmid','Virus','HEK.Cells'))
{
  p <- ggbarplot(Ginidf, 'Samples','Gini'
                 , fill = 'SampleTypes', color = 'SampleTypes', label = TRUE
                 , lab.nb.digits = 3,title = 'Gini')+ 
    theme(text = element_text(size = 16))
  png('Gini.png',width = 1100, height = 500)
  print(p)
  dev.off()
}

##############################Fraction of reads that make up top 1% of all UVIs###########################
"Fraction top 1%"
top.1.pc.UVIs <-unlapply(ls.n.vec, fr.topX.pc) *100
top.1.pc.UVIsdf <- tibble(top.1.pc.UVIs, Samples, SampleTypes)
top.1.pc.UVIsdf$SampleTypes <- factor(top.1.pc.UVIsdf$SampleTypes, 
                             levels = c('Plasmid','Virus','HEK.Cells'))
{
  p <- ggbarplot(top.1.pc.UVIsdf, 'Samples','top.1.pc.UVIs'
                 , fill = 'SampleTypes', color = 'SampleTypes', label = TRUE
                 , lab.nb.digits = 1,title = 'top.1.pc.UVIs')+ 
    theme(text = element_text(size = 16))
  png('top.1.pc.UVIs.png', width = 1100, height = 500)
  print(p)
  dev.off()
}

"Fraction top 10%"
top.10.pc.UVIs <-unlapply(ls.n.vec, fr.topX.pc, quantile.thr = .9) *100
top.10.pc.UVIsdf <- tibble(top.10.pc.UVIs, Samples, SampleTypes)
top.10.pc.UVIsdf$SampleTypes <- factor(top.10.pc.UVIsdf$SampleTypes, 
                                      levels = c('Plasmid','Virus','HEK.Cells'))
{
  p <- ggbarplot(top.10.pc.UVIsdf, 'Samples','top.10.pc.UVIs'
                 , fill = 'SampleTypes', color = 'SampleTypes', label = TRUE
                 , lab.nb.digits = 1,title = 'top.10.pc.UVIs')+ 
    theme(text = element_text(size = 16))
  png('top.10.pc.UVIs.png', width = 1100, height = 500)
  print(p)
  dev.off()
}
}; create_set_Original_OutDir()

##################################################################################################

#Scatter Plots
Depth <- unlapply(ls.n.vec, sum)
Depth.VS.N.Barcodes.Detecteddf <- tibble(Depth, N.Barcodes.Detected, Samples, SampleTypes)
Depth.VS.N.Barcodes.Detecteddf$SampleTypes <- factor(Depth.VS.N.Barcodes.Detecteddf$SampleTypes, 
                                       levels = c('Plasmid','Virus','HEK.Cells'))
{
  create_set_Original_OutDir(); create_set_SubDir("ScatterPlots")
{
s <- ggscatter(Depth.VS.N.Barcodes.Detecteddf, x = 'Depth', y = 'N.Barcodes.Detected',
          xlab = 'Sequencing Depth',
          ylab = 'Number of Barcodes Detected',
          label = 'Samples',
          color = 'SampleTypes',
          show.legend = F,
          show.legend.text = F,
          font.label = c(16, 'plain'))+
  # setting legend for color https://github.com/kassambara/ggpubr/issues/111
  scale_color_discrete(
    name = '', 
    labels = c('Plasmid','Virus','HEK Cells'))+
  grids(linetype = 'solid')+ 
    theme(text = element_text(size = 20))

  png('scatter.png', width = 1500, height = 500)
  print(s)
  dev.off()
}
  }; create_set_Original_OutDir()



#-----------------------------------------Distribution Plots------------------------------
(ls.frac <- lapply(ls.df, FUN = fr.UVI.df2))
names(ls.frac) <- Samples

ls.frac.H <- ls.frac[Samples[SampleTypes == 'HEK.Cells']]
ls.frac.V <- ls.frac[Samples[SampleTypes == 'Virus']]
ls.frac.Pl <- ls.frac[Samples[SampleTypes == 'Plasmid']]

{
  create_set_Original_OutDir(); create_set_SubDir("Distribution_Plots")
  
  # RABIES HEK Distribution Plot---------------------------------------------------------------------------------
{
  (dfl.fraction <- tibble(reshape2::melt(data = ls.frac.H, value.name = c("Fraction"), id.vars = "rank"))[, -2])
  colnames(dfl.fraction) <- c("UVI Rank", "UVI Fraction", "HEK experiment")
  
  
  plot.UVI.Freq.HEK <- ggscatter(dfl.fraction, x = 'UVI Rank', y = 'UVI Fraction', color = 'HEK experiment'
                                 , size = 0.1, title = 'RabV UVI Distribution HEK Cell' 
                                 , ylim = c(1e-7,1)) +
    xscale('log10') +
    yscale('log10') +
    grids(color = 'grey') +
    stat_function(fun = pwr, aes(color='Power Law (x^-1.7)')) 
  #stat_smooth(method = 'nls', formula = 'y~a*exp(b*x)', method.args = list(start=c(a=1, b=1)), se=FALSE, aes(color="exponential"), size =0.5) + # https://stackoverflow.com/questions/38378161/how-to-plot-non-linear-regression-lines-within-groups-and-total-data-in-ggplot2
  #stat_smooth(method = 'nls', formula = 'y~a*x^b', method.args = list(start=c(a=1,b=1)), se=FALSE, aes(color="pwr.law"), size =0.5)
  plot.UVI.Freq.HEK
  png('HEK_Freq.png', width = 1000, height = 500)
  print(plot.UVI.Freq.HEK)
  dev.off()
  
}
  
  {(dfl.fraction <- tibble(reshape2::melt(data = ls.frac.V, value.name = c("Fraction"), id.vars = "rank"))[, -2])
    colnames(dfl.fraction) <- c("UVI Rank", "UVI Fraction", "Experiment")
    
    
    plot.UVI.Freq.Virus <- ggscatter(dfl.fraction, x = 'UVI Rank', y = 'UVI Fraction', color = 'Experiment'
                                     , size = 0.1, title = 'RabV UVI Distribution in Library Virus' #RN Changed title
                                     , ylim = c(1e-7,1)) +
      xscale('log10') +
      yscale('log10') +
      grids(color = 'grey') +
      stat_function(fun = pwr, aes(color='Power Law (x^-1.7)')) 
    #stat_smooth(method = 'nls', formula = 'y~a*exp(b*x)', method.args = list(start=c(a=1, b=1)), se=FALSE, aes(color="exponential"), size =0.5) + # https://stackoverflow.com/questions/38378161/how-to-plot-non-linear-regression-lines-within-groups-and-total-data-in-ggplot2
    #stat_smooth(method = 'nls', formula = 'y~a*x^b', method.args = list(start=c(a=1,b=1)), se=FALSE, aes(color="pwr.law"), size =0.5)
    plot.UVI.Freq.Virus
    png('Virus_Freq.png', width = 1000, height = 500)
    print(plot.UVI.Freq.Virus)
    dev.off()
    }  
  
  {ls.frac.Pl <- lapply(ls.frac.Pl, head, 5e5)
    "Take first 500K BCs (most freq)."
    
    (dfl.fraction <- tibble(reshape2::melt(data = ls.frac.Pl, value.name = c("Fraction"), id.vars = "rank"))[, -2])
    colnames(dfl.fraction) <- c("UVI Rank", "UVI Fraction", "Plasmid Prep")
    
    
    
    plot.UVI.Freq.Plasmid <- ggscatter(dfl.fraction, x = 'UVI Rank', y = 'UVI Fraction', color = 'Plasmid Prep'
                                       , size = 0.1, title = 'RabV UVI Distribution in Plasmid Libraries'
                                       , ylim = c(1e-7, 1e0)
    ) +
      xscale('log10') +
      yscale('log10') +
      grids(color = 'grey') +
      stat_function(fun = pwr, aes(color='Power Law (x^-1.7)')) 
    #stat_smooth(method = 'nls', formula = 'y~a*exp(b*x)', method.args = list(start=c(a=1, b=1)), se=FALSE, aes(color="exponential"), size =0.5) + # https://stackoverflow.com/questions/38378161/how-to-plot-non-linear-regression-lines-within-groups-and-total-data-in-ggplot2
    #stat_smooth(method = 'nls', formula = 'y~a*x^b', method.args = list(start=c(a=1, b=1)), se=FALSE, aes(color="pwr.law"), size =0.5)
    plot.UVI.Freq.Plasmid}  
  png('Plasmid_Freq.png', width = 1000, height = 500)
  print(plot.UVI.Freq.Plasmid)
  dev.off()

}; create_set_Original_OutDir()
#--------------------------------Pie Plots------------------------------------------------
{
create_set_Original_OutDir(); create_set_SubDir("Pie_Plots")


UVI.dis.in.pCrimGG14 <- (pCrimGG14$n)
qpie(vec = UVI.dis.in.pCrimGG14)

  
}; create_set_Original_OutDir()


