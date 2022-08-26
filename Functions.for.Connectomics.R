######################################################################
# Functions.for.Connectomics.R
######################################################################
# source('~/GitHub/Projects/CON/functions/Functions.for.Connectomics.R')
# rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)


any.duplicated <- function(vec, summarize = TRUE) { # How many entries are duplicated?.
  y = sum(duplicated(vec))
  if (summarize & y) {
    x = table(vec); x = x[x > 1] - 1;
    print("The following elements have  > 1 extra copies:")
    print(x) # table formatting requires a separate entry
  }
  return(y)
}




#  ------------------------

FixUnderscores <- function(string = "stairway__to_heaven_", trimFinal = TRUE) { #
  string <- gsub(x = string, pattern = '_+', replacement = '_')
  LastChr <- substr(string, nchar(string), nchar(string))
  if (trimFinal && LastChr == "_") {
    iprint('LastChr: ', LastChr)
    string = substr(string, 1, (nchar(string)-1))
  }
  return(string)
}



#  ------------------------



Calc.Network.Size <- function(obj = combined.obj, col = 'UVI.assignment') {
  obj[[col]] %>% 
    na.omit.strip() %>% 
    table
}


Assign.Network.Size.Metadata <- function(obj = combined.obj, plot.Size.Distr= TRUE
                                         , input.col = 'UVI.assignment', output.col = 'Network.Size' ) {
  NetworkSizeTable <- Calc.Network.Size(obj = obj, col = input.col)  
  if (plot.Size.Distr) qhistogram(NetworkSizeTable)
  NetwSizePerCell <- as.numeric(translate(vec = deframe(obj[[input.col]])
                                          , oldvalues = names(NetworkSizeTable)
                                          , newvalues = NetworkSizeTable))
  obj[[output.col]] <- NetwSizePerCell
  obj
}


#  ------------------------

fr.topX.pc <- function(vec = UVI.table.147135$n, quantile.thr = .99) {
  nr.top.X.pc <- round(l(vec) *(1-quantile.thr))
  sum.topX <- sum(sort(vec, decreasing = T)[1:nr.top.X.pc])
  sum.topX / sum(vec)
}
# fr.topX.pc()

fr.UVI.df <- function(dfx = UVI.table.147135, col = "n") {
  dfx$'percent' <- dfx[,col] / sum(dfx[,col])
  dfx
}
# head(fr.UVI.df())
fr.UVI.df2 <- function(dfx = UVI.table.147135, col = "n") {
  dfx <- dfx[,col] / sum(dfx[,col])
  dfx
}


fr.UVI.df2 <- function(dfx = UVI.table.147135, col = "n") {
  dfx$'rank' <- seq(1:nrow(dfx))
  dfx$'percent' <- dfx[,col] / sum(dfx[,col])
  dfx$'cumul.frac' <- cumsum(dfx$'percent')
  tbl <- tibble(dfx[, c('rank', 'percent', 'cumul.frac')])
  rownames(tbl) <- rownames(dfx)
  tbl
}

pwr <- function (x, exponent = -1.7)   x^exponent


merge_1col_dfs_by_rn <- function(list_of_dfs, FILLwith = 0,columnUSE= 'n') {
  all.rn <- sort(union.ls(lapply(list_of_dfs, rownames)))
  iprint("n rownames:",l(all.rn))
  xx <- data.frame(matrix(data = FILLwith, nrow = l(all.rn), ncol = l(list_of_dfs)), row.names = all.rn)
  colnames(xx) <-  names(list_of_dfs)
  for (i in 1:l(list_of_dfs)) {
    print(i)
    indf <- list_of_dfs[[i]]
    xx[rownames(indf),i] <- indf[,columnUSE]
  }
  xx
}
# merge_1col_dfs_by_rn()

#  ------------------------


#  ------------------------


#  ------------------------



theme_sky<- function (base_size = 16, base_family = ""){
  theme_void() %+replace% 
    theme(
      line = element_line(colour = "white", size = 0.5, linetype = 1, 
                          lineend = "butt") 
      , text = element_text(family = base_family, 
                            face = "plain", colour = "white", size = base_size,
                            angle = 0, lineheight = 0.9, hjust = 0, vjust = 0)
      , plot.background = element_rect(colour = '#151E3D', fill = '#151E3D')
      , plot.title = element_text(size = rel(1.2))
      , panel.grid.major = element_line(colour = "grey20", size = 0.2)
      , panel.grid.minor = element_line(colour = "grey5", size = 0.5)
      , strip.background = element_rect(fill = "grey30", colour = "grey30")
    )
}

#  ------------------------
AddGOGeneList.manual <- function(obj = combined.obj, GO = 'GO:0034976', web.open=F  # Add GO terms via Biomart package.
                                 , genes =  c("A0A140VKG3", "ARX", "CNTN2", "DRD1", "DRD2", "FEZF2", "LHX6")) {
  print(head(genes, n = 15))
  genes <- IntersectWithExpressed(obj = obj, genes = genes)
  
  if (is.null(obj@misc$GO)) obj@misc$GO <- list()
  obj@misc$GO[[make.names(GO)]] <- genes
  MarkdownReportsDev::iprint("Genes in", GO, "are saved under obj@misc$GO$", make.names(GO))
  if (web.open) system(paste0("open https://www.ebi.ac.uk/QuickGO/search/", GO))
  return(obj)
}