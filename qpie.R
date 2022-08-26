qpie <- function(vec, ext = "png", plot = TRUE, save = TRUE, mdlink = FALSE
                 , suffix = NULL
                 , plotname = sppp(substitute(vec), suffix)
                 , LegendSide = T, LegendTitle = as.character(substitute(vec)), NoLegend = F
                 , pcdigits = 2, NamedSlices =F
                 , custom.order = F
                 # , custom.margin = F
                 , palette_use = 'jco'
                 , max.names = 30
                 , max.categories = 1000
                 , w = 5, h = w, ...) {
  
  print(plotname)
  
  if (length(vec) > max.categories) {
    iprint("Warning, there are more than", max.categories, "categories. Only the top", max.categories - 1, "items are show, the rest is added up.")
    sv <- sort(vec, decreasing = T)
    vec.new <- sv[1:(max.categories - 1)]
    idx.remaining <- max.categories:length(vec)
    sum.of.remaining <- sum(sv[idx.remaining])
    fr.sum <- percentage_formatter(sum.of.remaining / sum(vec))
    iprint("The remaining", length(idx.remaining), "values make up", fr.sum,"of the data.")
    
    vec.new[max.categories] <- sum.of.remaining
    vec <- vec.new
    
  }
  
  
  
  if (is_null(names(vec))) { names(vec) <- as.character(1:length(vec)) }
  
  
  
  df <- qqqCovert.named.vec2tbl(namedVec = vec, thr = max.names)
  print(df)
  nrCategories.DFcol1 <- length(unique(df[,1])); stopif( nrCategories.DFcol1 > max.categories)
  
  
  
 print(nrCategories.DFcol1)
  
  
  
  pcX <- df$"value" / sum(df$"value")
  labs <- paste(100 * signif(pcX, pcdigits), "%", sep = "")
  if (NamedSlices) labs <- paste(df$names, "\n", labs)
  if (custom.order != F) df$'names' <- factor(df$'names', levels = custom.order)
  
  
  
  p <- ggpubr::ggpie(data = df, x = "value", label = labs
                     , fill = "names", color = "white"
                     , title = plotname
                     , palette = palette_use, ...)
  if (LegendSide) p <- ggpubr::ggpar(p, legend = "right", legend.title = LegendTitle)
  # if (custom.margin) p <- p + theme(plot.margin = unit(custom.margin, "cm"))
  # p <- if (NoLegend) p + NoLegend() else p
  p <- if (NoLegend) p + theme(legend.position = "none", validate = TRUE) else p
  fname = Stringendo::kpp(plotname, "pie", ext)
  if (save) qqSave(ggobj = p, title = plotname, fname = fname, ext = ext, w = w, h = h)
  if (mdlink & save) qMarkdownImageLink(fname)
  if (plot) p
}




