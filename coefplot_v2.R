library(gdata)
library(DOQTL)
library(GenomicRanges)
library(VariantAnnotation)
library(foreach)
library(doParallel)
#library(abind)

prob.plot <- function(pheno, pheno.col, probs, marker) {
  ## Remove NAs from pheno data
  keep = rownames(pheno[!is.na(pheno[,pheno.col]),])
  if (length(keep) < 3) {
    stop("Too many missing phenotype values!")
  }
  probmat = probs[keep, , marker]
  
  image(z=probmat, zlim=c(0, 2), col=rev(gray(seq(from=0, to=1, by=0.01))), axes=FALSE)
  states <- colnames(probmat)
  nk <- ncol(probmat)
  axis(2, at=seq(from=0, to=1, len=nk), labels=states, las=1, tck=0)
  
  ## Calculate z-scores
  y = pheno[keep, pheno.col]
  z <- (y - mean(y))/sd(y)
  
  axis(3, at=ecdf(z)(pretty(z)), labels=pretty(z))
  mtext("z-score scale", side=3, at=0.5, line=2)
  xpos <- c(0, mean(y<mean(y)), 1)
  xtxt <- signif(c(min(y), mean(y), max(y)), 4)
  xtxt <- ifelse(abs(xtxt) < 1e-8, 0, xtxt)
  axis(1, at=xpos, labels=xtxt)
  box()
}
attr(prob.plot, 'help') = "
This function creates a probability plot for the mapping population.

Parameters:
pheno: The phenotype dataframe
pheno.col: The column in the pheno dataframe to use as the phenotype
probs: The 3D probability array
marker: The ID of the marker to plot

"




coefplot_v2 = function(doqtl, chr = 1, start=NULL, end=NULL, stat.name = "LOD", remove.outliers=50, conf.int = FALSE, legend = TRUE, colors = "DO", sex, ...) {
  
  old.par = par(no.readonly = TRUE)
  
  cross = attr(doqtl, "cross")
  if(is.null(cross)) {
    if(colors[1] == "DO") {    
      colors = do.colors
    } else if(colors[1] == "HS") {
      colors = hs.colors
    } # else if(colors[1] == "HS")
  } else {
    if(cross == "DO") {    
      colors = do.colors
    } else if(cross == "HS") {
      colors = hs.colors
    } # else if(cross == "HS")
  } # else
  
  num.founders = nrow(colors)
  call = match.call()
  
  # Keep only the founder coefficients from the coef.matrix.
  lod = NULL
  coef = NULL
  if(chr == "X") {
    if(missing(sex)) {
      stop("Sex (either M or F) must be specified on X chromosome.")
    } # if(missing(sex))
    lod  = doqtl$lod$X
    if (is.null(start) | is.null(end)) {
      lod = lod[lod[,2] == chr, ]
    } else {
      lod = lod[lod[,2] == chr & lod[,3] >= start & lod[,3] <= end,]
    }
    coef = doqtl$coef$X
    if(sex == "F") {
      columns = match(paste("F", colors[,1], sep = "."), colnames(coef))
    } else {
      columns = match(paste("M", colors[,1], sep = "."), colnames(coef))
    } # else
    columns = columns[!is.na(columns)]
    coef = coef[,c(1,columns)]
    colnames(coef)[1] = "A"
    colnames(coef) = sub("^[MF]\\.", "", colnames(coef))
    coef = coef[rownames(coef) %in% lod[,1],]
    
    ## Remove outliers
    coef[abs(coef) > remove.outliers | is.na(coef)] = 0
    
    # Center the coefficient values.
    coef[,2:ncol(coef)] = coef[,2:ncol(coef)] + coef[,1]
    coef = coef - rowMeans(coef)
  } else {
    lod = doqtl$lod$A
    if (is.null(start) | is.null(end)) {
      lod = lod[lod[,2] == chr, ]
    } else {
      lod = lod[lod[,2] == chr & lod[,3] >= start & lod[,3] <= end,]
    }
    ## Remove markers with outlier effects
    #coef_subset = doqtl$coef$A[rownames(doqtl$coef$A) %in% lod[,1], 3:9]
    #markers_to_remove = apply(coef_subset, 1, function(x) {any(abs(x) > remove.outliers | is.na(x))})
    #lod = lod[!markers_to_remove,]
    
    intercept = doqtl$coef$A[,1]
    coef = doqtl$coef$A[,(ncol(doqtl$coef$A)-num.founders+1):ncol(doqtl$coef$A)]
    
    ## Remove outliers
    coef[abs(coef) > remove.outliers | is.na(coef)] = 0
    
    coef[,1] = intercept
    colnames(coef)[1] = "A"
    coef = coef[rownames(coef) %in% lod[,1],]
    
    # Center the coefficient values.
    coef[,2:ncol(coef)] = coef[,2:ncol(coef)] + coef[,1]
    coef = coef - rowMeans(coef)
  } # else 
  # Verify that the SNP IDs in the lod & coef matrices match.
  if(!all(lod[,1] == rownames(coef))) {
    stop(paste("The SNP IDs in column 1 of the qtl data frame must match",
               "the SNP IDs in the rownames of the coef matrix."))
  } # if(!all(lod[,1] == rownames(coef)))
  # Verify that the coefficient column names are in the colors.
  if(!all(colnames(coef) %in% colors[,1])) {
    stop(paste("The founder names in the colnames of the coefficient matrix",
               "must be in column 1 of the colors matrix."))
  } # if(!all(colnames(coef) %in% colors[,1]))
  # Convert the chromosome locations to Mb.
  if(max(lod[,3], na.rm = TRUE) > 200) {
    lod[,3] = lod[,3] * 1e-6
  } # if(max(lod[,3], na.rm = TRUE) > 200)
  # Set the layout to plot the coefficients on top and the p-values on the 
  # bottom.
  layout(mat = matrix(1:2, 2, 1), heights = c(0.66666667, 0.3333333))
  par(font = 2, font.lab = 2, font.axis = 2, las = 1, plt =
        c(0.12, 0.99, 0, 0.85), xaxs = "i", lwd = 2)
  # Plot the coefficients.
  plot(lod[,3], coef[,colors[1,1]], type = "l", col = colors[1,3], lwd = 2,
       ylim = c(min(coef, na.rm = TRUE), max(coef * 2, na.rm = TRUE)), xlab = 
         paste("Chr", chr, "(Mb)"), ylab = "Founder Effects", axes = FALSE, ...)
  #abline(v = 0:20 * 10, col = "grey80")
  for(i in 1:nrow(colors)) {
    points(lod[,3], coef[,colors[i,1]], type = "l", col = colors[i,3],
           lwd = 2)
  } # for(i)
  # Draw a legend for the founder names and colors.
  if(legend) {
    legend.side = "topleft"
    if(which.max(lod[,7]) < nrow(lod) / 2) {
      legend.side = "topright"
    } # if(which.max(apply(coef, 1, max)) < nrow(lod) / 2)
    legend(legend.side, colors[,2], col = colors[,3], lty = 1, lwd = 2,
           x.intersp = 0.75, y.intersp = 0.75, bg = "white", cex = 0.8)
  } # if(legend)
  # Add the axis.
  axis(2)
  # Plot a rectangle around the plot.
  par(xpd = NA)
  usr = par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], lwd = 2)
  par(xpd = FALSE)
  # Plot the mapping statistic.
  par(plt = c(0.12, 0.99, 0.35, 1))
  # Single chromosome plot.
  plot(lod[,3], lod[,7], type = "l", lwd = 2, xlab = "",
       ylab = stat.name, ...)
  #abline(v = 0:20 * 10, col = "grey80")
  points(lod[,3], lod[,7], type = "l", lwd = 2)
  # Shade the confidence interval.
  if(conf.int) {
    interval = bayesint(doqtl, chr = chr)
    usr = par("usr")
    rect(interval[1,3], usr[3], interval[3,3], usr[4], col = rgb(0,0,1,0.1), 
         border = NA)
  } # if(!is.na(conf.int))
  mtext(paste("Chr", chr, "(Mb)"), 1, 2)
  usr = par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], lwd = 2)
  par(old.par)
}