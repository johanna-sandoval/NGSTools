################################################################################
# Copyright (C) 2016 Johanna Sandoval
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# If not, see <http://www.gnu.org/licenses/>.
################################################################################



library("BiSeq")
#library(gqSeqUtils)

results_dir="BiSeqDMR"
data_files=list.files(path = results_dir, pattern = "DRM_stats.RData$", all.files = T, full.names = T, recursive = T)        
groups=gsub("/","",gsub("BiSeqDMR", "", gsub("DRM_stats.RData", "", data_files)))

# Just copied the function and eliminated the 1 region validation that sucks
plotMethMapTop <- function(object, region, groups, intervals=FALSE,...) {

  strand(object) <- "*"
  object <- sort(object)
  ind <- overlapsAny(rowRanges(object), region)
  if (sum(ind) == 0) {
    stop("No data to plot within the given region/samples.")
  }

  md <- methLevel(object)[ind, ]
  pos <- start(rowRanges(object))[ind]
  rownames(md) <- as.character(pos)
  md <- t(md)
  ind <- !apply(is.na(md), 2, all)
  md <- md[, ind]
  pos <- pos[ind]

  if (intervals) {
    fullMd <- matrix(NA, nrow=nrow(md), ncol=pos[length(pos)]-pos[1]+1)
    rownames(fullMd) <- rownames(md)
    colnames(fullMd) <- as.character(pos[1]:pos[length(pos)])
    fullMd[, colnames(md)] <- md
    md <- fullMd
    rm(fullMd)
  }
  
  ind.na <- apply(md, 1, function(x) all(is.na(x)))
  if (sum(ind.na) > 0) {
    names.na <- names(which(ind.na))
    md <- md[!ind.na, ]
    object <- object[, !ind.na]
    if (is.element("RowSideColors", names(args))) {
      args$RowSideColors <- args$RowSideColors[!ind.na]
    }
    warning("Sample(s) ", names.na, " omitted due to missing data within \"region\".")
  }

  if (intervals) {
    cn <- as.character(pretty(colnames(md), n=5))
    colnames(md)[!is.element(colnames(md), cn)] <- ""
  }
  
  args <- list(...)
  
  # set RowSideColors if groups is given
  if (!missing(groups) & !is.element("RowSideColors", names(args))) {
    pDat <- as.character(groups)
    uPDat <- unique(pDat)
    cols <- .categorialColors(length(uPDat))
    names(cols) <- uPDat
    args <- c(list(RowSideColors=cols[pDat]), args)
  }

  # default scale="none"
  if (!is.element("scale", names(args))) {
    args <- c(args, list(scale="none"))
  }

  # default zlim=c(0,1)
  if (!is.element("zlim", names(args))) {
    args <- c(args, list(zlim=c(0,1)))
  }
  
  # default Colv = NA
  if (!is.element("Colv", names(args))) {
    args <- c(args, list(Colv=NA))
  }

  # default cexCol=0.8 if intervals=TRUE and no labCol given
  if (intervals &
      !is.element("labCol", names(args)) &
      !is.element("cexCol", names(args))) {
    args <- c(args, list(cexCol=0.8))
  }

  # default colors green - black - red
  if (!is.element("col", names(args))) {
    colF <- colorRampPalette(colors=c("green", "black", "red"))
    args <- c(args, list(col=colF(64)))
  }

  args <- c(list(x=md), args)
  do.call(heatmap, args)
}

# Load Data
# Methylation plots (best DM regions and heatplot with top N diff methylated regions)
plot.top.regions=100
for (fi in 1:length(data_files)){
    output=gsub("RData$","tsv", data_files[fi])
    group=groups[fi]
    load(data_files[fi])
    o <- order(-abs(elementMetadata(DMRs)$median.meth.diff), elementMetadata(DMRs)$median.p)
    sorted.DMRs=DMRs[o,]
    pdf ( file = gsub(".tsv$", ".pdf",output) , width = 10 , height = 5 )
    # COverage boxplots
    covBoxplots(rrbs.clust.lim, las=2, outline=F)    
    title(main = "Sample wise coverage distribution")    
    rowCols <- unlist(lapply(colData(predictedMeth)[group], function(x) {if(inherits(as.numeric(x),"try-error")) {return(NA)} else{return(c("red", "blue")[as.numeric(x)])} })) 
    # top DRM heatmap plots
    plotMethMapTop(predictedMeth,
               region = sorted.DMRs[1:plot.top.regions],
               groups = as.factor(colData(predictedMeth)[,group]),
               intervals = FALSE,
               zlim = c(0,1),
               RowSideColors = rowCols,
               labCol = "", margins = c(0, 6),
               main=paste("Smoothed methylation values within the top ", plot.top.regions, "methylated regions, ","\n","sorted by absolute median meth. difference (decreasing) and median.pvalue(increasing) ")
               )    
    # plots of specific regions, smoothed curves and heatmaps of smoothed methylation values
    for (re in 1:max(1,plot.top.regions)){
        plotSmoothMeth(object.rel = predictedMeth, 
                    region = sorted.DMRs[re] + 1000, 
                    groups = colData(predictedMeth)[,group], 
                    col=c("magenta", "blue"))
        title(main = paste( "Methylation levels for region ", paste(sorted.DMRs[re], collapse=" " )))    
        legend("topright", lty=1, levels(as.factor(colData(predictedMeth)[,group])), col=c("magenta", "blue"))        
        rowCols <- unlist(lapply(colData(predictedMeth)[group], function(x) {if(inherits(as.numeric(x),"try-error")) {return(NA)} else{return(c("red", "blue")[as.numeric(x)])} }))
        try(plotMethMap(predictedMeth,
               region = sorted.DMRs[re],
               groups = as.factor(colData(predictedMeth)[,group]),
               intervals = FALSE,
               zlim = c(0,1),
               RowSideColors = rowCols,
               labCol = "", margins = c(0, 6),
               main=paste("Smoothed methylation values within a detected DMR in region ", paste(sorted.DMRs[re], collapse=" "))
               ))
    }
    dev.off()
    
    # Using our R packages for a different heatmap
#     ind <- overlapsAny(rowRanges(predictedMeth), sorted.DMRs[1:plot.top.regions])
#     md <- methLevel(predictedMeth)[ind, ]
#     pos <- apply(cbind(as.data.frame(seqnames(rowRanges(predictedMeth))[ind]), start(rowRanges(predictedMeth))[ind], end(rowRanges(predictedMeth))[ind]),1, function(x) paste(x, collapse="_"))
#     rownames(md) <- as.character(pos)
#     rc.eset =  md
#     rc.eset = ExpressionSet( as.matrix(rc.eset))     
#     pheatmap.eset(rc.eset, file=gsub(".tsv$",paste0("_top_",plot.top.regions,".pdf"),output))    
}
