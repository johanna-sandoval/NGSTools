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
library(gqSeqUtils)

results_dir="BiSeqDMR"
data_files=list.files(path = results_dir, pattern = "DRM_stats.RData$", all.files = T, full.names = T, recursive = T)        
groups=gsub("/","",gsub("BiSeqDMR", "", gsub("DRM_stats.RData", "", data_files)))

exploratoryAnalysisBiSeq <-function(
    predictedMeth,
    output.dir     = "exploratory",
    reset.index = TRUE,
    df=NULL, 
    ...
)
{
    ## Create output dir
    dir.create(output.dir,showWarnings = FALSE)    
    md <- round(methLevel(predictedMeth),2)
    rownames(md) <- with(rowRanges(predictedMeth),paste(seqnames,start, sep="_"))        
    print(summary(md))    
    # Logit transformation 
    md [md == 0] <- 1*10^-10
    md [md == 1] <- 1-1*10^-10
    meth.odds = log(md / (1-md))
    print(summary(meth.odds ))
    rc.eset = ExpressionSet( na.exclude(as.matrix(meth.odds)))
    print(rc.eset)
    l = list(        
        # Boxplot of converted raw values (no variablity columns eliminated)
        uA.boxplot(rc.eset, path = output.dir, fn = "boxplot_gene_raw_methLevel.pdf", is.log2 = T, 
                   alt.colnames = pData(rc.eset)[['.label']], ylab = "Methylation odds = log(pred.meth / (1-pred.meth) )", xlab = "Sample", 
                   desc = "Boxplot visualization of methylation levels", 
                   ...)
        
        # Standard plots: Correlation distance clustering 
        ,uA.cordist.hclust(rc.eset, path = output.dir,fn = "cordist_hclust_methLevel.pdf",
                           alt.colnames = pData(rc.eset)[[".label"]], is.log2=T, xlab='Samples',
                           desc = "Hierachical clustering based on the correlation distance, methylation odds = log(pred.meth / (1-pred.meth) )", 
                           method='ward.D', 
                           ...)
        ,uA.pca(rc.eset, path = output.dir, fn = "pca_methLevel.pdf", is.log2 = T, 
                alt.colnames = pData(rc.eset)[['.label']], pdata=df,
                desc = "PCA (first two components) of the methylation odds = log(pred.meth / (1-pred.meth) )", 
                ...)
        ,uA.mds(rc.eset, path = output.dir, fn = "mds_methLevel.pdf", is.log2 = T, 
                alt.colnames = pData(rc.eset)[['.label']], pdata=df,
                desc = "Multidimensional Scaling (2D) of the methylation odds = log(pred.meth / (1-pred.meth) )", 
                ...)
        ,uA.pairwise.manhattan.hclust(rc.eset, path = output.dir, fn = "pw_mean_abs_diff_methLevel.pdf", is.log2 = T, 
                                      alt.colnames = pData(rc.eset)[['.label']],
                                      lab = "Samples",desc = "Clustering of the mean absolute difference distance for methylation odds = log(pred.meth / (1-pred.meth) ) "
                                      , ...)         
        ) 
        # Most variable sites by Standard devitation of Methylation level
        fn = file.path(output.dir, "top_sd_heatmap_methLevel.pdf")
        pheatmap.eset(rc.eset, file = fn, scale = "row", col.labels = pData(rc.eset)[[".label"]], cols.ann = colnames(df))
        l = c(l, list(data.frame(Description = "Heat map of most varying sites by methylation odds = log(pred.meth / (1-pred.meth) ) standard deviation", 
                                 path = fn)))
}

# Load Data
for (fi in 1:length(data_files)){
    output=gsub("RData$","tsv", data_files[fi])
    group=groups[fi]
    load(data_files[fi])
    try(exploratoryAnalysisBiSeq(rawToRel(rrbs.clust.lim), output.dir= file.path("exploratory_raw",group), main="Raw methylation levels") )
    try(exploratoryAnalysisBiSeq(predictedMeth, output.dir = file.path("exploratory_smoothed",group), main="Smoothed methylation levels") )
    
}
