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

library(magrittr)
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library("BiSeq")
library(ggplot2)

options(stringsAsFactors=F)

results_dir="BiSeqDMR/"
results_list=list.files(path = results_dir, pattern = "DRM_stats.RData$", all.files = T, full.names = T, recursive = T)        
# Compare and annotate regions
cutoff.p.val=0.05
# Annotate all results
contrast_names=unlist(lapply(results_list, function(x) {pl=unlist(strsplit(x,"/")); return(pl[(length(pl)-1)])}))
gr.list=list()
for (fname in results_list){
    pl=unlist(strsplit(fname,"/"));     
    load(fname)
    # Convert raw object to rel
    rrbs<-rawToRel(rrbs.clust.lim)
    # raw  and smoothed methylation levels
    annotated.gr=intersect(rowRanges(rrbs),rowRanges(predictedMeth))
    tmp=data.frame(methLevel(rrbs)); colnames(tmp)=paste(names(tmp),"raw",sep=".")
    values(annotated.gr)=cbind(values(annotated.gr),tmp)
    tmp=data.frame(methLevel(predictedMeth)); colnames(tmp)=paste(names(tmp),"smooth",sep=".")
    values(annotated.gr)=cbind(values(annotated.gr),tmp)
    print(annotated.gr)
    annotated.df=as.data.frame(values(annotated.gr))
    # Mean of methylation levels
    p=list()    
    for (item in  c("raw","smooth") )  {
        item.values=grep(item,names(values(annotated.gr)))
        toplot<-melt(annotated.df[,item.values])
#         p[item]=list(ggplot(toplot) + geom_density(aes(x = value, colour = variable)) +
#                      theme(legend.position = "right") + labs(title = paste0("Histogram of ",item, " methylation values")))
        p[paste0(item,".q")]=list(ggplot(toplot, aes(value, colour = variable)) +
                                  geom_freqpoly(binwidth = 0.05) +  xlim(0.75, 1) +
                                  theme(legend.position = "right") + labs(title = paste0("Histogram of ",item, " methylation levels"))
                                  )
    }    
    pdf(gsub(".RData$", ".methLevels.stats.pdf",fname ))
    invisible(lapply(p, print))
    dev.off()
    #write.table(apply(as.data.frame(values(annotated.gr)),2,function(x) summary(as.numeric(na.exclude(x)))),gsub(".RData$", ".methLevels.stats.tsv",fname ),quote=F,sep="\t", col.names=T,row.names=T)    
}
