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



suppressPackageStartupMessages(library(gqSeqUtils))
library(magrittr)
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(VariantAnnotation)
library(org.Hs.eg.db)
library("TxDb.Hsapiens.UCSC.hg19.knownGene")

options(stringsAsFactors=F)

results_dir="BiSeqDMR"
targets_file="SeqCap_EPI_4M_CPGiant.targets.bed"
results_list=list.files(path = results_dir, pattern = "DRM_stats.tsv$", all.files = T, full.names = T, recursive = T)
# Read db to extract genes/transcript db

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

cutoff.median.p=0.05


# Select genes or transcripts overlapping genomic regions
itemRanges <- function(db, item="genes",column="ENTREZID")
{
    if(item == "genes"){
        g <- genes(db, columns=column)
    }else{
        if(item == "transcripts"){    
            g <- transcripts(db, columns=column)    
            }        
    }
#     col <- mcols(g)[[column]]
    genes <- granges(g)[rep(seq_along(g), elementLengths(col))]
    mcols(genes)[[column]] <- as.character(unlist(col))
    genes
}
    
# Find which gene coordinates overlap which copy number variant coordinates, and then split the column of gene identifiers into lists corresponding to the regions of overlap
splitByOverlap <-function(query, subject, column="ENTREZID", ...)
{
    best.rank=c("promoter","spliceSite","coding","fiveUTR","threeUTR","intron","intergenic")    
    olaps <- findOverlaps(query, subject, ...)
    f1 <- factor(subjectHits(olaps),
                levels=seq_len(subjectLength(olaps)))
    lapply(splitAsList(mcols(query)[[column]][queryHits(olaps)], f1), function(x) { tmp=unique(unlist(x)); if(length(tmp)>1) { tmp=na.exclude(tmp[match(best.rank,tmp)])}; return(tmp[1]); })
}
        
# Compare and annotate regions        
# Annotate all results
contrast_names=unlist(lapply(results_list, function(x) {pl=unlist(strsplit(x,"/")); return(pl[(length(pl)-1)])}))
gr.list=list()
i=1
for (fname in results_list){
    pl=unlist(strsplit(fname,"/")); 
    cnv.df = read.delim(fname,stringsAsFactors=F)
    print(summary(cnv.df))
    print(paste("Notice: before filtering"))
    cnv.df = subset(cnv.df , median.p <= cutoff.median.p)
    print(paste("Notice: after filtering"))
    print(summary(cnv.df))
    gr.list[i]=makeGRangesFromDataFrame(cnv.df)
    # add metadata (if needed to filter)
    #mcols(gr.list[pl[(length(pl)-1)]])=cnv.df[,3:ncol(cnv.df)]
    rm(cnv.df)    
    i=i+1
}

# Read target regions
unique.gr=unique(c(gr.list[[1]], gr.list[[2]]))
cnv.df = read.delim(targets_file,stringsAsFactors=F, header=F)
head(cnv.df)
cnv.df[,2]=cnv.df[,2]+1
names(cnv.df)=c("seqnames","start","end")
gr.targets=makeGRangesFromDataFrame(cnv.df)

# Compute common and set differences
set.op.gr=list(intersect(gr.list[[1]], gr.list[[2]],type="within"), setdiff(gr.list[[1]], gr.list[[2]]), setdiff(gr.list[[2]], gr.list[[1]]), gr.list[[1]], gr.list[[2]], gr.targets)
names(set.op.gr)=c("shared", paste("set.diff",contrast_names,sep="."), contrast_names, "Target_regions")
#loc.region <- lapply(set.op.gr, function(x) locateVariants(unique(set.op.gr), txdb, AllVariants()))
loc.region <- locateVariants(unique.gr, txdb, AllVariants())
locate.targets=locateVariants(gr.targets, txdb, AllVariants())
genomic.context=lapply(set.op.gr[-length(set.op.gr)], function(x) {table(unlist(splitByOverlap(loc.region, x, "LOCATION", type="equal")), useNA="always")})
genomic.context["Target_regions"]=list(table(unlist(splitByOverlap(locate.targets, gr.targets, "LOCATION", type="equal")), useNA="always"))
genomic.context.counts= Reduce(function(x, y) merge(x, y, all=TRUE, by=1), genomic.context)
names(genomic.context.counts)=c('context', names(set.op.gr))
genomic.context.freq=cbind(context=genomic.context.counts$context, data.frame(apply(genomic.context.counts[,2:ncol(genomic.context.counts)],2, function(x) 100*x/sum(x))))
toplot=melt(genomic.context.freq,id.vars="context")

write.table(merge(genomic.context.counts,genomic.context.freq, by="context", suffixes = c(".counts",".freq"), all=F), "report/GenomicContextofDMRs.tsv", sep="\t", quote=F,row.names=F,col.names=T)    

pdf(file="report/GenomicContextofDMRs.pdf", paper="letter")
ggplot(data=toplot, aes(x=context, y=value, fill=variable)) + 
ggtitle(paste(paste("Genomic context of DMRs", "filtered by median.p <= ", cutoff.median.p, sep=" "), paste(names(set.op.gr),paste("=",unlist(lapply(set.op.gr, length))), collapse = ",\n"), sep="\n")) +
ylab("Percent of DMRs") +
geom_bar(stat="identity", position=position_dodge()) 

dev.off()  
