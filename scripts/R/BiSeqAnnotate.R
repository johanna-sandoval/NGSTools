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
library("BiSeq")

options(stringsAsFactors=F)

results_dir="BiSeqDMR/"
results_list=list.files(path = results_dir, pattern = "DRM_stats.RData$", all.files = T, full.names = T, recursive = T)        

# return a list containing metadata for ovelapping subjects
splitColumnsByOverlap <-
function(query, subject, columns=c("ENTREZID"), ...)
{
    olaps <- findOverlaps(query, subject, ...)
    f1 <- factor(subjectHits(olaps),
                 levels=seq_len(subjectLength(olaps)))
    splitAsList(mcols(query)[[columns]][queryHits(olaps)], f1)
}        
# Read db to extract genes/transcript db
txdb.UCSC <- TxDb.Hsapiens.UCSC.hg19.knownGene
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
    # Mean of methylation levels
    write.table(apply(as.data.frame(values(annotated.gr)),2,function(x) summary(as.numeric(na.exclude(x)))),gsub(".RData$", ".methLevels.stats.tsv",fname ),quote=F,sep="\t", col.names=T,row.names=T)    
    # subset sites in DM regions 
    annotated.gr=subsetByOverlaps(annotated.gr, DMRs)
    print(annotated.gr)        
    write.table(as.data.frame(annotated.gr),gsub(".RData$", ".methLevels.tsv",fname ),quote=F,sep="\t", col.names=T,row.names=F)
    # Annotate DMRs    
    loc.region.UCSC=locateVariants(DMRs, txdb.UCSC, AllVariants());  
    #loc.region=locateVariants(DMRs, txdb, AllVariants());  
    # They asked to add lines par location, we will give a new line per annotation :-\
    #cols <- c("ENTREZID","ALIAS", "ENSEMBL", "ENSEMBLPROT", "ENSEMBLTRANS", "GENENAME", "GO", "GOALL", "ONTOLOGY", "ONTOLOGYALL" , "PFAM", "UNIPROT")        
    cols <- c("ENTREZID","SYMBOL","ENSEMBL")
    keys <- unlist(na.exclude(unique(loc.region.UCSC$GENEID)));  
    gene.info<-select(org.Hs.eg.db, keys, cols, keytype="ENTREZID");
    gene.splitted<-data.frame(do.call("rbind",lapply(split(gene.info, as.factor(gene.info$ENTREZID)), function(x) return(list(ENTREZID=paste(unique(x[,1]),collapse=","), SYMBOL=paste(unique(x[,2]),collapse=","), ENSEMBL=paste(unique(x[,3]),collapse=","))))))
    # annotate transcripts
    vals <- list("gene_id"=loc.region.UCSC$GENEID, "tx_id"=loc.region.UCSC$TXID)
    cols=c("GENEID", "TXID", "TXNAME")
    tx=transcripts(txdb.UCSC,vals,columns=cols)
    # merge gene annotations + transcripts
    results=merge(as.data.frame(loc.region.UCSC,row.names=NULL),gene.splitted,by.x="GENEID",by.y="ENTREZID", sort=F, all.x=T, all.y=F)
    results=merge(as.data.frame(results),as.data.frame(tx,row.names=NULL),by=c("GENEID", "TXID"), sort=F, all.x=T, all.y=F,suffixes=c("_DMR","_transcript"))
    # Read Regression p-values to count the number of significant sites per DMRs
    betaResults.gr<-with(betaResults, 
                             GRanges(
                                 seqnames=chr,
                                 ranges=IRanges(start=pos, width=1),
                                 p.val=p.val,
                                 cluster.id=cluster.id,
                                 significant=ifelse(p.val<=cutoff.p.val,1,0)))
    #results["nb.significant"]=unlist(lapply(splitColumnsByOverlap(betaResults.gr, loc.region.UCSC, c("significant")), function(x) sum(unlist(x)) ))                             
    #coco=unlist(lapply(splitColumnsByOverlap(betaResults.gr, loc.region.UCSC, c("significant")), function(x) sum(unlist(x)) ))                             
    mcols(loc.region.UCSC)=results
    #Overlap DMRs and annotations
    hits <- findOverlaps(loc.region.UCSC, DMRs)
    values <- as.data.frame(DMRs[unique(subjectHits(hits)),c("median.p","median.meth.group1","median.meth.group2","median.meth.diff")])    
    # merge transcript annotations
    annotated.regions=merge( values, results, sort=F, all=T, by.x=c("seqnames","start","end"),by.y=c("seqnames_DMR","start_DMR","end_DMR"))    
    head(annotated.regions)
    print(nrow(annotated.regions))
    o <- order(-abs(annotated.regions$median.meth.diff), annotated.regions$median.p)    
    annotated.regions=annotated.regions[o,]    
    print(head(annotated.regions))
    # Correct some terms causing problems (\n dans les annotations, integer(0) et character(0) for missing values, seigneur jesus)
    annotated.regions[annotated.regions=="integer(0)"]<-NA
    annotated.regions[annotated.regions=="character(0)"]<-NA    
    annotated.regions[annotated.regions=="<NA>"]<-NA    
    annotated.regions[annotated.regions=="NA"]<-NA    
    # also the ugly \n right in the middle of annotations :-\
    annotated.regions$PRECEDEID=gsub("\\n","", unlist(lapply(annotated.regions$PRECEDEID, function(x) paste(x,collapse=","))))
    annotated.regions$FOLLOWID=gsub("\\n","", unlist(lapply(annotated.regions$FOLLOWID, function(x) paste(x,collapse=","))))
    annotated.regions$SYMBOL=unlist(lapply(annotated.regions$SYMBOL, function(x) paste(x,collapse=",")))
    annotated.regions$ENSEMBL=unlist(lapply(annotated.regions$ENSEMBL,function(x) paste(x,collapse=",")))
    write.table(annotated.regions,gsub(".RData$", ".DMRs_annotated.tsv",fname ),quote=F,sep="\t",col.names=T,row.names=F, fileEncoding = 'UTF-8',na="")
}
