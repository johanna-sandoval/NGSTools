suppressPackageStartupMessages(library(gqSeqUtils))
library(magrittr)
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)

options(stringsAsFactors=F)

results_dir="freec"
results_list=list.files(path = results_dir, pattern = "_CNVs.p.value.txt$", all.files = T, full.names = T, recursive = T)

# select genes or transcripts overlapping genomic regions
itemRanges <- function(db, item="genes",column="ENTREZID")
{
    if(item == "genes"){
        g <- genes(db, columns=column)
    }else{
        if(item == "transcripts"){    
            g <- transcripts(db, columns=column)    
        }        
    }
    col <- mcols(g)[[column]]
    genes <- granges(g)[rep(seq_along(g), elementLengths(col))]
    mcols(genes)[[column]] <- as.character(unlist(col))
    genes
}

# Find which gene coordinates overlap which copy number variant coordinates, and then split the column of gene identifiers into lists corresponding to the regions of overlap
splitByOverlap <-function(query, subject, column="ENTREZID", ...)
{
    olaps <- findOverlaps(query, subject, ...)
    f1 <- factor(subjectHits(olaps),
                 levels=seq_len(subjectLength(olaps)))
    splitAsList(mcols(query)[[column]][queryHits(olaps)], f1)
}

# Read db to extract genes/transcript db
txdb=makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",dataset="celegans_gene_ensembl", transcript_ids=NULL, circ_seqs=DEFAULT_CIRC_SEQS, filters="",id_prefix="ensembl_",
                         host="oct2014.archive.ensembl.org",port=80, miRBaseBuild=NA) 
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="celegans_gene_ensembl", host="oct2014.archive.ensembl.org",port=80,
                archive=FALSE)   
mart<-useDataset("celegans_gene_ensembl",mart)

# Annotate all cnv results
for (fname in results_list){
    cnv.df = read.delim(fname) 
    cnv = makeGRangesFromDataFrame(cnv.df)        
    gns = itemRanges(txdb, item="genes", column="gene_id")
    cnv.df = cbind(cnv.df, gene_overlap=unstrsplit(splitByOverlap(gns, cnv, "gene_id"), ": ") )
    tx =  itemRanges (txdb, item="transcripts", column="tx_name")
    cnv.df = cbind(cnv.df, transcript_overlap=unstrsplit(splitByOverlap(tx , cnv, "tx_name"), ": "))
    output_name = gsub(".txt",".annotated.tsv", fname)
    write.table(cnv.df, output_name, sep="\t", quote=F, row.names=F)
    output_genes_name = gsub(".txt",".annotated_genes.tsv", fname)
    results<-getBM(attributes=c("ensembl_gene_id", "chromosome_name", "start_position","end_position","ensembl_transcript_id", "description" ,"go_id", "name_1006", "definition_1006"),filters="chromosomal_region",values=as.data.frame(cnv), mart=mart) 
    write.table(results, output_genes_name, sep="\t", quote=F, row.names=F)
}
  
