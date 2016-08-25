###
##### edgeR
# cd  $curDir &&  module load mugqic_dev/R_Bioconductor && R --no-save --no-restore
# myrc = function(x){x %>% DNAStringSet %>% reverseComplement %>% as.character}
library(edgeR)
library(magrittr)
library(data.table)
library(Biostrings)
library(stringr)
library(dplyr)
options(stringsAsFactors=F)

###
##### edgeR
# cd  $curDir &&  module load mugqic_dev/R_Bioconductor && R --no-save --no-restore
# myrc = function(x){x %>% DNAStringSet %>% reverseComplement %>% as.character}
library(edgeR)
library(magrittr)
library(data.table)
library(Biostrings)
library(stringr)
library(dplyr)
options(stringsAsFactors=F)


# Usage
usage=function(errM) {
    cat("\nUsage : Rscript deseq.R [option] <Value>\n")
    cat("       -o        : list if directories having Process amplicons results files (y_long.Rdata and pdata.Rdata)\n")
    cat("       -m        : Prefix for merged results file\n")
    cat("       -h        : this help\n\n")
    stop(errM)
}
set.seed(123456789)

## default arg values
fastq=NULL
hairpinEnd=NULL
out_path="./"
barcodefile=NULL
hairpinfile=NULL
merge.results=FALSE
merged.prefix=NULL
## get arg variables
ARG = commandArgs(trailingOnly = T)
for (i in 1:length(ARG)) {
    if (ARG[i] == "-o") {
        out_path=ARG[i+1]
    } else if (ARG[i] == "-m") {
        merge.results=TRUE
        merged.prefix=ARG[i+1]
    } else if (ARG[i] == "-h") {
        usage("")
    }
}


amplicons_fn=unlist(strsplit(out_path, ",", fixed=T))


### Sumarize results

summarize_results<-function(amplicons_fn, merged.prefix){
    merged.counts=data.frame()
    for (k in 1:length(amplicons_fn)){
        output.path=amplicons_fn[[k]]
        print(paste0("NOTICE: loading ",file.path(output.path,"y_long.RData")))
        load(file=file.path(output.path,"y_long.RData")) 
        print(paste0("NOTICE: loading ",file.path(output.path,"data.RData")))
        load(file=file.path(output.path,"data.RData"))
        dat = pdata
        dat$"edgeR_count" = y$samples[rownames(dat),"lib.size"]        
        write.csv(dat,file.path(output.path,"samples.csv"),row.names=F)
        if(k==1){
            merged.pdata=dat                        
        }else{
            merged.pdata$"edgeR_count"=merged.pdata$"edgeR_count" + dat[rownames(merged.pdata),"edgeR_count"]            
        }
        suffix=output.path        
        dat = y$count
        dat = cbind(dat,fdata[rownames(dat),])
        write.csv(dat,file.path(output.path,"counts.csv"),row.names=F)        
        merged.counts=rbind(merged.counts, dat)        
    }
    if(!is.null(merged.prefix)){
        write.csv(merged.counts,paste0(merged.prefix,".counts.csv"),row.names=F) 
        write.csv(merged.pdata,paste0(merged.prefix,".samples.csv"),row.names=F) 
    }        
    
}
summarize_results(amplicons_fn, merged.prefix)    