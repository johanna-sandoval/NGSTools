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
    cat("       -r        : read file\n")
    cat("       -b        : barcodes file\n")
    cat("       -g        : hairpin/genes file\n")
    cat("       -s        : barcode starts at\n")
    cat("       -e        : barcode starts at\n")
    cat("       -t        : hairpin/gene starts at\n")
    cat("       -n        : hairpin/gene ends at\n")
    cat("       -o        : output directory\n")
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
## get arg variables
ARG = commandArgs(trailingOnly = T)
for (i in 1:length(ARG)) {
    if (ARG[i] == "-r") {
        fastq=ARG[i+1]
    } else if (ARG[i] == "-b") {
        barcodefile=ARG[i+1]    
    } else if (ARG[i] == "-g") {
        hairpinfile=ARG[i+1]
    } else if (ARG[i] == "-o") {
        out_path=ARG[i+1]
    } else if (ARG[i] == "-s") {
        barcodeStart=as.numeric(ARG[i+1])
    } else if (ARG[i] == "-e") {
        barcodeEnd=as.numeric(ARG[i+1])
    } else if (ARG[i] == "-t") {
        hairpinStart=as.numeric(ARG[i+1])
    } else if (ARG[i] == "-n") {
        hairpinEnd=as.numeric(ARG[i+1])
    } else if (ARG[i] == "-h") {
        usage("")
    }
}


amplicons_fn=list()
amplicons_fn[[hairpinfile]]=list(readfile=fastq, readfile2=NULL, barcodefile=barcodefile,hairpinfile=hairpinfile, barcodeStart=barcodeStart,barcodeEnd=barcodeEnd, hairpinStart=hairpinStart, barcode2Start=NULL, barcode2End=NULL)


for (k in 1:length(amplicons_fn)){
     output.path=file.path(out_path, paste(names(amplicons_fn)[k],"out",sep="_"))
     system(paste0("mkdir -p ", output.path))
     ### Read barcodes & sample annotation
     pdata = read.table(amplicons_fn[[k]]$"barcodefile",header=T,comment.char="", sep="\t",stringsAsFactors=F)
     rownames(pdata) = pdata$"ID"
    ### shRNA library definition        
     fdata = read.delim(amplicons_fn[[k]]$"hairpinfile")
     table(nchar(fdata$Sequences)) 
     fdata %<>% group_by(Sequences) %>% summarise_each(funs(toString)) %>% as.data.frame # this collapses duplicate Sequences
     rownames(fdata) = fdata$"ID"    
     #save data
     save(pdata, fdata,file=file.path(output.path, "data.RData"))    
     #tmp.fname=tempfile(pattern = "file", tmpdir = getwd(), fileext = ".fastq")
     #unlink(tmp.fname)
     #file.symlink(file.path(getwd(), amplicons_fn[[k]]$"readfile"),  tmp.fname)
     tmp.fname=amplicons_fn[[k]]$"readfile"
     ### process
     Sys.time()
     #system("head -n 10 hairpinfile.tsv > test_hairpinfile.tsv && head -n 40000000 raw/reads_R1.fastq > test.fastq")
     #system("head -n 100001 hairpinfile.tsv > test_hairpinfile.tsv")
     y = processAmplicons(   
        readfile = tmp.fname,
        readfile2 = amplicons_fn[[k]]$"readfile2",
        barcodefile = amplicons_fn[[k]]$"barcodefile",
        hairpinfile = amplicons_fn[[k]]$"hairpinfile",
        barcodeStart=amplicons_fn[[k]]$"barcodeStart",
        barcodeEnd=amplicons_fn[[k]]$"barcodeEnd", 
        hairpinStart=amplicons_fn[[k]]$"hairpinStart",
        hairpinEnd = amplicons_fn[[k]]$"hairpinStart" + unique(nchar(fdata$"Sequences")) - 1 ,
        allowShifting=FALSE, 
        shiftingBase=3,
        allowMismatch=TRUE,  
        barcodeMismatchBase=0, 
        hairpinMismatchBase=1,
        allowShiftedMismatch=FALSE, 
        verbose=F,
        barcode2Start = amplicons_fn[[k]]$"barcode2Start",
        barcode2End = amplicons_fn[[k]]$"barcode2End"      
        )
    Sys.time()
    #save(y,file="y.RData")    
    save(y,file=file.path(output.path,"y_long.RData"))
    print (y)
    ## load("y.RData")
}
