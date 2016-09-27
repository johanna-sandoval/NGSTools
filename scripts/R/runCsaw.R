options(digits=3)


## This code was originally taken from bioconductor's csaw lab
## at https://www.bioconductor.org/help/course-materials/2015/BioC2015/csaw_lab.html

suppressPackageStartupMessages({
  library(csaw)
  library(edgeR)
  library(Gviz)
  library(GenomicFeatures)
})


usage=function(errM) {
    cat("\nUsage : Rscript runCsaw.R -d design.dbr.tsv.csv [option] <Value>\n")
    cat("       -d        : design file\n")
    cat("       -f        : fragment length\n")
    cat("       -w        : window width\n")
    cat("       -s        : spacing\n\n")
    cat("       -q        : minimum mapping quality\n\n")
    cat("       -a        : annotation databases generated from UCSC as TxDb objects (default is TxDb.Hsapiens.UCSC.hg19.knownGene)\n\n")
    cat("       -o        : gene based information for the organism (default is org.Hs.eg.db) \n\n")
    cat("       -t        : type of sequencing (none=single end; both=paired end; first=paired end, use only first read; second=paired end, use only second read.)\n\n")
    cat("       -h        : this help\n\n")    
    stop(errM)
}

## default arg values

bam.files<- list.files("alignment", "sorted.filtered.bam$", full.names=T, recursive=T)
frag.len <- 250
max.frag <- 500
window.width <- 10
spacing <- 50
minq <- 20
OrgDb <- "org.Hs.eg.db"
TxDb <- "TxDb.Hsapiens.UCSC.hg19.knownGene"
BSgenome <- "BSgenome.Hsapiens.UCSC.hg19"
design_file<-file.path("./","design.dbr.tsv.csv")
output_path="differential_binding"
ARG = commandArgs(trailingOnly = T)
bin_width=10000
pe="none"

# By default, only proper pairs are used in which the two paired reads
# are on the same chromosome, face inward and are no more than max.frag apart max.frag=400
ARG = commandArgs(trailingOnly = T)

## get arg variables
for (i in 1:length(ARG)) {
    if (ARG[i] == "-d") {
        design_file=ARG[i+1]
    } else if (ARG[i] == "-f") {
        frag.len=as.numeric(ARG[i+1])
    } else if (ARG[i] == "-o") {        
        output_path=ARG[i+1]
    } else if (ARG[i] == "-g") {
        grouping=eval(parse(text=ARG[i+1]))      
    } else if (ARG[i] == "-q") {
        minq=as.numeric(ARG[i+1])        
    } else if (ARG[i] == "-w") {
        window.width =as.numeric(ARG[i+1])
    } else if (ARG[i] == "-s") {
        spacing=ARG[i+1]
    } else if (ARG[i] == "-a") {
        TxDb <- ARG[i+1]        
    } else if (ARG[i] == "-t") {
        pe <- ARG[i+1]        
    } else if (ARG[i] == "-o") {
        OrgDb=ARG[i+1]
        
    } else if (ARG[i] == "-h") {
        usage("")
    }    
}

stopifnot(file.exists(design_file))

# Load annotation libraries
eval(parse(text=paste("library(", TxDb,")")))
eval(parse(text=paste("library(", OrgDb,")")))
eval(parse(text=paste("TxDb=", TxDb)))
eval(parse(text=paste("OrgDb=", OrgDb)))


perform_dbr<-function(bam.files,              
              pe,        
              grouping,        
              frag.len, 
              window.width , 
              spacing, 
              minq, 
              OrgDb, 
              TxDb, 
              BSgenome,
              out_path,              
              contrast_name
){
    ## ------------------------------------------------------------------------
    system.file("doc", "sra2bam.sh", package="csaw")

    ## ------------------------------------------------------------------------
    ## Paired-end data can also be treated as single-end by specifying pe="first" or     "second" in the     readParam    constructor.
    ## ------------------------------------------------------------------------    
    param <- readParam(minq=minq, dedup=TRUE, pe=pe)
    data <- windowCounts(bam.files, ext=frag.len, width=window.width, spacing=spacing, param=param)
    data
    ## ------------------------------------------------------------------------    
    ##   Cross-correlation plots can be generated directly from BAM les using the correlateReads
    ##   function. This provides a measure of the immunoprecipitation (IP) efficiency of a ChIP-seq
    ##   experiment [Kharchenko et al., 2008]. Efficient IP should yield a smooth peak at a delay
    ##   distance corresponding to the average fragment length
    ##
    
    max.delay <- 2 * frag.len
    dedup.on <- readParam(dedup=TRUE, minq=minq)
    
    pdf(file.path(out_path, "differential_binding_results.pdf"), width=11, height=8.5 )
    
    ## plot cross correlation
    x <- correlateReads(bam.files, max.delay, param=dedup.on)
    plot(0:max.delay, x, type="l", ylab="CCF", xlab="Delay (bp)", main=paste0("cross-correlation plot, design",contrast_name))
    
    ## QC fragment size
    if(pe == "both"){
        for (pe.bam in bam.files){
            out <- getPESizes(pe.bam)
            frag.sizes <- out$sizes[out$sizes<=max.frag]
            hist(frag.sizes, breaks=50, xlab="Fragment sizes (bp)", ylab="Frequency", main=paste0("Sample: ",pe.bam))
            abline(v=max.frag, col="red")
            print(paste0("NOTICE: fragment sizes too large for samples bam: ",pe.bam))
            print(c(out$diagnostics, too.large=sum(out$sizes > max.frag)))
        }
    }

    ## ------------------------------------------------------------------------
    ## Load data, count on bins 
    binned <- windowCounts(bam.files, bin=TRUE, width=bin_width, param=param)

    ## ------------------------------------------------------------------------
    ## Normalize bins
    normfacs <- normOffsets(binned)
    print(paste0("NOTICE: Normalization factors: ",normfacs))    

    ## ------------------------------------------------------------------------
    ## MA plots
    adj.counts <- cpm(asDGEList(binned), log=TRUE)
    for (i in 1:(length(bam.files)-1)) {
        cur.x <- adj.counts[,1]
        cur.y <- adj.counts[,1+i]
        smoothScatter(x=(cur.x+cur.y)/2+6*log2(10), y=cur.x-cur.y,
                      xlab="A", ylab="M", main=paste(bam.files[1], "vs", bam.files[i+1]))
        all.dist <- diff(log2(normfacs[c(i+1, 1)]))
        abline(h=all.dist, col="red")
    }
    
    ## ------------------------------------------------------------------------
    ## Filters 
    ## Filtering away low-abundance windows
  
    abundance <- aveLogCPM(asDGEList(data))
    keep <- abundance > -1
    original <- data
    data <- data[keep,]
    print(paste0("Notice: after filtering, still ",nrow(data), " points."))
    
    ## ------------------------------------------------------------------------
    binned.2 <- windowCounts(bam.files, bin=TRUE, width=bin_width/5, param=param)
    filter.stat.bg <- filterWindows(original, background=binned.2, type="global")
    background.keep <- filter.stat.bg$filter > log2(3)
    print(paste0("NOTICE: Testing a different threshold for abundance "))
    summary(background.keep)

    ## ------------------------------------------------------------------------
    ## We can mimic peak callers like MACS, where the background is estimated locally for each window. 
    ## In this example, we define the local background as the 2 kbp interval around each window.
    surrounds <- 2000
    neighbour <- suppressWarnings(resize(rowRanges(original), surrounds, fix="center"))
    wider <- regionCounts(bam.files, regions=neighbour, ext=frag.len, param=param)

    ## ------------------------------------------------------------------------
    filter.stat.loc <- filterWindows(original, wider, type="local")
    local.keep <- filter.stat.loc$filter > log2(3)
    summary(local.keep)

    ## ------------------------------------------------------------------------
    ## Testing for differential binding
    
    design <- model.matrix(~0 + grouping)
    colnames(design) <- paste(contrast_name, levels(as.factor(grouping)), sep="_")    
    y <- asDGEList(data, norm.factors=normfacs)
    y <- estimateDisp(y, design)

    ## The QL dispersions account for window-specific variability, while the NB dispersions model biological variability between replicates.
    fit <- glmQLFit(y, design, robust=TRUE)

    ## ------------------------------------------------------------------------
    contrast <- makeContrasts(paste(colnames(design),collapse=" - "), levels=design)
    results <- glmQLFTest(fit, contrast=contrast)

    ## ------------------------------------------------------------------------
    ## TODO If no replicates a model with fixed dispersion should be used   
#     norep.fit <- glmFit(y, design, dispersion=0.05)
#     norep.results <- glmLRT(norep.fit, contrast=contrast)

    ## ------------------------------------------------------------------------
    bin.adjc <- cpm(asDGEList(binned.2), log=TRUE)
    plotMDS(bin.adjc, labels=grouping)

    ## ------------------------------------------------------------------------
    plotBCV(y)

    ## ------------------------------------------------------------------------
    plotQLDisp(fit)

    ## ------------------------------------------------------------------------
    clustered <- mergeWindows(rowRanges(data), tol=1000)
    clustered$region

    ## ------------------------------------------------------------------------    
    ## Correcting per multiple combineTests
    
    tabcom <- combineTests(clustered$id, results$table)
    head(tabcom)
    
    ## ------------------------------------------------------------------------
    tab.best <- getBestTest(clustered$id, results$table)
    print("Notice: differential binding, top results")
    head(tab.best)

    write.table(tab.best, file.path(out_path,"BestTest.tsv"), sep="\t", row.names=T,col.names=T,quote=F)

    ####
    ## Annotate
    ## ------------------------------------------------------------------------
    gene.bodies <- genes(TxDb)
    prom <- promoters(gene.bodies, upstream=3000, downstream=1000)
    head(prom)

    ## ------------------------------------------------------------------------
    olap <- findOverlaps(prom, rowRanges(data))
    print("NOTICE: Overlap with gene model")
    print(olap)

    ## ------------------------------------------------------------------------
    tabbroad <- combineOverlaps(olap, results$table)
    head(tabbroad[!is.na(tabbroad$PValue),])

    ## ------------------------------------------------------------------------
    anno <- detailRanges(clustered$region, orgdb=OrgDb,
                        txdb=TxDb)
    head(anno$overlap)

    ## ------------------------------------------------------------------------
    head(anno$left)

    ## ------------------------------------------------------------------------
    ## Merge annotations, adjusted  p-values and general results
    
    all.results <- data.frame(as.data.frame(clustered$region)[,1:3], tabcom, anno)
    all.results <- all.results[order(all.results$PValue),]
    write.table(all.results, file=file.path(out_path,"differential_binding_results.tsv"), row.names=FALSE, quote=FALSE, sep="\t")

    ## ------------------------------------------------------------------------
    all.regions <- clustered$region
    elementMetadata(all.regions) <- tabcom
    print(all.regions)
    dev.off()
}

# Read design, process bam files
design = read.csv2(design_file, header=T, sep = "\t", na.strings = "0", check.names=F,colClasses =c('character',rep('numeric',unique(count.fields(design_file))-1)))
print(design)

name_sample= as.character(as.vector(design[,1]))

# Iterate over each design
for (i in 2:ncol(design)) {
    
    name_folder = paste(output_path,names(design[i]),sep="/")     
    # Create output directory              
    if (!file.exists(name_folder)) { 
        dir.create(name_folder, showWarnings=F, recursive=T)
    }
    current_design=design[,i]
    subsampleN=name_sample[!(is.na(current_design))]
    group = as.character(current_design)[!(is.na(current_design))]
    groupN = unique(group)
    
    cat("Processing design\n")
    cat(paste("Name folder: ",name_folder,"\n",sep=""))
    cat(paste("Design : ",paste(subsampleN, group,sep="=",collapse=" ; "),"\n",spe=""))
    samples_bam_files=c()
    for (s in subsampleN){
        samples_bam_files=c(bam.files=samples_bam_files, list.files(file.path("alignment",s), "sorted.dup.bam$", full.names=T, recursive=T))
    }    
    names(samples_bam_files)=subsampleN
#     bam.files=samples_bam_files;   
#     grouping=group;        
#     pe="none";
#     frag.len=frag.len; 
#     window.width=window.width; 
#     spacing=spacing; 
#     minq=minq; 
#     OrgDb=OrgDb; 
#     TxDb=TxDb; 
#     BSgenome=BSgenome;
#     out_path=file.path(output_path, name_folder)
    try(perform_dbr(bam.files=samples_bam_files,
                grouping=group,       
                pe=pe,
                frag.len=frag.len, 
                window.width=window.width, 
                spacing=spacing, 
                minq=minq, 
                OrgDb=OrgDb, 
                TxDb=TxDb, 
                BSgenome=BSgenome,
                out_path=name_folder,
                contrast_name=names(design[i])
    ))
    # perform_dbr(samples_bam_files, group, frag.len, window.width , spacing, minq, OrgDb, TxDb, BSgenome, output_path=name_folder)
}
