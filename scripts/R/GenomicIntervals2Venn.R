vcfFiles=NULL
gtfFile=NULL
# vcfFiles = paste(list.files("peak_call", pattern="peaks.xls$", full.names=T,recursive=T),collapse=",")
vcfFiles = "peak_call/Basal/Basal_peaks.xls,peak_call/IL6/IL6_peaks.xls,peak_call/IL6_R1881/IL6_R1881_peaks.xls,peak_call/R1881/R1881_peaks.xls"
type="macs2"
gtfFile="../genome/nDi.2.2.2.gtf"
outFile="report/UniqueGenomicIntervalsAnnotated"
seq.alias="human"
comparisons=eval(parse(text="list('all_groups'=c('Basal_peaks','R1881_peaks','IL6_peaks','IL6_R1881_peaks'))"))
overlap_type ="any"
reduce_unique=TRUE

usage=function(errM) {
    cat("\nUsage : Rscript GenomicIntervals2Venn.R [option] \n")
    cat("       -v        : GenomicIntervals file(s)\n")
    cat("       -t        : Genomic Intervals type (bed=bed file, macs2=standard macs2 peaks xls file, vcf=standard vcf file, other=other tab separated file)\n")
    cat("       -g        : gtf file \n")
    cat("       -o        : output prefix\n")
    cat("       -c        : comparisons\n")
    cat("       -y        : find overlaps types, (any=any overlap , start, end, within, equal\n")
    cat("       -r        : reduce list of unique genomic ranges, default=False\n")
    cat("       -h        : this help\n\n")
    
    stop(errM)
}

#module load mugqic_dev/R_Bioconductor && Rscript bin/GenomicIntervals2Venn.R -v "peak_call/Basal/Basal_peaks.xls,peak_call/IL6/IL6_peaks.xls,peak_call/IL6_R1881/IL6_R1881_peaks.xls,peak_call/R1881/R1881_peaks.xls" -t macs2 -o comparisons -c "list('all_groups'=c('Basal_peaks','R1881_peaks','IL6_peaks','IL6_R1881_peaks'))"
 
ARG = commandArgs(trailingOnly = T)
## get arg variables
for (i in 1:length(ARG)) {
if (ARG[i] == "-v") {
 vcfFiles=ARG[i+1]
} else if (ARG[i] == "-g") {
 gtfFile=ARG[i+1]
} else if (ARG[i] == "-t") {
  type = ARG[i+1]
} else if (ARG[i] == "-o") {
 outFile=ARG[i+1]
} else if (ARG[i] == "-c") {
comparisons=eval(parse(text=ARG[i+1]))            
} else if (ARG[i] == "-y") {
    overlap_type = ARG[i+1]
} else if (ARG[i] == "-r") {
    reduce_unique=TRUE        
} else if (ARG[i] == "-h") {        
 usage("")
}

}


# Read Genomic Variants File
readGVFile<-function(file, type="bed" , unique.id="gnames"){    
    if (type=="bed"){
        reads <- import.bed(con=file, asRangedData=F)
    }else if (type=="macs2"){        
        data=read.delim(file,sep="\t",comment.char = "#",stringsAsFactors=F)
        data=cbind(data,strand="+")
        colnames(data)[1:3] <- c('chr','start','end')
        reads <- with(data, GRanges(chr, IRanges(start, end), strand, score=fold_enrichment, id=name))
        mcols(reads)=data[,-c(1:3)]
        # unique ID for macs peaks is macs2 peaks name
        names(reads)=apply(data[,1:3], 1, function(x) gsub(" ","",paste(x,collapse="_",sep="_") ))
        mcols(reads)[unique.id]=names(reads)
    }else if (type=="vcf"){
        reads=readVcf(fi,seq.alias)
        mcols(reads)[unique.id]=geno(reads)$GT
    }else{
        data=read.delim(file,sep="\t")        
        colnames(data)[1:3] <- c('chr','start','end')
        reads <- with(data, GRanges(chr, IRanges(start, end)))
        names(reads)=apply(data[,1:3], 1, function(x) gsub(" ","",paste(x,collapse="_",sep="_") ))
    }    
    return(reads)
        
}

# Read variants or genomic Ranges files
## Extract the ranges from the VCF/bed/other file 
i=1
common.gr=NULL
merged.gr=NULL
sample.names=c()
# fi=unlist(strsplit(vcfFiles,","))[1]
# Forget grangesList, a simple list did the job
grList=list()
for (fi in unlist(strsplit(vcfFiles,","))){
    sname=gsub("xls|scalpel|filtered|vcf|bed|bedGraph|[.]", "", basename(fi),fixed=F,perl=T)
    sample.names=c(sample.names,sname)
    #hdr <- scanVcfHeader(file)
    #vcf <- readVcf(file)
    if ( i == 1 ){      
        if(type=="vcf"){
            vcf0=gr0=readGVFile(fi, type=type, unique.id=sname)
            gr0<- rowRanges(vcf0)
        }else{
            vcf0=gr0=readGVFile(fi, type=type, unique.id=sname)
        }
        unique.gr=gr0        
        merged.gr<-gr0
        grList[sname]<-gr0        
    }else{
        if(type=="vcf"){
            vcf=readGVFile(fi, type=type, unique.id=sname)
            gr<- rowRanges(vcf)
        }else{
            vcf=gr=readGVFile(fi, type=type, unique.id=sname)                
        }
        grList[sname]<-gr        
        # Keep only common columns in intersection        
        common.names=setdiff(intersect(names(mcols(unique.gr)),names(mcols(gr))),c("seqnames", "ranges", "strand", "seqlevels", "seqlengths", "isCircular", "start", "end", "width", "element"))
        mcols(gr)=mcols(gr)[,common.names]
        mcols(unique.gr)=mcols(unique.gr)[,common.names]
        print(paste("NOTICE: common columns of input files are restricted to : ", common.names))
        unique.gr=unique(c(unique.gr,gr))        
        
    }
    if (i==2){
        common.gr=intersect(gr0,gr)
        merged.gr<-union(gr0,gr)
    }else if (i>2){
        common.gr=intersect(common.gr,gr)
        merged.gr<-union(merged.gr,gr)
    }    
    i=i+1    
}
if(is.null(common.gr) && !is.null(gr0)){
    common.gr=gr0
}
names( grList)<-sample.names
print("NOTICE: merged.gr:")
print(merged.gr)
print("NOTICE: unique.gr:")
print(unique.gr)
print("NOTICE: common.gr:")
print(common.gr)
if(reduce_unique){
    unique.gr=reduce(unique.gr)
    names(unique.gr)=apply(data.frame(seqnames=seqnames(unique.gr), start=start(unique.gr), end=end(unique.gr)), 1, function(x) gsub(" ","",paste0(x,collapse="_")))
    print("NOTICE: Reduced unique.gr:")
    print(unique.gr)    
}

## Venn diagrams
# comparisons=list(
#    TRS_vs_resistant=c("TRS_vs_JYD_survivor","TRS_vs_Tootie-RES","TRS_vs_Jojo-RES"),
#    Susceptible_vs_resistant=c("susceptible_vs_JYD_survivor","susceptible_vs_Tootie-RES","susceptible_vs_Jojo-RES"),
#    All_samples=c("TRS_vs_JYD_survivor","TRS_vs_Tootie-RES","TRS_vs_Jojo-RES","susceptible_vs_JYD_survivor","susceptible_vs_Tootie-RES","susceptible_vs_Jojo-RES")
# )

# splitColumnByOverlap <-
# function(query, subject, column="gnames", ...)
# {
#     olaps <- findOverlaps(query, subject, ...)
#     f1 <- factor(subjectHits(olaps),
#                  levels=seq_len(subjectLength(olaps)))
#     splitAsList(mcols(query)[[column]][queryHits(olaps)], f1)
# }
# # to get genotype of any 
# subsetByOverlapsAnnotateQuery<-function(query, subject, select.cols=c("val"), ...){
#     hits <- cbind(findOverlaps(query, subject, ignore.strand=T, select="first", ...),c(1:length(query)))    
#     # select one hit per query
#     subject.idx<-unlist(lapply(1:length(subject), function(x) if(x %in% hits[,1]) {return(subset(hits, hits[,1]==x,2)[1,])} else{return(NA)}))
#     values<-unlist(lapply(subject.idx, function (x) {if (is.na(x)) {return(rep(NA,length(select.cols)))} else {return (mcols(query)[x,select.cols])}}))        
#     return(values)
# }

## Get unique identifiers (names(unique.gr) of intervals overlapping)
splitColumnsByOverlap <-
function(query, subject, column="gnames", ...)
{
    olaps <- findOverlaps(query, subject, ...)
    f1 <- factor(subjectHits(olaps),
                 levels=seq_len(subjectLength(olaps)))
    
    splitAsList(names(subject)[subjectHits(olaps)], f1)
}
    
if(overlap_type =="equal"){       
    elementsList=lapply(grList[comparisons[[l]]], names)
}else{
    elements.df=NULL
    for (i in 1:length(grList)){
        if(is.null(elements.df)){ 
            elements.df=data.frame(dummy=unlist(lapply(splitColumnsByOverlap(query=grList[[i]], subject=unique.gr, column=names(grList)[i], type=overlap_type), function(x) paste(sort(unlist(x)),collapse=""))))
        }else{
            elements.df[names(grList)[i]]=unlist(lapply(splitColumnsByOverlap(query=grList[[i]], subject=unique.gr, column=names(grList)[i], type=overlap_type), function(x) paste(sort(unlist(x)),collapse="")))
        }
    }
    colnames(elements.df)=names(grList)
    elementsList=as.list(as.data.frame(elements.df))
    # remove blanks 
    elementsList=lapply(elementsList, function(x) x[which(!is.na(x))])
    names(elementsList)=names(grList)
    #elementsList=lapply(grList, function(x) unlist(splitColumnByOverlap(query=x, subject=unique.gr, column=names(grList)[i], type=overlap_type)))    
}    
    
    

print("Notice: comparing the following groups: ")
print(comparisons)

for (l in names(comparisons)){
    v=Venn(elementsList)
    dnames=attr(attr(v,"IntersectionSets"), "names")    
    header=cbind("Variant",t(dimnames(attr(v,"IndicatorWeight"))[[2]][1:dim(attr(v,"IndicatorWeight"))[[2]]-1]))
    write.table(header,paste0(outFile,l,".tsv"),row.names=F,col.names=F,sep="\t",quote=F,append=F)
    count=0
    for (i in attr(v,"IntersectionSets")){
            count=count+1
            line=cbind(i,paste(unlist(strsplit(dnames[count],"")),collapse="\t"))
            write.table(line,paste0(outFile,l,".tsv"),row.names=F,col.names=F,sep="\t",quote=F,append=T)
    }   
    pdf(paste0(outFile,l,".pdf"),paper="letter")  
    #vennplot(grList,"gplots")    
    C2 <- compute.Venn(v, doWeights = FALSE)
    grid.newpage()
    plot(C2)
    dev.off()
                    
}

