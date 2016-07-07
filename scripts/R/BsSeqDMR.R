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

# using the BSseq library to test for differentially methylated regions
#source("http://bioconductor.org/biocLite.R")
#biocLite("bsseq")
#biocLite("bsseqData")

usage=function(errM) {
     cat("\nUsage : Rscript BsSeqDMR.R [option] \n")
     cat("       -i        : input file(s), a comma separated string for path to input files (methylated profiles containting chromosome, position, coverage and methyl counts, tab separated, header=T) \n")
     cat("       -s        : sampleNames, a comma separated string (used if multiple files and each file correspond to a specific sample) \n")
     cat("       -o        : output prefix\n")
     cat("       -c        : coverage counts columns in input files\n")
     cat("       -m        : methylated counts in input files\n\n")
     cat("       -d        : path to the design file\n")
     cat("       -g        : name of the group that must be used as control\n")
     cat("       -v        : variance estimation for BSmooth.tstat (same=variance in the two groups are the same, paired=for paired samples, group2=only estimate the variance based on group 2) \n")
     cat("       -h        : this help\n\n")     
     stop(errM)
}
variance_estimation=NULL 
ARG = commandArgs(trailingOnly = T)
 ## get arg variables
for (i in 1:length(ARG)) {
 if (ARG[i] == "-i") {
 input_files=ARG[i+1]
 } else if (ARG[i] == "-s") {
 sampleNames=ARG[i+1]
 } else if (ARG[i] == "-o") {
 outFile=ARG[i+1]
 } else if (ARG[i] == "-c") {
   coverage.counts.cols=as.numeric(split(ARG[i+1],","))
 } else if (ARG[i] == "-m") {
   methylated.counts.cols=as.numeric(split(ARG[i+1],","))
 } else if (ARG[i] == "-d") {
       design_file=ARG[i+1]
} else if (ARG[i] == "-g") {
    group_as_control=ARG[i+1]       
} else if (ARG[i] == "-v") {
    variance_estimation=ARG[i+1]       
 } else if (ARG[i] == "-h") {        
    usage("Help:")
 }
}
#   list_files=list.files(path="hg19/methylation", pattern="profile.cg_strand_combined.csv$", full.names=T, include.dirs=T);
#   input_files=paste(list_files,collapse=",")
#   sampleNames=paste(unlist(lapply(basename(list_files), function(x) gsub(".profile.cg_strand_combined.csv","", x, fixed=T))),collapse=",")
#   outFile="DMR"
#   coverage.counts.cols=11
#   methylated.counts.cols=10
#   group_as_control="2"
#   design_file="design.tsv"
#  variance_estimation="same"

library(bsseq) 
library(data.table)

sample_names=unlist(strsplit(sampleNames,","))
input_files=unlist(strsplit(input_files,","))


# Minimun coverage for methylated sites
minCoverage=2

# Default value for variance estimation
if(is.null(variance_estimation)){ variance_estimation='same'}

if(length(input_files)<length(sample_names)) { 
    if(length(input_files)==1){
        warn("ERROR: Sample names and input files are not from the same length, only one input file? "); 
        sampleNames=lapply(sample_names, function(x) sampleNames[[x]]=sample_names)
    }            
}else{
    sampleNames=sample_names
}
# Read design file
# Create BSSeq object, one per sample or one per input file if 1 input file
print (paste ("NOTICE: Creating a BSSeq object from ", input_files, collapse=","))
BSlist <-list()
for ( i in 1:length(input_files)){
         df0=fread(input_files[i], header=T)
         M = matrix(as.numeric(unlist(df0[,methylated.counts.cols,with=FALSE])), ncol=length(methylated.counts.cols))
         Cov = matrix(as.numeric(unlist(df0[,coverage.counts.cols,with=FALSE])), ncol=length(coverage.counts.cols))
         #BSlist[sample_names[i]]=BSseq(chr = unlist(df0[,1,with=FALSE]), pos = as.numeric(unlist(df0[,2,with=FALSE])),M=M, Cov=Cov, sampleNames = sampleNames[i], pData=subset(pData,pData[,1]==sampleNames[i]))
         BSlist[sample_names[i]]=BSseq(chr = unlist(df0[,1,with=FALSE]), pos = as.numeric(unlist(df0[,2,with=FALSE])),M=M, Cov=Cov, sampleNames = sampleNames[i])
}    
# Combine and collapse samples if needed
print(BSlist)

# Just got killed right here... tried in lm nodes
print (paste ("NOTICE: Combining BSseq objects resulting from ", input_files, collapse=","))
BS.all=combineList(BSlist)

print(BS.all)


# Smoothing all
BS.fit=BSmooth(BS.all, parallelBy='chromosome')

# Analyzing WGBS with the bsseq package 
# Computing t-statistics
# Before computing t-statistics, we will remove CpGs with little or no coverage. 
# Add phenotype data at this point, combining lists of bsseq objects is always giving errors related to the pData
pData=read.csv2(design_file, header=T, sep = "\t", na.strings = "0", check.names=F, colClasses =c('character',rep('numeric',unique(count.fields(design_file))-1)))
pData(BS.fit)<-pData
BS.cov<- getCoverage (BS.fit)    
keepLoci.ex<-which(rowSums(BS.cov[,1:ncol(BS.cov)] >=    minCoverage ) >=  minCoverage)
print (paste ("NOTICE: Loci kept ", length(keepLoci.ex), " from ", length(BS.fit)))
BS.fit<-BS.fit[keepLoci.ex,]
# controls<-unlist(subset(pData,Type==group_as_control,1))
# treatment<-unlist(subset(pData,Type != group_as_control,1))
# We are now ready to compute t-statistics, using controls to estimate variance
# But again, the fit object is missing the sampleNames attribute :-( 
if(is.null(sampleNames(BS.fit))) { sampleNames(BS.fit)=sampleNames(BS.all) }

perform_DRM<-function(BS.fit, controls, treatment, variance_estimation, cutoff.quantiles=c(0.025, 0.975), output="BsSeqDRM.tsv", plot.top.regions=50)
{
    # We are now ready to compute t-statistics, using controls to estimate variance
    BS.fit.tstat <- BSmooth.tstat (BS.fit, group1 = treatment, group2=controls, estimate.var = variance_estimation , local.correct = FALSE , verbose = TRUE )
    # Compute manually quantiles as there's an  error if missing values 
    stat="tstat"
    tstats.quantile<-quantile(BS.fit.tstat@stats[, stat], cutoff.quantiles, na.rm=T)
    print(paste("NOTICE: Computing DM regions using tstats cutoff = ", "[",paste(tstats.quantile, collapse=" , "),"]" ))
    dmrs0 <- dmrFinder(BS.fit.tstat , stat=stat , verbose = TRUE, cutoff=tstats.quantile)
    # sort by absolut mean difference    
    dmrs<- dmrs0[order(abs(dmrs0[,"meanDiff"]), decreasing=T),] 
    write.table(dmrs, output ,quote=F,row.names=F,col.names=T,sep="\t")    
    pdf ( file = gsub(".tsv$", ".pdf",output) , width = 10 , height = 8.5 )     
    try(plotManyRegions (BS.fit, dmrs0 [ 1 : max(nrow(dmrs0),plot.top.regions),], extend = 5000 , addRegions = dmrs) )
    dev.off ()
}

# Now ready to find differentially methylated regions in groupd of samples grouped and contrasted in a design file 
name_sample= as.character(as.vector(pData[,1]))
                          
for (i in 2:ncol(pData)) {
    name_folder = paste(outFile, names(pData[i]), sep="/")
     # Create output directory     
    if (!file.exists(name_folder)) { 
        system(paste("mkdir -p",name_folder,sep=" "))
    }
    current_pData=pData[,i]
    subsampleN=name_sample[!(is.na(current_pData))]
    group = as.character(current_pData)[!(is.na(current_pData))]
    controls = unlist(subset(pData, pData[,i] == group_as_control, 1))
    treatment = unlist(subset(pData, pData[,i] != group_as_control, 1))
    cat("NOTICE: Processing contrast \n")
    cat(paste("Name folder: ",name_folder,"\n",sep=""))
    cat(paste("Design : ",paste(subsampleN, group, sep="=",collapse=" ; "),"\n",spe=""))

    # Perform Differentially methylated regions
    perform_DRM(BS.fit, controls, treatment, variance_estimation, output=paste(name_folder, "DRM_stats.tsv", sep="/"))
}
