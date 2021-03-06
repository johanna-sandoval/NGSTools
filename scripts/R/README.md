# Running R scripts

## DMR using BiSeq

Code to process CpG methylation counts generated by the EMC-bwa pipeline

```
# Run Biseq using a capture bed file (SeqCap_EPI_4M_CPGiant.targets.bed)
infiles=`find hg19/methylation | grep profile.cg_strand_combined.csv | tr '\n' ',' | sed 's/,$//'`
sampleNames=`find hg19/methylation | grep profile.cg_strand_combined.csv  | awk -F"/" '{split($3,a,"."); print a[1]}'| tr '\n' ',' | sed 's/,$//'`
module load mugqic_dev/R_Bioconductor && mkdir -p DMR/BiSeqDMR && Rscript BiSeqDMR.R -i $infiles -s $sampleNames -d design.tsv -o DMR/BiSeqDMR -c 11 -m 10 -p 8 -t SeqCap_EPI_4M_CPGiant.targets.bed

# Plot BiSeq best DMR 
module load mugqic_dev/R_Bioconductor && R --no-save < BiSeqPlotHeatmap.R 

# Add annotations and redo reports
module load mugqic_dev/R_Bioconductor/3.2.2_3.1 mugqic_dev/mugqic_R_packages/master && R --no-save < BiSeq_annotate.R

# Run exploratory analysis
module load mugqic_dev/R_Bioconductor/3.2.2_3.1 mugqic_dev/mugqic_R_packages/master && R --no-save < bin/BiSeqExploratory.R

```

## barcode/haipins counts using EDGER
Code to process hairpin reads from Illumina sequencer. Here we assume basic fixed structure of read:
# Barcode + Common sequence + Hairpin sequence
fastq is a path to a fastq file , must be located in the working directory
The input barcode file and hairpin/sgRNA files are tab-separated text files with at least two columns
(named ’ID’ and ’Sequences’) containing the sample or hairpin/sgRNA ids and a second column indicating 
the sample index or hairpin/sgRNA sequences to be matched.
These files, along with the fastq file/(s) are assumed to be in the current working directory. 
Any additional column in the hairpin (fi) file will be added to the counts output

```
# Run edgeR function (take ~7 hours for 80K hairpins 1 barcode :-\ )
Rscript edgerProcessAmplicons.R -r ${fastq} -b ${barcodes} -g ${fi} -s 1 -e 8 -t 32
   
# Collect results
outdirs=`ls -d *out/*out | tr '\n', ',' | sed 's/.$//'`
module load mugqic_dev/R_Bioconductor && Rscript bin/collectResultsProcessAmplicons.R -o $outdirs 

# Report using markdown when 

rm edgerProcessAmpliconsReport.Rmd && module load mugqic/R_Bioconductor/3.2.3_3.2 mugqic/pandoc/1.15.2 && echo "library(rmarkdown); report_dir='report'; source_dir='Genome-wide_mouse_lentiviral_CRISPR_gRNA_dilution.MPS12300165-C03.M03555_0138.1.single_hairpins.tsv_out/hairpins.tsv_out/,Genome-wide_mouse_lentiviral_CRISPR_gRNA.MPS12300153-A01.M03992_0084.1.single_hairpins.tsv_out/hairpins.tsv_out/'; log_files='job_output/edgerProcessAmplicons/edgerProcessAmplicons.Genome-wide_mouse_lentiviral_CRISPR_gRNA_dilution.MPS12300165-C03.M03555_0138.1.single_hairpins_2016-08-23T11.23.25.o,job_output/edgerProcessAmplicons/edgerProcessAmplicons.Genome-wide_mouse_lentiviral_CRISPR_gRNA.MPS12300153-A01.M03992_0084.1.single_hairpins_2016-08-23T11.23.25.o'; file.copy('bin/edgerProcessAmpliconsReport.Rmd', './'); rmarkdown::render(input='edgerProcessAmpliconsReport.Rmd', output_format = c('html_document'), output_dir = 'report');" | R --vanilla

```


### Detecting differential binding in ChIP-seq data with csaw
A fairly modular process, as explained [here](https://www.bioconductor.org/help/course-materials/2015/BioC2015/csaw_lab.html)

```
module load mugqic_dev/R_Bioconductor && Rscript runCsaw.R -d design.dbr.tsv
```
