# Python scripts

Companion python scripts used to process data generated from Next Generation Sequencing pipelines in 
[mugqic_pipelines](https://bitbucket.org/mugqic/mugqic_pipelines)

Require python 2.7.8 and python libs provided by MUGQIC, check python installation scripts

[python.sh](https://bitbucket.org/mugqic/mugqic_pipelines/src/a0c0188ae4a576b03e02fac7cdb53ed90e5f4d42/resources/modules/python.sh?at=master&fileviewer=file-view-default)
[python_libs.sh](https://bitbucket.org/mugqic/mugqic_pipelines/src/a0c0188ae4a576b03e02fac7cdb53ed90e5f4d42/resources/modules/python_lib.sh?at=master&fileviewer=file-view-default)


## Usage

Get help about usage and arguments with 

```
#!bash
NGSTools/scripts/python/<script_name>.py --help
```
         

## fastaAnnotate.py

Annotate a fasta file using tab or comma separated annotation files

```
#!bash
python bin/NGSTools/scripts/python/fastaAnnotate.py -f path_to_transctriptome/Trinity.fa -i path_to_transcriptome/differential_expression/isoforms.TMM.fpkm.matrix.tmp path_to_transcriptome/trinotate/trinotate_annotation_report.tsv.isoforms_blast.tsv -c tId transcript_id -o Trinity_annotated  -x tId transcript_id

# Will create a Trinity_annotated.fa file with fields comming from input files, default separator is \t

```
