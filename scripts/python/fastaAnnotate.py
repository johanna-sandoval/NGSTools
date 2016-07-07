#!/usr/bin/env python

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

import csv
import argparse
import string
import sys
import os
import collections
import warnings
import re
from Bio import SeqIO


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", "--fasta_file", help="Fasta file to annotate" , type=str, required=True)
    parser.add_argument("-i", "--input_files", help="List of annotation files that will be merged to the fasta file" , type=str, nargs="+", required=True)
    parser.add_argument("-d", "--delimiter", help="Delimiter for input annotation files" , type=str, nargs="+", required=False, default='\t')    
    parser.add_argument("-o", "--output", help="Output File prefix", type=str , required=True)
    parser.add_argument("-c", "--common", help="Columns in input files containing the fasta IDs. If more than one, use comma separated column names, i.e. gene,transcript ", nargs="+",  type=str , required=False)
    parser.add_argument("-s", "--subset", help="A subset of column names to add to the fasta annotation", nargs="+",  type=str , required=False, default=None)
    parser.add_argument("-x", "--exclude", help="A subset of column names to exclude from  the fasta annotation", nargs="+",  type=str , required=False, default=None)
    parser.add_argument("-l", "--left", help="If selected, left outer join of input annotation files", action="store_true")
    parser.add_argument("-n", "--make_names", help="If any value, then the names of the columns in the input files that are duplicated are adjusted so that they are all preserved. If a list of names is provided, this list will be used in the final file", nargs="+",  type=str , required=False, default=None )
    
    
    
    args = parser.parse_args()
    
    #error found when testing, don't ask me if system size is exceeded
    csv.field_size_limit(sys.maxsize)
    
    # Parameters
    infiles=args.input_files
    outfile=args.output
    subset=args.subset
    exclude=args.exclude
    
    # Separator for output file is \t
    out_sep='\t'
    in_sep=[]
    
    # Separator , common must have the same length as infiles
    if len(args.delimiter) == len(infiles) and len(args.delimiter) > 1:
        in_sep=args.delimiter
    elif len(args.delimiter) == 1 :
        in_sep=[args.delimiter[0]] * len(infiles)
    else:
        raise Exception("Error: Two or more input files are required " + args.delimiter )
 
    if not args.common is None and len(args.common) == len(infiles) and len(args.common) > 1:
        key=args.common
    elif not args.common is None  and len(args.common) == 1:
        key=[args.common[0]] * len(infiles)
    else:
        key=[None] * len(infiles)
    if args.make_names:
        if len(args.make_names) == len(infiles) :
            suffix_names=args.make_names
        else:
            suffix_names=range(len(infiles))
            
    print "NOTICE: merging " + " ".join(infiles) + ", delimited by: " + str(in_sep) + ", with keys: " +  str(key) 


    # Open input files, parse and merge using the specified column
    data = collections.OrderedDict()
    #data = {}
    fieldnames = []
    replacement_names={}   
    i=0
    for filename in infiles:
        with open(filename, "rb") as fp: # python 2
            if in_sep[i] in ['\t',"\\t"]:                
                reader = csv.DictReader(fp, quoting=csv.QUOTE_NONE, dialect=csv.excel_tab)
            else:
                dialect = csv.Sniffer().sniff(fp.read(1024))
                fp.seek(0)
                reader = csv.DictReader(fp, delimiter=in_sep[i], dialect=dialect)
            if args.make_names:
                new_names = dict((rn,rn + "_" + suffix_names[i]) if ( rn in fieldnames and rn not in key[i].split(",") ) else (rn,rn) for rn in reader.fieldnames)
                #print str([new_names.values()])
                #print str([new_names.keys()])
                if i==0 :
                    replacement_names = { v:k for k, v in new_names.items() if k in fieldnames }
                else:
                    replacement_names.update( { v:k for k, v in new_names.items() if k in fieldnames } )
                reader.fieldnames = [new_names[fn] for fn in reader.fieldnames]
            fieldnames.extend(reader.fieldnames)
            for row in reader:
                key_field = out_sep.join(row[k] for k in key[i].split(",") ) if key[i] else reader.fieldnames[0]
                # Default behavior is cross join (common elements of all tables are added)
                # if left outer join, preserve only keys of the first file
                if (args.left and i == 0 ) or not args.left:
                    data.setdefault(key_field , {}).update(row)
                elif data.has_key(key_field) :
                    data[key_field].update(row)
            del reader
        i+=1
    
    # Convert to OrderedDict
    # data=collections.OrderedDict(data)
    # Subset columns if required
    if not subset is None:
        diff = [x for x in subset if x not in fieldnames]
        if len(diff) > 0 :
            error_subset = " ".join(diff)        
            warnings.warn("WARNING: Columns to include are not in input files: all columns will be included " + error_subset )
        else:    
            fieldnames = subset
            if args.make_names:
                fieldnames.extend([k for k,v in replacement_names.items() if v in fieldnames])
    else:
      fieldnames = list(data.fromkeys(fieldnames))
    if exclude is None:
        diff = []        
    else:
        diff = [x for x in exclude if x not in fieldnames]
    
    # Exclude columns if required
    if exclude and len(diff) > 0 :
        error_exclude = " ".join(diff)        
        warnings.warn("WARNING: Columns to exclude are not in input files: " + error_exclude )
    elif exclude is None:
        pass
    else:        
        if args.make_names:
            exclude.extend([k for k,v in replacement_names.items() if v in exclude])
        fieldnames = [x for x in fieldnames if x not in exclude ]
    print "NOTICE: will preserve the following headers: " + str(fieldnames)
    
    # Select unique IDs
    if args.make_names:
        unique_ids=list(collections.OrderedDict.fromkeys(fieldnames))
    else:
        unique_ids=fieldnames
        
    # Read Input fasta 
    output_file=open(args.output + ".fa",'w')
    seq_records_FF=[]
    for seq_record in SeqIO.parse(args.fasta_file, "fasta"):
        if data.has_key(seq_record.id) :
            seq_record.description = seq_record.description  + out_sep + "len=" + str(len(seq_record.seq)) + out_sep + out_sep.join(map(lambda (x,y): str(x) + "=" + str(y), zip(unique_ids, [data[seq_record.id].get(field, '') for field in unique_ids])))
        # Write Fasta and tabular output
        SeqIO.write(seq_record, output_file, 'fasta')

    output_file.close()
