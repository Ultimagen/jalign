#!/usr/bin/env python3
 
import argparse
import subprocess
import os
from os.path import join as pjoin
from os.path import dirname, abspath
import pysam
import tempfile
from joblib import Parallel, delayed

def jalign_cnv_realign(input_cram,range_bed,ref_fasta,output_prefix,min_mismatches):
    cmd = f"python /jalign/dup_cnv_realign.py {input_cram} {range_bed} {ref_fasta} {output_prefix} DUP {min_mismatches}"
    print(cmd)
    subprocess.check_output(cmd, shell=True)


parser = argparse.ArgumentParser(prog="annotate_with_jalign.py", 
                                description="adds jump alignment annotations to the VCF")
parser.add_argument("--input_cram", help="input cram file",required=True,type=str)
parser.add_argument("--ref_fasta", help="reference fasta file", required=False, type=str)
parser.add_argument("--input_vcf", help="input_file", required=True, type=str)
parser.add_argument("--output_vcf", help="output_file", required=True, type=str)
parser.add_argument("--min_mismatches", help="minimum mismatches for jalign", required=False, type=int, default=5)
parser.add_argument("--num_jobs", help="number of jobs run in parallel", required=True, type=str)
parser.add_argument("-v","--verbose",type=bool,required=False,default=True,help="""Whether to print debug messages (default: True)""",)

args = parser.parse_args()

out_folder = dirname(abspath(args.output_vcf))

if os.path.exists(out_folder):
    print(f"output folder exists: {out_folder}")
else:
    os.makedirs(out_folder)
    print(f"created output folder: {out_folder}")

bed_files = []
temp_bed_file = None

# Create a temporary directory that persists
tmpdirname = tempfile.mkdtemp(dir=out_folder)
print(f"Using temporary directory: {tmpdirname}")
with pysam.VariantFile(args.input_vcf, "r") as vcf_in:
    with pysam.VariantFile(args.output_vcf, "w", header=vcf_in.header) as vcf_out:

        for nrecord,record in enumerate(vcf_in):
            if nrecord % 200 == 0:
                if temp_bed_file is not None:
                    temp_bed_file.close()
                    bed_files.append(temp_bed_file.name)
                
                temp_bed_file = tempfile.NamedTemporaryFile(dir=tmpdirname, mode='w+', suffix=".bed", delete=False)
            if temp_bed_file is None:
                raise RuntimeError("temp_bed_file is None - should not happen")
            temp_bed_file.write(f"{record.chrom}\t{record.pos-1}\t{record.stop}\n")
        
        # Close the last bed file
        if temp_bed_file is not None:
            temp_bed_file.close()
            bed_files.append(temp_bed_file.name)

# run jalign
num_jobs = args.num_jobs

Parallel(n_jobs=num_jobs, verbose=10)(delayed(jalign_cnv_realign)(args.input_cram,range_bed,args.ref_fasta,pjoin(tmpdirname,os.path.basename(range_bed)+'.jalign'),args.min_mismatches) 
                          for range_bed in bed_files)
with open(pjoin(tmpdirname, "jalign.bed"),'w') as merged_bed_file:
    for bed_file in bed_files:
        with open(f"{bed_file}.jalign.bed",'r') as bf:
            for line in bf:
                merged_bed_file.write(line)

annotated_df = open(merged_bed_file.name).readlines()
annotated_df_columns = dict(zip(['chrom','start','end','JALIGN_DUP_SUPPORT','JALIGN_DEL_SUPPORT','JALIGN_DUP_SUPPORT_STRONG','JALIGN_DEL_SUPPORT_STRONG'], range(7)))

with pysam.VariantFile(args.input_vcf, "r") as vcf_in:
    header = vcf_in.header
    header.info.add("JALIGN_DUP_SUPPORT", 1, 'Integer', 'Number of reads supporting the duplication via jump alignment')
    header.info.add("JALIGN_DEL_SUPPORT", 1, 'Integer', 'Number of reads supporting the deletion via jump alignment')
    header.info.add("JALIGN_DUP_SUPPORT_STRONG", 1, 'Integer', 'Number of reads strongly supporting the duplication via jump alignment')
    header.info.add("JALIGN_DEL_SUPPORT_STRONG", 1, 'Integer', 'Number of reads strongly supporting the deletion via jump alignment')    
    with pysam.VariantFile(args.output_vcf, "w", header=vcf_in.header) as vcf_out:
        for nrecord,record in enumerate(vcf_in):
            # get corresponding annotations
            ann_row = annotated_df[nrecord].strip().split("\t")
            # add annotations to record
            record.info['JALIGN_DUP_SUPPORT'] = int(ann_row[annotated_df_columns['JALIGN_DUP_SUPPORT']])
            record.info['JALIGN_DEL_SUPPORT'] = int(ann_row[annotated_df_columns['JALIGN_DEL_SUPPORT']])
            record.info['JALIGN_DUP_SUPPORT_STRONG'] = int(ann_row[annotated_df_columns['JALIGN_DUP_SUPPORT_STRONG']])
            record.info['JALIGN_DEL_SUPPORT_STRONG'] = int(ann_row[annotated_df_columns['JALIGN_DEL_SUPPORT_STRONG']])
            vcf_out.write(record)

