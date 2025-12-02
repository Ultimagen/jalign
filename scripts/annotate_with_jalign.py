import argparse
import subprocess
import os
from os.path import join as pjoin
import pysam
import tempfile
from joblib import Parallel, delayed


def jalign_cnv_realign(input_cram,range_bed,ref_fasta,output_prefix,mode,min_mismatches):
    cmd = f"python /jalign/dup_cnv_realign.py {input_cram} {range_bed} {ref_fasta} {output_prefix} DUP {min_mismatches}"
    print(cmd)
    subprocess.check_output(cmd, shell=True)


parser = argparse.ArgumentParser(prog="annotate_with_jalign.py", 
                                description="adds jump alignment annotations to the VCF")
parser.add_argument("--input_cram", help="input cram file",required=True,type=str)
parser.add_argument("--ref_fasta", help="reference fasta file", required=False, type=str)
parser.add_argument("--sample_name", help="sample_name", required=True, type=str)
parser.add_argument("--input_vcf", help="input_file", required=True, type=str)
parser.add_argument("--output_vcf", help="output_file", required=True, type=str)

parser.add_argument("--num_jobs", help="number of jobs run in parallel", required=True, type=str)
parser.add_argument("-v","--verbose",type=bool,required=False,default=True,help="""Whether to print debug messages (default: True)""",)

args = parser.parse_args()

if os.path.exists(args.out_folder):
    print(f"output folder exists: {args.out_folder}")
else:
    os.makedirs(args.out_folder)
    print(f"created output folder: {args.out_folder}")

bed_files = []
temp_bed_file = None

# Create a temporary directory that persists
tmpdirname = tempfile.mkdtemp()
print(f"Using temporary directory: {tmpdirname}")
with pysam.VariantFile(args.input_vcf, "r") as vcf_in:
    with pysam.VariantFile(args.output_vcf, "w", header=vcf_in.header) as vcf_out:

        for nrecord,record in enumerate(vcf_in):
            if nrecord % 1000 == 0:
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
Parallel(n_jobs=num_jobs)(delayed(jalign_cnv_realign)(args.input_cram,range_bed,args.ref_fasta,pjoin(tmpdirname,os.path.basename(range_bed)+'.jalign'),args.min_mismatches) for range_bed in bed_files)
with open(pjoin(tmpdirname, "jalign.bed"),'w') as merged_bed_file:
    for bed_file in bed_files:
        with open(bed_file,'r') as bf:
            for line in bf:
                merged_bed_file.write(line)

# out_DEL_jalign_merged_results_folder = pjoin(args.out_folder,'DEL_jalign_merged_results')
# os.makedirs(out_DEL_jalign_merged_results_folder)

# # aggregate bed file
# cmd = f"cat {args.out_folder}/*.jalign.bed > {out_DEL_jalign_merged_results_folder}/{args.sample_name}.DEL.jalign.bed"
# print(cmd)
# subprocess.check_output(cmd, shell=True)

# output = subprocess.check_output(f"ls -1 {args.out_folder}/")
# print(output)
# # merge bam files
# header_bam = subprocess.check_output(f"ls -1 {args.out_folder}/*.bam | head -1",shell=True, text=True).strip()
# print(header_bam)

# cmd = f"samtools view -H {header_bam} > {out_DEL_jalign_merged_results_folder}/{args.sample_name}.DEL.jalign.sam"
# print(cmd)
# subprocess.check_output(cmd, shell=True)

# cmd = f"cat {args.out_folder}/*sam | grep -v -E \"^@\" >> {out_DEL_jalign_merged_results_folder}/{args.sample_name}.DEL.jalign.sam"
# print(cmd)
# subprocess.check_output(cmd, shell=True)

# cmd = f"samtools sort --threads {args.num_jobs} {out_DEL_jalign_merged_results_folder}/{args.sample_name}.DEL.jalign.sam > {out_DEL_jalign_merged_results_folder}/{args.sample_name}.DEL.jalign.bam"
# print(cmd)
# subprocess.check_output(cmd, shell=True)

# cmd = f"samtools index {out_DEL_jalign_merged_results_folder}/{args.sample_name}.DEL.jalign.bam"
# print(cmd)
# subprocess.check_output(cmd, shell=True)

# cmd = f"rm {out_DEL_jalign_merged_results_folder}/{args.sample_name}.DEL.jalign.sam"
# print(cmd)
# subprocess.check_output(cmd, shell=True)

# print("output files are:")
# print(f"{out_DEL_jalign_merged_results_folder}/{args.sample_name}.DEL.jalign.bed")
# print(f"{out_DEL_jalign_merged_results_folder}/{args.sample_name}.DEL.jalign.bam")
# print(f"{out_DEL_jalign_merged_results_folder}/{args.sample_name}.DEL.jalign.bam.bai")
