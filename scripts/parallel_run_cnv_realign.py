import argparse
import subprocess
import os
import glob
from os.path import join as pjoin
import sys
from joblib import Parallel, delayed


def jalign_cnv_realign(input_cram,range_bed,ref_fasta,output_prefix):
    cmd = f"python /jalign/cnv_realign.py {input_cram} {range_bed} {ref_fasta} {output_prefix}"
    print(cmd)
    subprocess.check_output(cmd, shell=True)


parser = argparse.ArgumentParser(prog="parallel_run_cnv_align.py", 
                                description="splits input bed file into smaller bed files and runs jalign on each of them in parallel. also merge results into a single output file per filte type")
parser.add_argument("--folder_with_cnv_del_bed_files", help="folder with cnv bed files. this script will run jalign for each *.bed file in the folder", required=True, type=str)
parser.add_argument("--input_cram", help="input cram file",required=True,type=str)
parser.add_argument("--ref_fasta", help="reference fasta file", required=False, type=str)
parser.add_argument("--out_folder", help="output folder", required=True, type=str)
parser.add_argument("--sample_name", help="sample_name", required=True, type=str)
parser.add_argument("--num_jobs", help="number of jobs run in parallel", required=True, type=str)
parser.add_argument("-v","--verbose",type=bool,required=False,default=True,help="""Whether to print debug messages (default: True)""",)

args = parser.parse_args()

if os.path.exists(args.out_folder):
    print(f"output folder exists: {args.out_folder}")
else:
    os.makedirs(args.out_folder)
    print(f"created output folder: {args.out_folder}")

range_bed_files = glob.glob(args.folder_with_cnv_del_bed_files+'/*.bed')
print(f"found {len(range_bed_files)} bed files in {args.folder_with_cnv_del_bed_files}")   

# run jalign
num_jobs = args.num_jobs
Parallel(n_jobs=num_jobs)(delayed(jalign_cnv_realign)(args.input_cram,range_bed,args.ref_fasta,pjoin(args.out_folder,os.path.basename(range_bed)+'.jalign')) for range_bed in range_bed_files)

out_DEL_jalign_merged_results_folder = pjoin(args.out_folder,'DEL_jalign_merged_results')
os.makedirs(out_DEL_jalign_merged_results_folder)

# aggregate bed file
cmd = f"cat {args.out_folder}/*.jalign.bed > {out_DEL_jalign_merged_results_folder}/{args.sample_name}.DEL.jalign.bed"
print(cmd)
subprocess.check_output(cmd, shell=True)

# merge bam files
header_bam = subprocess.check_output(f"ls -1 {args.out_folder}/*.bam | head -1",shell=True, text=True).strip()
print(header_bam)

cmd = f"samtools view -H {header_bam} > {out_DEL_jalign_merged_results_folder}/{args.sample_name}.DEL.jalign.sam"
print(cmd)
subprocess.check_output(cmd, shell=True)

cmd = f"cat {args.out_folder}/*sam | grep -v -E \"^@\" >> {out_DEL_jalign_merged_results_folder}/{args.sample_name}.DEL.jalign.sam"
print(cmd)
subprocess.check_output(cmd, shell=True)

cmd = f"samtools sort --threads {args.num_jobs} {out_DEL_jalign_merged_results_folder}/{args.sample_name}.DEL.jalign.sam > {out_DEL_jalign_merged_results_folder}/{args.sample_name}.DEL.jalign.bam"
print(cmd)
subprocess.check_output(cmd, shell=True)

cmd = f"samtools index {out_DEL_jalign_merged_results_folder}/{args.sample_name}.DEL.jalign.bam"
print(cmd)
subprocess.check_output(cmd, shell=True)

cmd = f"rm {out_DEL_jalign_merged_results_folder}/{args.sample_name}.DEL.jalign.sam"
print(cmd)
subprocess.check_output(cmd, shell=True)

print(f"output files are:")
print(f"{out_DEL_jalign_merged_results_folder}/{args.sample_name}.DEL.jalign.bed")
print(f"{out_DEL_jalign_merged_results_folder}/{args.sample_name}.DEL.jalign.bam")
print(f"{out_DEL_jalign_merged_results_folder}/{args.sample_name}.DEL.jalign.bam.bai")
