import argparse
import subprocess
import os
from os.path import join as pjoin
import sys
from joblib import Parallel, delayed
import pysam

def vcf_to_bed(vcf_record):
    """Convert a VCF record to bed format (chrom, start, end)"""
    chrom = vcf_record.chrom
    start = vcf_record.pos - 1  # VCF is 1-based, BED is 0-based
    # Try to get END from INFO field, otherwise use POS + len(REF) - 1
    end = vcf_record.stop
    return chrom, start, end


def split_vcf_into_bed_chunks(input_vcf, num_chunks, output_dir):
    """Split VCF file into multiple BED chunk files and return list of BED chunk files"""
    vcf_in = pysam.VariantFile(input_vcf, 'r')
    
    # Read all records
    records = list(vcf_in.fetch())
    total_records = len(records)
    
    if total_records == 0:
        print("Warning: No records found in VCF file")
        return []
    
    chunk_size = max(1, (total_records + num_chunks - 1) // num_chunks)
    chunk_bed_files = []
    
    for i in range(num_chunks):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, total_records)
        
        if start_idx >= total_records:
            break
            
        chunk_records = records[start_idx:end_idx]
        
        # Create BED file directly for this chunk
        chunk_bed_file = pjoin(output_dir, f"chunk_{i}.bed")
        chunk_bed_files.append(chunk_bed_file)
        
        with open(chunk_bed_file, 'w') as bed_out:
            for record in chunk_records:
                chrom, start, end = vcf_to_bed(record)
                bed_out.write(f"{chrom}\t{start}\t{end}\n")
    
    vcf_in.close()
    return chunk_bed_files


def jalign_cnv_realign_chunk(input_cram, chunk_bed, ref_fasta, output_prefix, mode, min_mismatches):
    """Run jalign on a BED chunk"""
    cmd = f"python /jalign/cnv_realign.py {input_cram} {chunk_bed} {ref_fasta} {output_prefix} {mode} {min_mismatches}"
    print(cmd)
    subprocess.check_output(cmd, shell=True)


def annotate_vcf_with_jalign_results(input_vcf, jalign_bed_file, output_vcf):
    """Annotate VCF file with jalign results from BED file"""
    vcf_in = pysam.VariantFile(input_vcf, 'r')
    
    # Add new INFO fields to header
    vcf_in.header.info.add("JALIGN_READS", "1", "Integer", "Number of reads supporting CNV from jalign")
    vcf_in.header.info.add("JALIGN_LOWSCORE", "1", "Integer", "Number of low-score reads from jalign")
    vcf_in.header.info.add("JALIGN_LESSREF", "1", "Integer", "Number of reads with score less than reference from jalign")
    vcf_in.header.info.add("JALIGN_STATS", ".", "String", "Additional jalign statistics")
    
    vcf_out = pysam.VariantFile(output_vcf, 'w', header=vcf_in.header)
    
    # Read jalign results into a dictionary keyed by (chrom, start, end)
    jalign_results = {}
    with open(jalign_bed_file, 'r') as bed_in:
        for line in bed_in:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 6:
                continue
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            # Store all additional columns as annotation
            annotation_data = fields[3:]
            jalign_results[(chrom, start, end)] = annotation_data
    
    # Annotate VCF records
    for record in vcf_in.fetch():
        chrom, start, end = vcf_to_bed(record)
        
        if (chrom, start, end) in jalign_results:
            annotation = jalign_results[(chrom, start, end)]
            
            # Parse annotation based on expected format
            # Expected: jump_read_written, jump_read_lowscore, jump_read_lessthanref, [stats...]
            if len(annotation) >= 3:
                record.info['JALIGN_READS'] = int(annotation[0])
                record.info['JALIGN_LOWSCORE'] = int(annotation[1])
                record.info['JALIGN_LESSREF'] = int(annotation[2])
                
            if len(annotation) > 3:
                record.info['JALIGN_STATS'] = ','.join(annotation[3:])
        
        vcf_out.write(record)
    
    vcf_in.close()
    vcf_out.close()


parser = argparse.ArgumentParser(prog="parallel_run_cnv_realign.py", 
                                description="Process VCF file with CNV deletions: splits VCF into chunks, runs jalign on each chunk in parallel, and annotates the VCF with results")
parser.add_argument("--input_vcf", help="Input VCF file with CNV records", required=True, type=str)
parser.add_argument("--input_cram", help="Input CRAM file", required=True, type=str)
parser.add_argument("--ref_fasta", help="Reference FASTA file", required=True, type=str)
parser.add_argument("--out_folder", help="Output folder", required=True, type=str)
parser.add_argument("--sample_name", help="Sample name", required=True, type=str)
parser.add_argument("--min_mismatches", help="The minimal mismatches to make a read be considered", required=False, type=int, default=0)
parser.add_argument("--mode", help="Jalign mode, one of: {DEL,DUP}", required=False, type=str, default="DEL")
parser.add_argument("--num_jobs", help="Number of jobs to run in parallel", required=True, type=int)
parser.add_argument("-v", "--verbose", type=bool, required=False, default=True, help="Whether to print debug messages (default: True)")


args = parser.parse_args()

# Create output folder
if os.path.exists(args.out_folder):
    print(f"Output folder exists: {args.out_folder}")
else:
    os.makedirs(args.out_folder)
    print(f"Created output folder: {args.out_folder}")

# Create temporary folder for chunks
chunks_folder = pjoin(args.out_folder, 'chunks')
os.makedirs(chunks_folder, exist_ok=True)

# Step 1: Split VCF into BED chunks
print(f"Splitting VCF file into {args.num_jobs} BED chunks...")
chunk_bed_files = split_vcf_into_bed_chunks(args.input_vcf, args.num_jobs, chunks_folder)
print(f"Created {len(chunk_bed_files)} BED chunks")

if len(chunk_bed_files) == 0:
    print("No chunks created. Exiting.")
    sys.exit(1)

# Step 2: Run jalign on each chunk in parallel
print("Running jalign on each chunk in parallel...")
Parallel(n_jobs=args.num_jobs)(
    delayed(jalign_cnv_realign_chunk)(
        args.input_cram,
        chunk_bed,
        args.ref_fasta,
        pjoin(chunks_folder, os.path.basename(chunk_bed) + '.jalign'),
        args.mode,
        args.min_mismatches
    ) for chunk_bed in chunk_bed_files
)

# Step 3: Merge results
out_jalign_merged_results_folder = pjoin(args.out_folder, f'{args.mode}_jalign_merged_results')
os.makedirs(out_jalign_merged_results_folder, exist_ok=True)

# Aggregate BED file
print("Merging BED results...")
cmd = f"cat {chunks_folder}/*.jalign.bed > {out_jalign_merged_results_folder}/{args.sample_name}.{args.mode}.jalign.bed"
print(cmd)
subprocess.check_output(cmd, shell=True)

# Merge BAM files
print("Merging BAM files...")
header_bam = subprocess.check_output(f"ls -1 {chunks_folder}/*.bam | head -1", shell=True, text=True).strip()
print(f"Using header from: {header_bam}")

cmd = f"samtools view -H {header_bam} > {out_jalign_merged_results_folder}/{args.sample_name}.{args.mode}.jalign.sam"
print(cmd)
subprocess.check_output(cmd, shell=True)

cmd = f"cat {chunks_folder}/*.sam | grep -v -E \"^@\" >> {out_jalign_merged_results_folder}/{args.sample_name}.{args.mode}.jalign.sam"
print(cmd)
subprocess.check_output(cmd, shell=True)

cmd = f"samtools sort --threads {args.num_jobs} {out_jalign_merged_results_folder}/{args.sample_name}.{args.mode}.jalign.sam > {out_jalign_merged_results_folder}/{args.sample_name}.{args.mode}.jalign.bam"
print(cmd)
subprocess.check_output(cmd, shell=True)

cmd = f"samtools index {out_jalign_merged_results_folder}/{args.sample_name}.{args.mode}.jalign.bam"
print(cmd)
subprocess.check_output(cmd, shell=True)

cmd = f"rm {out_jalign_merged_results_folder}/{args.sample_name}.{args.mode}.jalign.sam"
print(cmd)
subprocess.check_output(cmd, shell=True)

# Step 4: Annotate VCF with jalign results
print("Annotating VCF with jalign results...")
output_vcf = pjoin(out_jalign_merged_results_folder, f"{args.sample_name}.{args.mode}.jalign.vcf")
merged_bed = pjoin(out_jalign_merged_results_folder, f"{args.sample_name}.{args.mode}.jalign.bed")
annotate_vcf_with_jalign_results(args.input_vcf, merged_bed, output_vcf)

print("\n=== Output files ===")
print(f"Annotated VCF: {output_vcf}")
print(f"BED file: {merged_bed}")
print(f"BAM file: {out_jalign_merged_results_folder}/{args.sample_name}.{args.mode}.jalign.bam")
print(f"BAM index: {out_jalign_merged_results_folder}/{args.sample_name}.{args.mode}.jalign.bam.bai")

