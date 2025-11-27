#!/usr/bin/env python3

# dup_cnv_realign.py - shell for jump alignment 

import sys
import os
import pysam
import subprocess
import re
import random

# predictive results
random.seed(0)

# this code only supports DUP
MODE = "DUP"

# aligner weights
MATCH_SCORE = 2
MISMATCH_SCORE = -8
OPEN_SCORE = -18
EXTEND_SCORE = -1
JUMP_SCORE = 0

# parametrization
MIN_MISMATCHES = 30
SOFTCLIP_THRESHOLD = 30
FETCH_READ_PADDING = 500
FETCH_REF_PADDING = 0
MIN_SEQ_LEN_JUMP_ALIGN_COMPONENT = 30
MIN_GAP_LEN = 30
MAX_READS_PER_CNV = 4000

# tmp file (prefix) for communicating w/ C aligner
# longer tmp names contain the cnv region itself in the filename
tmp = "/tmp/jump_align_input." + str(os.getpid())
LONG_TMP_NAME = False

# alignment tool (last one overrides)
TOOL = "para_jalign"
TOOL = "jump_align"

# template for C aligner command line
JUMP_ALIGN_CMD = [TOOL, str(MATCH_SCORE), str(MISMATCH_SCORE), str(OPEN_SCORE), str(EXTEND_SCORE), "1", str(JUMP_SCORE), tmp]

# some adjustments for running/debugging locally (outside of a docker)
LOCAL_JUMP_ALIGN = "jump_align/" + TOOL
if os.path.exists(LOCAL_JUMP_ALIGN):
    JUMP_ALIGN_CMD[0] = LOCAL_JUMP_ALIGN
    LONG_TMP_NAME = True

# parse commandline
if len(sys.argv) < 5:
    print("usage: " + sys.argv[0] + " <input-cram> <range-bed> <ref-fasta> <output-prefix> [<mode>] [<min-mismatches] [<fetch-read-padding>]\n")
    sys.exit(-1)

# commandline invocation
IN_CRAM = sys.argv[1]
CNV_BED = sys.argv[2]
REF_FASTA = sys.argv[3]
OUT_SAM = sys.argv[4] + ".sam"
if len(sys.argv) >= 6:
    MODE = sys.argv[5]
if len(sys.argv) >= 7:
    MIN_MISMATCHES = int(sys.argv[6])
if len(sys.argv) >= 8:
    FETCH_READ_PADDING = int(sys.argv[7])

assert MODE == "DUP"

# set up filename for output bed file (main result)
# sam result is supplemental (debugging aid)
OUT_BED = OUT_SAM.replace(".sam", ".bed")
OUT_LOG = OUT_SAM.replace(".sam", ".log")

# open input files
reads_file = pysam.AlignmentFile(IN_CRAM, "rb", reference_filename=REF_FASTA)
fasta_file = pysam.FastaFile(REF_FASTA)
chrom_sizes = dict(zip(fasta_file.references, fasta_file.lengths))

# Runs an external command, captures its output and exit code, and 
# returns the command’s stdout along with its return status.
def run_process(command, flog):
  print(command)
  try:
    process = subprocess.Popen(command, 
                              stdin=subprocess.PIPE, 
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE, 
                              text=True) 

    stdout, stderr = process.communicate()
    returncode = process.returncode
    #print(stderr)

  except subprocess.CalledProcessError as e:
    print(f"Error executing command: {e}")
    return None, e.returncode

  return stdout, returncode


def count_md_mismatches(read):
    """
    Count the number of mismatches in a pysam.AlignedSegment read using the MD tag.
    Indels are not counted.
    """
    try:
        md_tag = read.get_tag("MD")
    except KeyError:
        return None

    # Find all letters in the MD string, which represent mismatches
    mismatches = re.findall(r"[A-Z]", md_tag)
    return len(mismatches)

def count_nm_mismatches(read):
    """
    Count the number of mismatches in a pysam.AlignedSegment read using the NM tag.
    Indels are counted.
    """
    try:
        nm_tag = read.get_tag("NM")
    except KeyError:
        return None

    return int(nm_tag)

def count_softclip_mismatches(read, reference):
    """
    Count mismatches in soft-clipped regions (both left and right) of a read.
    `read` is a pysam.AlignedSegment
    `reference` is a pysam.FastaFile
    """
    if read.is_unmapped:
        return 0

    seq = read.query_sequence
    mismatches = 0
    ref_name = read.reference_name
    start = read.reference_start
    end = read.reference_end

    cigartuples = read.cigartuples
    # CIGAR operations
    SOFT_CLIP = 4

    # Left soft clip
    if cigartuples[0][0] == SOFT_CLIP:
        clip_len = cigartuples[0][1]
        clipped_bases = seq[:clip_len]
        ref_start = max(0, start - clip_len)
        ref_bases = reference.fetch(ref_name, ref_start, start)
        for rb, qb in zip(ref_bases, clipped_bases):
            if rb.upper() != qb.upper():
                mismatches += 1

    # Right soft clip
    if cigartuples[-1][0] == SOFT_CLIP:
        clip_len = cigartuples[-1][1]
        clipped_bases = seq[-clip_len:]
        ref_bases = reference.fetch(ref_name, end, end + clip_len)
        for rb, qb in zip(ref_bases, clipped_bases):
            if rb.upper() != qb.upper():
                mismatches += 1

    return mismatches

# is read soft clipped
def is_softclipped(read):
    return read.cigartuples[0][0] == pysam.CSOFT_CLIP or read.cigartuples[-1][0] == pysam.CSOFT_CLIP

def is_substential_softclipped(read):
    return (read.cigartuples[0][0] == pysam.CSOFT_CLIP and read.cigartuples[0][1]) >= SOFTCLIP_THRESHOLD \
                or (read.cigartuples[-1][0] == pysam.CSOFT_CLIP and read.cigartuples[-1][1] >= SOFTCLIP_THRESHOLD)

# a filtering function to see if a read is worth "accepting" into the set to reads to consider
def accept_read(read):
    if MIN_MISMATCHES <= 0:
        return True
    sc = count_softclip_mismatches(read, fasta_file);
    nm = count_nm_mismatches(read)
    return (sc + nm) >= MIN_MISMATCHES

# process a single cnv region (essentially a line from an input bed file)
def process_cnv(chrom, start, end, flog):

    if flog:
        flog.write(">>> %s:%d-%d\n" % (chrom, start, end))

    # get all reads that cross the two cnv edges
    reads = dict()
    reads_in_ref = [set(), set()]
    refs = []
    refs_extents = []
    ref_id = 0
    for loc in [start, end]:
        rmin = max(0, loc - FETCH_READ_PADDING)
        rmax = loc + FETCH_READ_PADDING
        # print("fetching ... ", max(0, loc - FETCH_READ_PADDING), loc + FETCH_READ_PADDING)
        for read in reads_file.fetch(chrom, max(0, loc - FETCH_READ_PADDING), loc + FETCH_READ_PADDING):
            if not is_substential_softclipped(read):
                continue
            if not accept_read(read):
                continue
            #if mode == "DUP" and read.is_supplementary:
            #    continue
            reads[read.qname] = read
            rmin = min(rmin, read.reference_start)
            rmax = max(rmax, read.reference_end)
            reads_in_ref[ref_id].add(read.qname)
        refs_extents.append([rmin, rmax])
        ref_id += 1
        
    # extend references before and after
    refs_extents[0][0] = max(0, refs_extents[0][0] - FETCH_REF_PADDING)
    refs_extents[1][1] = refs_extents[1][1] + FETCH_REF_PADDING

    # get references
    for extents in refs_extents:
        rmin, rmax = extents
        ref = fasta_file.fetch(chrom, rmin, rmax)
        refs.append([rmin, ref])

    # create input file for jump aligner
    ref_emitted = False
    reads_in_order = []
    subsample_ratio = 1.0
    if len(reads) > MAX_READS_PER_CNV:
        subsample_ratio = MAX_READS_PER_CNV / len(reads)
        print("subsample_ratio", subsample_ratio)
    nlines = 0
    jalign_input = tmp
    if LONG_TMP_NAME:
        jalign_input += "_" + chrom + ":" + str(start) + "-" + str(end)
    with open(jalign_input, 'w') as f:
        for read in reads.values():
            if subsample_ratio < 1.0:
                if random.random() > subsample_ratio:
                    continue
            if not accept_read(read):
                continue
            reads_in_order.append(read)
            if not ref_emitted:
                line = read.qname + "\t" + read.seq + "\t" + refs[1][1] + "\t" + refs[0][1] + "\n"
                ref_emitted = True
            else:
                line = read.qname + "\t" + read.seq + "\t=\n"
            f.write(line)
            nlines += 1
            if flog:
                flog.write(line)

    # run jump_align
    JUMP_ALIGN_CMD[-1] = jalign_input
    if flog:
        flog.write("<<< %s\n" % (JUMP_ALIGN_CMD))
    alignments = run_process(JUMP_ALIGN_CMD, flog)
    header_seen = False
    realignments = []
    rheader = []
    for alignment, read in zip(alignments[0].split("\n"), [None, *reads_in_order]):
        if flog:
            flog.write(alignment)
        if not header_seen:
            rheader = alignment.split("\t")
            header_seen = True;
        else:
            a = alignment.split("\t")
            in1 = read.qname in reads_in_ref[0]
            in2 = read.qname in reads_in_ref[1]
            realignments.append([read, refs[0][0], refs[1][0], a, in1, in2])
    return (rheader, realignments, nlines)    

# main loop starts here
# read each line from the bed file, call process_cnv, write output bed file
with open(OUT_BED, "w") as out_bed, open(OUT_LOG, "w") as flog:
    with open(CNV_BED) as f:
        for line in f:

            # parse line
            if line.startswith("#"):
                continue
            bed_line = line.strip().split()
            bed_chrom, bed_start, bed_end = bed_line[:3]
            bed_start = int(bed_start)
            bed_end = int(bed_end)

            # check for valid end locus 
            if bed_end + FETCH_READ_PADDING > chrom_sizes[bed_chrom]:
                continue

            # process this line
            rheader, realignments, nlines = process_cnv(bed_chrom, bed_start, bed_end, flog)

            # collect better (DUP) jump alignment and (DEL) jump alignments
            # variables with natural names relate to DUP
            # variables with a d prefix relate to DEL
            jump_better = 0
            djump_better = 0

            for realignment in realignments:
                in_ref = [False, False]
                read, ref1_start, ref2_start, ainfo, in_ref[0], in_ref[1] = realignment
    
                # decode alignment info
                if TOOL == "jump_align":
                    readName,\
                        score,jumpInsertSize,jumpRange,jbegin1,japath1,jreadlen1,jreflen1,jbegin2,japath2,jreadlen2,jreflen2,\
                        dscore,djumpInsertSize,djumpRange,djbegin1,djapath1,djreadlen1,djreflen1,djbegin2,djapath2,djreadlen2,djreflen2,\
                        score1,begin1,apath1,readlen1,reflen1,\
                        score2,begin2,apath2,readlen2,reflen2 = ainfo
                    score,jumpInsertSize,jumpRange,jbegin1,jreadlen1,jreflen1,jbegin2,jreadlen2,jreflen2,\
                        dscore,djumpInsertSize,djumpRange,djbegin1,djreadlen1,djreflen1,djbegin2,djreadlen2,djreflen2,\
                        score1,begin1,readlen1,reflen1,\
                        score2,begin2,readlen2,reflen2 \
                        = map(int, (score,jumpInsertSize,jumpRange,jbegin1,jreadlen1,jreflen1,jbegin2,jreadlen2,jreflen2,dscore,djumpInsertSize,djumpRange,djbegin1,djreadlen1,djreflen1,djbegin2,djreadlen2,djreflen2,score1,begin1,readlen1,reflen1,score2,begin2,readlen2,reflen2))

                    # jump score better?
                    if score > 0 and score > max(score1, score2) + MIN_SEQ_LEN_JUMP_ALIGN_COMPONENT:
                        if min(jreadlen1, jreadlen1) >= MIN_SEQ_LEN_JUMP_ALIGN_COMPONENT:
                            jump_better += 1

                    if dscore > 0 and dscore > max(score1, score2) + MIN_SEQ_LEN_JUMP_ALIGN_COMPONENT:
                        if min(djreadlen1, djreadlen1) >= MIN_SEQ_LEN_JUMP_ALIGN_COMPONENT:
                            djump_better += 1

                if TOOL == "para_align":
                    qname1, better, score, score1, score2, jgain, size1, size2, \
                    dscore, dscore1, dscore2, djgain, dsize1, dsize2 = ainfo              
                    score, score1, score2, size1, size2 = map(int, (score, score1, score2, size1, size2))
                    dscore, dscore1, dscore2, dsize1, dsize2 = map(int, (dscore, dscore1, dscore2, dsize1, dsize2))

                    # jump score better?
                    if score > 0 and score > max(score1, score2) + MIN_SEQ_LEN_JUMP_ALIGN_COMPONENT:
                        if min(size1, size2) >= MIN_SEQ_LEN_JUMP_ALIGN_COMPONENT:
                            jump_better += 1

                    if dscore > 0 and dscore > max(dscore1, dscore2) + MIN_SEQ_LEN_JUMP_ALIGN_COMPONENT:
                        if min(dsize1, dsize2) >= MIN_SEQ_LEN_JUMP_ALIGN_COMPONENT:
                            djump_better += 1

            
            outline = line[:-1] + ("\t%d\t%d\n" % (jump_better, djump_better))
            out_bed.write(outline)
            print(outline.strip())
            del outline
