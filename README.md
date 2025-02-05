# jalign
Jump Alignment tool assisting in SV calls

## cnv_realign

This tool attempts to suggest alternative alignments to reads by using two reference areas and allowing the read's alignment to span both. Please see data folder for examples.

The core of the tool is a jump aligner borrowed from the manta project. It is wrapped and called by a python notebook.

### building

<pre>
$ docker build . --tag jalign
</pre>

### running

<pre>
$ docker run -it --rm \
    -v .:/here -v /tmp:/out -v ~/tmp/ref:/ref \
    jalign \
    \
    ./cnv_realign.py \
    /here/data/chr22_1M.cram \
    /here/data/chr22.bed \
    /ref/Homo_sapiens_assembly38.fasta \
    /out/cnv_realign_out

$ ls -1 /tmp/cnv_realign_out.*
/tmp/cnv_realign_out.bam
/tmp/cnv_realign_out.bam.bai
/tmp/cnv_realign_out.bed
/tmp/cnv_realign_out.sam

$ head -3 /tmp/cnv_realign_out.bed
chr22	22762500	22763500	CN1	TP	0	4	19 ...
chr22	22765000	22768000	CN1	TP	0	0	0 ...
chr22	22771500	22773000	CN1	TP	0	0	2 ...
</pre>

### output

Two main files are generated: .bed and .bam

bed file is a copy of the input bed file, with several additional columns added:
* jalign_written – the number of jump alignment reads emitted for this region (i.e. with JT:i:3). This is the main result.
* jalign_low_score – the number of discarded jump alignment that had a low score (compared to the potential score of an ideal alignment)
* jalign_less_than_ref – the number of discarded jump alignments that had a score lower than a simple alignment to the reference (with some clearance).
* jump_size_min - minimal jump size expressed in *_JUMP reads. The jump size is expressed in bases (bp).
* jump_size_max - maximal ...
* jump_size_avg - averege ...
* jump_landing_min - minimal point of landing of the jump on the input region. The landing is the point at which the jump “lands” on the region. 0.0 is the start and 1.0 is the end. For example 1.2 on a 1K region means it lands 200 bp after the end of the region. 0.8 means 200 before the end.
* jump_landing_max - maximal ...
* jump_landing_avg - average ...

bam file contains reads around realignment regions. Each read can potentially appear 4 times:
* original read (with original qname)
* read aligned to left reference (around the left edge of a region). qname will be suffixed by _REF1
* read aligned to right reference (around the right edge of a region). qname will be suffixed by _REF2
* read aligned to both references, with a jump between them. qname will be suffixed by _JUMP

The four types are also marked by the JT tag, which has the value of 0, 1, 2, or 3, respectively.
Additional values for JT indicate jump alignments which were considered of low relevance or quality (see section below).

Each read (except JT=0) also contains the tag JS, which indicates the relevant alignment score.

To process whole genome it is best to shard the input bed into smaller sections (10 lines or so) and to run the tool on each smaller bed file - following with a merge stage. This can be done by WDL, bash scripts (see scripts folder) or other frameworks. 

### performance

On a 128 core machine, a complete genome run of about 10K regions takes about 30 minutes.

### additional HT values

Apart from the 0,1,2,3 value, jump alignment that were discarded are  emitted to the bam file as well with a JT that has several additive values (similar to a SAM flag, but decimal):
* 50 - jump alignment discarded because one of the 2 alignment sections is smaller than 30 bp (either on read or on ref)
* 100 - was discarded because the jump alignment score was too low compared to the ideal (all match) potential score of the read (at this point jump score has to be at least 0.85 of ideal score.
* 200 - was discarded because it scored less that a simple alignment to one of the references (without a jump)

Note that these can be added - so a JT of 250 means that the alignment had a section shorter than 30 and that it lacked in score.

