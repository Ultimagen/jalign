ls -1 out/beds/x*.bed \
	| tr "." "/" \
	| cut -d '/' -f 3 \
	| awk '
		{
			printf("sudo docker run --rm -v .:/src -v ~/src/hdata:/ref -v ~/data:/data jump_align_dev ./cnv_realign.py /data/jump_align/NA24385-Z0027.cram /data/jump_align/runs/r1/out/beds/%s.bed /ref/Homo_sapiens_assembly38.fasta /data/jump_align/runs/r1/out/%s\n", $1, $1)
		}
		'
