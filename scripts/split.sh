rm -rf out
mkdir -p out/beds
cat ../../cnmops.win500.NA24385.cnvs.with_label.bed \
	| grep -E "(CN0)|(CN1)" \
	| split -l 10 -d --additional-suffix .bed
mv x*.bed out/beds
