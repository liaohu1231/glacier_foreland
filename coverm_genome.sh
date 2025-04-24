#!/bin/sh
if [ ! -d $2 ];then mkdir $2;fi
for R1 in $(ls seqtk_sample/*R1.fastq.gz);do
i=${R1#*/};j=${i%%.*}
R2=$(echo $R1|sed 's/R1\.fastq/R2\.fastq/')
echo $2'/coverm-genome.'$j'.bam'
#samtools sort -@ 20 -o $2'/coverm-genome.'$i'.sort.bam' $2'/coverm-genome.'$i'.bam'
coverm genome -d $1 -x fa -1 $R1 -2 $R2 \
--min-read-aligned-percent 0.95 -t 128 \
-m relative_abundance \
--bam-file-cache-directory $2 \
--discard-unmapped > $2/$j'.abu'
done 

