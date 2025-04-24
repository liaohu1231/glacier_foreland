#!/bin/sh
#minimap2 -d $1.mmi $1; # make index
out=vamb_out_bam
if [ ! -d "$out" ];then mkdir $out;fi
for i in $(ls seqtk_sample/*R1.fastq.gz);do
R2=$(echo $i|sed "s/R1\.fastq/R2\.fastq/");
echo $R2
k=${i#*/};k=${k%%.*};
minimap2 --split-prefix=foo -t 64 -ax sr $1'.mmi' $i $R2 >$out/$k.sam
samtools view -@ 64 -bS $out/$k.sam > $out/$k.bam
coverm filter -b $out/$k.bam -o $out/$k"_filtered.bam" --min-read-aligned-percent 0.95 -t 64
samtools sort -@ 64 -o $out/$k.sorted $out/$k"_filtered.bam"
rm $out/$k.sam
rm $out/$k.bam
rm $out/$k"_filtered.bam"
done
jgi_summarize_bam_contig_depths --outputDepth $out/all.jgi.depth $out/*.sorted
vamb --fasta $1 --jgi $out/all.jgi.depth --outdir vamb -p 49
