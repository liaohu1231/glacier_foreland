#!/bin/sh
#source activate metawrap
for i in $(ls */*1500.fasta);do j=${i%/*};k=${j%_*}
if [ ! -d "$j"/refinement ];then
bowtie2-build $i bowtie2_database/$k
R1=$(ls seqtk_sample/${k}*.R1.4800.fastq.gz)
R2=$(ls seqtk_sample/${k}*.R2.4800.fastq.gz)
bowtie2 -p 64 --non-deterministic -x bowtie2_database/$k --very-sensitive -1 ${R1} -2 ${R2} -S $j/$k.sam
samtools view -@ 64 -bS $j/$k.sam > $j/$k.bam
coverm filter -b $j/$k.bam -o $j/$k"_filtered.bam" --min-read-aligned-percent 0.99 -t 64
samtools sort -@ 64 -o $j/$k.sorted $j/$k"_filtered.bam"
metawrap binning -a $i -l 1500 -o $j -m 512 --metabat1 --metabat2 --concoct --interleaved -t 30 $j/$k'.sorted'
metawrap bin_refinement -o $j/refinement -A $j/concoct_bins \
-B $j/metabat1_bins -C $j/metabat2_bins -c 50 -t 64
fi
done

#source activate gtdbtk
#gtdbtk classify_wf --genome_dir bins_refinement/ --extension fa --cpus 40 --pplacer_cpus 40 --out_dir gtdbtk_bins
