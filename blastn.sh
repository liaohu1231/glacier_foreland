#!/bin/sh
#blastn -query $1 -db $2 -evalue 1 -task blastn-short -outfmt 6 -num_threads 64 -out $3
makeblastdb -in $1 -dbtype nucl -out $1
blastn -query $2 -db $1 -evalue 1 -outfmt "6 qseqid sseqid pident qcovs length evalue qlen slen qstart qend sstart send bitscore" \
-num_threads 64 -out $3
