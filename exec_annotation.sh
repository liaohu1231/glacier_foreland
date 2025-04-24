#!/bin/sh
source activate plasflow
out=${1%.*}
exec_annotation -k /share/database/kfamscan/ko_list \
-p /share/database/kfamscan/profiles -f detail-tsv -o $out'_exec.tsv' --cpu 64 $1
cat $out'_exec.tsv'|sed -n '/^\*/p' > $out'_sig_exec.tsv'
