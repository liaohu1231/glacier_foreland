#!/bin/sh
gtdbtk classify_wf --genome_dir $1 \
--extension fa --cpus 64 --pplacer_cpus 64 \
--out_dir $2
