#!/bin/sh
source activate iphop_env
iphop predict --fa_file $1 --db_dir $2 --out_dir $3 --num_threads 90
