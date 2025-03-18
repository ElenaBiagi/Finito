#!/bin/bash

for folder in ecoli nanopore metagenome; do
    if [ $folder == "ecoli" ]; then
        data="coli3682"
    elif [ $folder == "metagenome" ]; then
        data="ERR5035349"
    elif [ $folder == "nanopore" ]; then
        data="SRR25689478"
    fi
    for k in 127 ; do
        echo $k
        ./benchmark search-fmin -i /home/scratch-hdd/ebiagi/Finimizers/double/${folder}/double_f_${data}-unitigs-k${k} -q /home/scratch-hdd/shared_data/finimizers/queries/data/${data}-positive.fa -o ./pos_results_${data}_${k}
        ./benchmark search-fmin -i /home/scratch-hdd/ebiagi/Finimizers/double/${folder}/double_f_${data}-unitigs-k${k} -q /home/scratch-hdd/shared_data/finimizers/queries/data/random.fa  -o ./pos_results_${data}_${k}
    done
done