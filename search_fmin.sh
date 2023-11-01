#!/bin/bash

#./search_fmin.sh true true nano (flipped unitigs nano)
#./search_fmin.sh false false  (unflipped eulertigs ecoli)

unitigs="unitigs"  # default unitigs
prefix="f_"  # default flipped
folder="ecoli"  # default ecoli
data="coli3682"
# unitigs or eulertigs
if [ $# -ge 1 ]; then
    if [ "$1" == "false" ]; then
        unitigs="eulertigs"
    else
        unitigs="unitigs"  # default unitigs
    fi
fi

# flipped or unflipped
if [ $# -ge 2 ]; then
    if [ "$2" == "false" ]; then
        prefix=""
    else
        prefix="f_"  # default flipped
    fi
fi

if [ $# -ge 3 ]; then
    if [ "$3" == "nano" ]; then
        folder="nanopore"
        data="SRR25689478"
    elif [ "$3" == "meta" ]; then
        folder="metagenome"
        data="ERR5035349"
    else
        folder="ecoli"  # default ecoli
        data="coli3682"
    fi
fi


for k in 21; do
    echo $k
    ./benchmark search-fmin -o ${folder}/${prefix}query-fmin-${unitigs}-${data}-k${k}_pos.txt -i ${folder}/${prefix}${data}-${unitigs}-k${k}.sbwt --lcs ${folder}/${prefix}${data}-${unitigs}-k${k}.sbwtLCS.sdsl -t 1 -f ${folder}/${prefix}${data}-${unitigs}-k${k}.sbwtFBV.sdsl -e ${folder}/${prefix}${data}-${unitigs}-k${k}.sbwtE.sdsl -g ${folder}/${prefix}${data}-${unitigs}-k${k}.sbwtO.sdsl -q /home/scratch-hdd/shared_data/finimizers/queries/data/${data}-positive.fa
    #./benchmark search-fmin -o ${folder}/${prefix}query-fmin-${unitigs}-${data}-k${k}_neg.txt -i ${folder}/${prefix}${data}-${unitigs}-k${k}.sbwt --lcs ${folder}/${prefix}${data}-${unitigs}-k${k}.sbwtLCS.sdsl -t 1 -f ${folder}/${prefix}${data}-${unitigs}-k${k}.sbwtFBV.sdsl -e ${folder}/${prefix}${data}-${unitigs}-k${k}.sbwtE.sdsl -g ${folder}/${prefix}${data}-${unitigs}-k${k}.sbwtO.sdsl -q /home/scratch-hdd/shared_data/finimizers/queries/data/random.fa
done

#k=31
#./benchmark search-fmin -o ${folder}/${prefix}query-fmin-${unitigs}-${data}-k${k}_pos.txt -i ${folder}/${prefix}${unitigs}-${data}-${k}.sbwt --lcs ${folder}/${prefix}${unitigs}-${data}-${k}.sbwtLCS.sdsl -t 1 -f ${folder}/${prefix}${unitigs}-${data}-${k}.sbwtFBV.sdsl -e ${folder}/${prefix}${unitigs}-${data}-${k}.sbwtE.sdsl -g ${folder}/${prefix}${unitigs}-${data}-${k}.sbwtO.sdsl -q /home/scratch-hdd/shared_data/finimizers/queries/data/coli3682-positive.fa
#./benchmark search-fmin -o ${folder}/${prefix}query-fmin-${unitigs}-${data}-k${k}_pos.txt -i ${folder}/${prefix}${unitigs}-${data}-${k}.sbwt --lcs ${folder}/${prefix}${unitigs}-${data}-${k}.sbwtLCS.sdsl -t 1 -f ${folder}/${prefix}${unitigs}-${data}-${k}.sbwtFBV.sdsl -e ${folder}/${prefix}${unitigs}-${data}-${k}.sbwtE.sdsl -g ${folder}/${prefix}${unitigs}-${data}-${k}.sbwtO.sdsl -q /home/scratch-hdd/shared_data/finimizers/queries/data/random.fa


