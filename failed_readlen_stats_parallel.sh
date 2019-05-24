#!/usr/bin/env bash

# ./failed_readlen_stats_parallel.sh failed.readlen.lst ./realen.stats 32

inputFP=${1}
outFP=${2}
numThreads=${3}

touch ${outFP}
echo fp,readlen > ${outFP}

runTask(){
    echo ${file}
    readlen="$(samtools view ${file} | cut -f10 | awk '{print length}' | sort | uniq -c | awk '{$1=$1;print}' | tr '\n' ';' | tr ' ' ':')"
    readlen=${readlen%?}

    echo ${file},${readlen} >> ${outFP}
}   

# for file in ${inputDir}/*/*cram ; do
    # ((i=i%numThreads)); ((i++==0)) && wait
        # runTask ${file} &
# done;

while read file; do
    ((i=i%numThreads)); ((i++==0)) && wait
        runTask ${file} &
done < ${inputFP}
