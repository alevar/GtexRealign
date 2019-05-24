#!/usr/bin/env bash

# ./readlen_stats.sh /ccb/salz4-4/gpertea/GTEX_align/by_tissue/ ./realen.stats 32

inputDir=${1}
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

for file in ${inputDir}/*/*cram ; do
    ((i=i%numThreads)); ((i++==0)) && wait
        runTask ${file} &
done;
