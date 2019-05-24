#!/usr/bin/env bash

# ./alignment_stats.sh /ccb/salz4-4/gpertea/GTEX_align/by_tissue/Spleen /ccb/salz3/avaraby/GtexRealign/outSpleen 10

inputDir=${1}
outputDir=${2}
threads=${3}

touch ${outputDir}/al_stats.csv
echo name,total,old_secondary,oldAligned,oldAlRate,new_secondary,newAligned,newAlRate > ${outputDir}/al_stats_py.csv

for file in ${outputDir}/*bam ; do
    baseName=$(basename ${file})
    sample=${baseName%.*}
    echo ${sample}

    # samtools flagstat --threads ${threads} ${inputDir}/${sample}.cram > ${outputDir}/${sample}.tmp
    oldTotal_tmp="$(cat ${outputDir}/${sample}.tmp | grep 'in total (QC-passed' | awk -F ' ' '{print $1}')"
    old_secondary="$(cat ${outputDir}/${sample}.tmp | grep 'secondary' | awk -F ' ' '{print $1}')"
    oldTotal=$((oldTotal_tmp - old_secondary))
    old_al_pair="$(cat ${outputDir}/${sample}.tmp | grep 'with itself and mate mapped' | awk -F ' ' '{print $1}')"
    old_al_single="$(cat ${outputDir}/${sample}.tmp | grep 'singletons' | awk -F ' ' '{print $1}')"
    oldAligned=$((old_al_pair + old_al_single))
    oldAlignmentRate=$(bc<<<"scale=4;${oldAligned}/${oldTotal}")

    ./compute_alrate.py -i ${file} > ${outputDir}/${sample}.tmp2
    new_secondary="$(grep 'secondary' ${outputDir}/${sample}.tmp2 | awk -F '\t' '{print $2}')"
    newAligned="$(grep 'unique' ${outputDir}/${sample}.tmp2 | awk -F '\t' '{print $2}')"
    newAlignmentRate=$(bc<<<"scale=4;${newAligned}/${oldTotal}")
    
    echo ${sample},${oldTotal},${old_secondary},${oldAligned},${oldAlignmentRate},${new_secondary},${newAligned},${newAlignmentRate} >> ${outputDir}/al_stats_py.csv

done;
