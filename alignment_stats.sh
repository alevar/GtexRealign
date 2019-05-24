#!/usr/bin/env bash

# ./alignment_stats.sh /ccb/salz4-4/gpertea/GTEX_align/by_tissue/Spleen /ccb/salz3/avaraby/GtexRealign/outSpleen 10

inputDir=${1}
outputDir=${2}
threads=${3}

touch ${outputDir}/al_stats.csv
echo name,total,old_secondary,old_al_pair,old_al_single,oldAlRate,new_secondary,new_al_pair,new_al_single,newAlRate > ${outputDir}/al_stats.csv

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
    
    # samtools flagstat --threads ${threads} ${file} > ${outputDir}/${sample}_new.tmp
    new_secondary="$(cat ${outputDir}/${sample}_new.tmp | grep 'secondary' | awk -F ' ' '{print $1}')"
    new_al_pair="$(cat ${outputDir}/${sample}_new.tmp | grep 'with itself and mate mapped' | awk -F ' ' '{print $1}')"
    new_al_single="$(cat ${outputDir}/${sample}_new.tmp | grep 'singletons' | awk -F ' ' '{print $1}')"
    newAligned=$((new_al_pair + new_al_single))
    newAlignmentRate=$(bc<<<"scale=4;${newAligned}/${oldTotal}")
    echo ${sample},${oldTotal},${old_secondary},${old_al_pair},${old_al_single},${oldAlignmentRate},${new_secondary},${new_al_pair},${new_al_single},${newAlignmentRate} >> ${outputDir}/al_stats.csv

done;
