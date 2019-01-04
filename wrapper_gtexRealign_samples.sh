#!/usr/bin/env bash

# realign all samples of a tissue from CRAM format

# ./wrapper_gtexRealign_samples.sh /ccb/salz4-4/gpertea/GTEX_align/by_tissue/Spleen ./test /ccb/salz3/avaraby/genomes/human/hg38/hg38_p12.fa asdf asdf 30

inputDir=${1}
outputDir=${2}
ref=${3}
transIDX=${4}
genomeIDX=${5}
threads=${6}

mkdir -p ${outputDir}
mkdir -p ${outputDir}/tmpData

for file in ${inputDir}/*cram ; do
    baseName=$(basename ${file})
    sample=${baseName%.*}
    echo ${sample}

    samtools bam2fq -1 ${outputDir}/tmpData/${sample}_r1.fq -2 ${outputDir}/tmpData/${sample}_r2.fq -s - -0 /dev/null -O -F 0x900 --reference ${ref} --threads ${threads} ${file} | awk 'BEGIN{RS="@";FS="\n";OFS="";ORS=""}; $1 in f {print "@"f[$1] > "'"${outputDir}/tmpData/${sample}_s1.fq"'"; print "@"$0 > "'"${outputDir}/tmpData/${sample}_s2.fq"'"} {f[$1]=$0}'

    transHISAT2.py align --m1 ${outputDir}/tmpData/${sample}_r1.fq,${outputDir}/tmpData/${sample}_s1.fq --m2 ${outputDir}/tmpData/${sample}_r2.fq,${outputDir}/tmpData/${sample}_s2.fq --db ${transIDX} -o ${outputDir}/${sample}.bam --type hisat --genome-db ${genomeIDX} --threads ${threads}

done;

rm -rf ${outputDir}/tmpData
