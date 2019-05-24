#!/usr/bin/env bash

inputDir=${1}
outputDir=${2}
ref=${3}
ann=${4}
threadsPerTissue=${5:-1}
threadsPerSample=${6:-1}
transIDX=${7}
genomeIDX=${8}

# ./wrapper_gtexRealign_tissues.sh /ccb/salz4-4/gpertea/GTEX_align/by_tissue/ test /ccb/salz3/avaraby/genomes/human/hg38/hg38_p12.fa data/hg38_p12_gencode29_refseq.gff3 4 10 data/transHisatIDX /ccb/salz3/avaraby/genomes/human/hg38/hg38_p12

# check that the reference genome exists
if [ ! ${ref} ]; then
    echo "REFERENCE GENOME DOES NOT EXIST. PLEASE VERIFY THE PATH"
    exit 1
fi

# check that the annotation of the reference genome exists
if [ ! ${ann} ]; then
    echo "ANNOTATION OF THE GENOME DOES NOT EXIST. PLEASE VERIFY THE PATH"
    exit 1
fi

# first create directories
mkdir -p ${outputDir}

# next check if the genome index exists and if it doesnt - ask and create
if [ -z "$genomeIDX"  ]; then
    # lets see if an index exists
    refBase=${ref%.*}
    if [ ! -f ${ref} ]; then
        echo "REFERENCE FILE DOES NOT EXIST"
        exit 1
    else
        echo "SEARCHING FOR A HISAT GENOME INDEX"
        if [ ! -f ${refBase}.1.ht2 ]; then
            echo "HISAT GENOME INDEX NOT FOUND"
            read -p "Do you want to create a hisat2 genome index for the provided reference genome?" yn
            case $yn in
                [Yy]* ) hisat2-build -p ${threadsPerSample} ${ref} ${refBase} ; break;;
                [Nn]* ) exit 1;;
                * ) echo "Please answer yes or no.";;
            esac
            genomeIDX="${refBase}"
        else
            echo "GENOME INDEX IS FOUND"
        fi
    fi
else
    echo "GENOME INDEX IS PROVIDED"
fi;

# next check if the transIDX exists, and if doesnt - ask and create
if [ -z "$transIDX"  ]; then
    if [ ! -f ${ann} ]; then
        echo "ANNOTATION NOT FOUND"
        exit 1
    fi
    if [ ! -f ${ref} ]; then
        echo "REFERENCE FILE NOT FOUND"
        exit 1
    fi

    echo "transHISAT2 INDEX NOT FOUND"
    read -p "Do you want to create a transHISAT2 index for the provided reference genome and annotation using hisat internally?" yn
    case $yn in
        [Yy]* ) transHISAT2.py build --gff ${ann} --ref ${ref} -o ${outputDir}/transHisat2IDX --threads ${threadsPerSample} --type hisat; break;;
        [Nn]* ) exit 1;;
        * ) echo "Please answer yes or no.";;
    esac
    transIDX="${outputDir}/transHisat2IDX"
else
    echo "TRANSCRIPTOME INDEX IS PROVIDED"
fi;

# lastly iterate over all tissues and run the sub-wrapper
for dirFP in ${inputDir}/*/ ; do
    curTissue=$(basename ${dirFP})
    echo ${dirFP} ${curTissue_tmp}
    ./wrapper_gtexRealign_samples.sh ${dirFP} ${outputDir}/${curTissue} ${ref} ${genomeIDX} ${transIDX} ${threadsPerSample} ${ann}
done;
