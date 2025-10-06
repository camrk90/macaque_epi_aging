#!/bin/sh

module load bedtools2-2.30.0-gcc-11.2.0

ANNO=$1
ANNO_BED=(${ANNO//./ })
BED_FILE=$2
FILE_NAME=(${BED_FILE//./ })


#This chunk converts a .txt file to .bed file if need be and runs the intersect STILL NEEDS TO BE PROPERLY CODED
if [[${ANNO} == *".txt"* ]]; then
    cat ${ANNO} > ${ANNO_BED[0]}.bed
    bedtools intersect -a ${BED_FILE} -b ${ANNO} -wo > ${FILE_NAME[0]}_intersect.txt
fi

#intersect with the lifted over chmm file
bedtools intersect -a ${BED_FILE} -b ${ANNO} -wo > ${FILE_NAME[0]}_intersect.txt

