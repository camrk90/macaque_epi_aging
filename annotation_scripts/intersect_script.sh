#!/bin/sh

module load bedtools2-2.30.0-gcc-11.2.0

#the glm bed file to be intersected
CPGS=regions_to_cpgs.txt
Y_CPGS=y_regions_to_cpgs.txt
ALL=all_regions_to_cpgs.txt

#Intersect autosome/x-chrom cpgs and regions
bedtools intersect -a cayo_cpgs2.bed -b regions.bed -wa -wb > ${CPGS}

#Intersect y-chrom cpgs and regions
bedtools intersect -a y_cpgs2.bed -b regions_y.bed -wa -wb > ${Y_CPGS}

#Concatenate autosomes/x-chrom/y-chrom cpgs
cat ${CPGS} ${Y_CPGS} > ${ALL}

#intersect with the repeats file
bedtools intersect -a repeats.bed -b ${ALL} -wo > region_repeat_intersect.txt

#intersect with the lifted over chmm file
bedtools intersect -a chmm_mmul.bed -b ${ALL} -wo > region_chmm_intersect.txt
