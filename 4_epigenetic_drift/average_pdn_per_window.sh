#!/bin/bash
#
#SBATCH --mail-user=ckelsey4@asu.edu
#SBATCH --mail-type=FAIL

module load bedtools2-2.30.0-gcc-11.2.0

#Map PDN per id to 200bp windows
lidpid=`sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/nsnyderm/cayo_rrbs/sort_bam/lid_pid`

bedtools map -F 0.51 -o count,mean -a /scratch/nsnyderm/cayo_rrbs/discordance/allregions_reference.bed \
        -b /scratch/ckelsey4/Cayo_meth/epigenetic_drift/pdn_out/${lidpid}.PDN_per_read.bed > /scratch/ckelsey4/Cayo_meth/epigenetic_drift/avg_pdn/${lidpid}_meanPDN.bed

echo Mapping ${lid_pid} complete
