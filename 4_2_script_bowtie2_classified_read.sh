#!/bin/bash
#SBATCH --array=0-8#(0-(n-1) où n: nombre de job en parallèle) ###
#SBATCH --job-name=Bowtie2_Tetraripis ###
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=0-12:00:00 ###
#SBATCH --partition=E5
#SBATCH --output=/home/tbessonn/stdout/%A_%a.out # standard output file format
#SBATCH --error=/home/tbessonn/stderr/%A_%a.err # error file format

module use /applis/PSMN/debian11/E5/modules/all
 


#bowtie2-build /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/4_Bowtie2/transcriptome_ref/tzet_transcriptome_trinity.fasta \
#    /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/4_Bowtie2/transcriptome_ref/index/index_tzet_transcriptome_ref



cd /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/4_Bowtie2/

PARAMS=(/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/4_Bowtie2/classified/*1classified_reads.fq)
f1="${PARAMS[$SLURM_ARRAY_TASK_ID]}"
f2=${PARAMS[$SLURM_ARRAY_TASK_ID]%%1classified_reads.fq}"2classified_reads.fq"
name=$(basename "${PARAMS[$SLURM_ARRAY_TASK_ID]}" "1classified_reads.fq") 


/usr/bin/bowtie2 \
     -x transcriptome_ref/index/index_tzet_transcriptome_ref \
     -1 $f1 \
     -2 $f2 \
     -S output_sam/$name"_output.sam" \
    --no-unal \
    --threads 16
