#!/bin/bash
#SBATCH --array=0-8 #(0-(n-1) où n: nombre de job en parallèle) ###
#SBATCH --job-name=Kraken2_Tetraripis ###
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=0-12:00:00 ###
#SBATCH --partition=E5
#SBATCH --output=/home/tbessonn/stdout/%A_%a.out # standard output file format
#SBATCH --error=/home/tbessonn/stderr/%A_%a.err # error file format

module use /applis/PSMN/debian11/E5/modules/all
 
cd /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/6_fastqc


PARAMS=(/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/6_fastqc/SortmeRNA_nonrRNA/*fwd.fq) ####
f1="${PARAMS[$SLURM_ARRAY_TASK_ID]}"
f2=${PARAMS[$SLURM_ARRAY_TASK_ID]%%fwd.fq}"rev.fq"
name=$(basename "${PARAMS[$SLURM_ARRAY_TASK_ID]}" "fwd.fq") 

fastqc \
	$f1 \
 	$f2 \
	-o fastqc2nd/
