#!/bin/bash
#SBATCH --array=0-8#(0-(n-1) où n: nombre de job en parallèle) ###
#SBATCH --job-name=Kallisto_Tetraripis ###
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --time=0-01:00:00 ###
#SBATCH --partition=Lake
#SBATCH --output=/home/tbessonn/stdout/%A_%a.out # standard output file format
#SBATCH --error=/home/tbessonn/stderr/%A_%a.err # error file format

cd /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/1_kallisto/

#create index pour kallisto
kallisto index -i /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/1_kallisto/Genome_ref/tz_transcripts.idx /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/1_kallisto/Genome_ref/trinity_tzet_complete_transcriptome.Trinity.fasta

cd /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/1_kallisto/

PARAMS=(/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/1_kallisto/unclassified_sormeRNA_after_kraken2/*fwd.fq)
f1="${PARAMS[$SLURM_ARRAY_TASK_ID]}"
f2=${PARAMS[$SLURM_ARRAY_TASK_ID]%%fwd.fq}"rev.fq"
name=$(basename "${PARAMS[$SLURM_ARRAY_TASK_ID]}" "fwd.fq") 


kallisto quant \
     -i Genome_ref/tz_transcripts.idx \
     -o kallisto_output/$name \
     -t 24 \
     $f1 \
     $f2 \
     -b 100 

