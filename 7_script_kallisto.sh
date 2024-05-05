#!/bin/bash
#SBATCH --array=0-8#(0-(n-1) où n: nombre de job en parallèle) ###
#SBATCH --job-name=Kallisto_Tetraripis ###
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00 ###
#SBATCH --partition=E5
#SBATCH --output=/home/tbessonn/stdout/%A_%a.out # standard output file format
#SBATCH --error=/home/tbessonn/stderr/%A_%a.err # error file format

#create index pour kallisto
#kallisto index -i /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/7_kallisto/Genome_ref/tz_transcripts.idx /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/7_kallisto/Genome_ref/tzet_transcriptome_trinityV2_kraken2.Trinity.fasta



cd /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/7_kallisto/

PARAMS=(/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/7_kallisto/unclassified_sormeRNA_after_kraken2/*fwd.fq)
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


