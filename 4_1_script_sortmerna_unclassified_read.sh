#!/bin/bash
#SBATCH --array=0-8 #(0-(n-1) où n: nombre de job en parallèle) ###
#SBATCH --job-name=Sortmerna_Tetraripis_unclassified ###
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --cpus-per-task=24
#SBATCH --time=3-00:00:00 ###
#SBATCH --partition=E5
#SBATCH --output=/home/tbessonn/stdout/%A_%a.out # standard output file format
#SBATCH --error=/home/tbessonn/stderr/%A_%a.err # error file format

module use /applis/PSMN/debian11/E5/modules/all

cd /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/3_kraken2 

#Number of task= 9
## Run sortmerna for unclassified
PARAMS=(/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/3_kraken2/unclassified/*1classified_reads.fq)
f1="${PARAMS[$SLURM_ARRAY_TASK_ID]}"
f2=${PARAMS[$SLURM_ARRAY_TASK_ID]%%1classified_reads.fq}"2classified_reads.fq"
name=$(basename "${PARAMS[$SLURM_ARRAY_TASK_ID]}" "1classified_reads.fq") 

/Xnfs/khila/software/sortmerna-4.3.6-Linux/bin/sortmerna \
	--ref /Xnfs/khila/abadiane/databases/sortmerna_databases/rfam-5.8s-database-id98.fasta \
	--ref /Xnfs/khila/abadiane/databases/sortmerna_databases/rfam-5s-database-id98.fasta \
	--ref /Xnfs/khila/abadiane/databases/sortmerna_databases/silva-euk-18s-id95.fasta \
	--ref /Xnfs/khila/abadiane/databases/sortmerna_databases/silva-euk-28s-id98.fasta \
	--reads $f1 \
	--reads $f2 \
  	--workdir /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/3_kraken2/sortmeRNA_unclassified/ \
	--fastx  \
	--paired_out \
	--out2 \
	--zip-out \
	--idx-dir /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/3_kraken2/sortmeRNA_unclassified/idx/$name"idx" \
	--kvdb /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/3_kraken2/sortmeRNA_unclassified/kvdb/$name"kvdb" \
	--readb /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/3_kraken2/sortmeRNA_unclassified/readb/$name"readb" \
	--aligned /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/3_kraken2/sortmeRNA_unclassified/$name"_sortmerna_rRNA" \
	--other /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/3_kraken2/sortmeRNA_unclassified/$name"_sortmerna_non_rRNA"




