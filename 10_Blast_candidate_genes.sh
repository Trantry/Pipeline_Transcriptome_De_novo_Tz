#!/bin/bash
#SBATCH --array=0#(0-(n-1) où n: nombre de job en parallèle) ###
#SBATCH --job-name=QC_complete_transcriptome_Tetraripis ###
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=150G
#SBATCH --time=1-00:00:00 ###
#SBATCH --partition=Lake
#SBATCH --output=/home/tbessonn/stdout/%A_%a.out # standard output file format
#SBATCH --error=/home/tbessonn/stderr/%A_%a.err # error file format


cd /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/9_Blast_identity/

# Make a blast database
makeblastdb -in Genome_ref/Transcriptome_whithout_isoforme.fasta -dbtype nucl -parse_seqids

# Blast the sequences
#tblastx -db tzet_transcriptome_trinity.fasta -query geisha_mogsha_genome.fa -out blast_output_geisha_mogsha.txt -outfmt "6"

# Create a file containe the IDs of the target sequences
#awk '{print$2}' blast_output_geisha_mogsha.txt > seqids_geisha_mogsha.txt

# Retrieve the sequences
blastdbcmd -db Genome_ref/Transcriptome_whithout_isoforme.fasta -dbtype nucl -entry_batch Genes_canditate_200.txt -out retrieved_sequence_Genes_candidate_200.txt

