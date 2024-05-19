#!/bin/bash
#SBATCH --array=0#(0-(n-1) où n: nombre de job en parallèle) ###
#SBATCH --job-name=QC_complete_transcriptome_Tetraripis ###
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=150G
#SBATCH --time=8-00:00:00 ###
#SBATCH --partition=Lake
#SBATCH --output=/home/tbessonn/stdout/%A_%a.out # standard output file format
#SBATCH --error=/home/tbessonn/stderr/%A_%a.err # error file format


#cd /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/9_Blast_identity/

# Make a blast database
#makeblastdb -in Genome_ref/Transcriptome_whithout_isoforme.fasta -dbtype nucl -parse_seqids

# Blast the sequences
#tblastx -db tzet_transcriptome_trinity.fasta -query geisha_mogsha_genome.fa -out blast_output_geisha_mogsha.txt -outfmt "6"

# Create a file containe the IDs of the target sequences
#awk '{print$2}' blast_output_geisha_mogsha.txt > seqids_geisha_mogsha.txt

# Retrieve the sequences 3743
#blastdbcmd -db Genome_ref/Transcriptome_whithout_isoforme.fasta -dbtype nucl -entry_batch Candidate_shrinkLFC_genes_3743_pvalue_0.01.txt  -out retrieved_sequence_Candidate_shrinkLFC_genes_3743_pvalue_0.01.fasta

# Retrieve the sequences for the 926 genes overxpressed
#blastdbcmd -db Genome_ref/Transcriptome_whithout_isoforme.fasta -dbtype nucl -entry_batch overexpressed_shrink_LFC_gene_926_ids.txt -out retrieved_sequence_overexpressed_shrink_LFC_gene_926_ids.fasta

#Same for underexpressed th 2817 genes
#blastdbcmd -db Genome_ref/Transcriptome_whithout_isoforme.fasta -dbtype nucl -entry_batch underexpressed_shrink_LFC_gene_2817_ids.txt -out retrieved_sequence_underexpressed_shrink_LFC_gene_2817_ids.fasta


##########################################################################################################
#Homologous genes in Rhagovelia as present in the 3743 differntial gens expressed in Tetraripis
#cd /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/9_Blast_identity/Rhago_homologous_Tetra/

# Make a blast database
#makeblastdb -in retrieved_sequence_Candidate_shrinkLFC_genes_3743_pvalue_0.01.fasta -dbtype nucl -parse_seqids

# Blast the sequences
#tblastx -db retrieved_sequence_Candidate_shrinkLFC_genes_3743_pvalue_0.01.fasta  -query Protein_homologous_Rhagovelia.fasta -out blast_output_Protein_homologous_Rhagovelia.txt -outfmt "6" -evalue 1e-20

# Create a file containe the IDs of the target sequences
#awk '{print$2}' blast_output_Protein_homologous_Rhagovelia.txt > seqids_Protein_homologous_Rhagovelia.txt

# Retrieve the sequences 
#blastdbcmd -db retrieved_sequence_Candidate_shrinkLFC_genes_3743_pvalue_0.01.fasta -dbtype nucl -entry_batch seqids_Protein_homologous_Rhagovelia.txt -out retrieved_sequence_seqids_Protein_homologous_Rhagovelia.txt
 



##########################################################################################################
#Blast whith NCBI data Protein nr Arthropoda present in the 3743 differntial gens expressed in Tetraripis
cd /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/9_Blast_identity/Blast_NCBI_protein_nr/
# Make a blast database
makeblastdb -in /Xnfs/khila/database/BLAST/ncbi_database/ncbi_nr_protein_database/nr -dbtype prot -parse_seqids

#####################################################################################
# Blast the sequences 926 overexpressed genes
blastx -db /Xnfs/khila/database/BLAST/ncbi_database/ncbi_nr_protein_database/nr  -query retrieved_sequence_overexpressed_shrink_LFC_gene_926_ids.fasta -out blast_output_sequence_overexpressed_shrink_LFC_gene_926_ids.txt -outfmt "6" -evalue 1e-20 

# Create a file containe the IDs of the target sequences
awk '{print$2}' blast_output_sequence_overexpressed_shrink_LFC_gene_926_ids.txt > seqids_Protein_blast_output_sequence_overexpressed_shrink_LFC_gene_926_ids.txt

# Retrieve the sequences 
blastdbcmd -db /Xnfs/khila/database/BLAST/ncbi_database/ncbi_nr_protein_database/nr -entry_batch seqids_Protein_blast_output_sequence_overexpressed_shrink_LFC_gene_926_ids.txt -out retrieved_sequence_seqids_Protein_overexpressed_shrink_LFC_gene_926_ids.fasta
 
#####################################################################################
# Blast the sequences 2817 underexpressed genes
blastx -db /Xnfs/khila/database/BLAST/ncbi_database/ncbi_nr_protein_database/nr  -query retrieved_sequence_underexpressed_shrink_LFC_gene_2817_ids.fasta -out blast_output_retrieved_sequence_underexpressed_shrink_LFC_gene_2817_ids.txt -outfmt "6" -evalue 1e-20 

# Create a file containe the IDs of the target sequences
awk '{print$2}' blast_output_retrieved_sequence_underexpressed_shrink_LFC_gene_2817_ids.txt > seqids_Protein_blast_output_sequence_underexpressed_shrink_LFC_gene_2817_ids.txt

# Retrieve the sequences 
blastdbcmd -db /Xnfs/khila/database/BLAST/ncbi_database/ncbi_nr_protein_database/nr -entry_batch seqids_Protein_blast_output_sequence_underexpressed_shrink_LFC_gene_2817_ids.txt -out retrieved_sequence_seqids_Protein_seqids_Protein_underexpressed_shrink_LFC_gene_2817_ids.fasta
 
