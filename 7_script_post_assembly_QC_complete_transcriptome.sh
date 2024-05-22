#!/bin/bash
#SBATCH --array=0#(0-(n-1) où n: nombre de job en parallèle) ###
#SBATCH --job-name=QC_complete_transcriptome_Tetraripis ###
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=150G
#SBATCH --time=8-00:00:00 ###
#SBATCH --partition=Lake
#SBATCH --output=/home/tbessonn/stdout/%A_%a.out # standard output file format
#SBATCH --error=/home/tbessonn/stderr/%A_%a.err # error file format

##-----------1) 
## to check that assembly is complete, we can align our reads on the assembled transcriptome.
## First, we must build a bowtie2 index
#bowtie2-build /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/2_Bowtie2/Genome_ref/trinity_tzet_complete_transcriptome.Trinity.fasta \
#    /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/2_Bowtie2/Genome_ref/index_tzet_transcriptome_trinity_ref

## Perform the alignment to capture the alignment statistics

#cd /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/2_Bowtie2/

#/usr/bin/bowtie2 \
#     -x Genome_ref/index_tzet_transcriptome_trinity_ref \
#     -1 /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/2_Bowtie2/read_use_for_transcriptom/ARNTzetteli_EKRN230034977-1A_HGK23DSX7_L4_uncontaminated_rRNA-free_reads_1.fq,/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/2_Bowtie2/read_use_for_transcriptom/Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7__sortmerna_non_rRNA_fwd.fq,/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/2_Bowtie2/read_use_for_transcriptom/Tzet_L2R1_EKRN230022263-1A_HFCCGDSX7__sortmerna_non_rRNA_fwd.fq,/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/2_Bowtie2/read_use_for_transcriptom/Tzet_L3R1_EKRN230022264-1A_HFCCGDSX7__sortmerna_non_rRNA_fwd.fq \
#     -2 /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/2_Bowtie2/read_use_for_transcriptom/ARNTzetteli_EKRN230034977-1A_HGK23DSX7_L4_uncontaminated_rRNA-free_reads_2.fq,/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/2_Bowtie2/read_use_for_transcriptom/Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7__sortmerna_non_rRNA_rev.fq,/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/2_Bowtie2/read_use_for_transcriptom/Tzet_L2R1_EKRN230022263-1A_HFCCGDSX7__sortmerna_non_rRNA_rev.fq,/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/2_Bowtie2/read_use_for_transcriptom/Tzet_L3R1_EKRN230022264-1A_HFCCGDSX7__sortmerna_non_rRNA_rev.fq \
#     -S output_sam/Tetraripis_transriptome_output.sam \
#    --no-unal \
#    --threads 24 \
#    -k 20 \
#	 2>align_stats.txt| samtools view -@ 10 -o bowtie2.bam

#conda activate BUSCO 
#cd /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/3_busco/
##-----------2)
## Let's run a BUSCO analysis on the transcriptome to check for completeness
#busco \
#	-i /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/2_Bowtie2/Genome_ref/trinity_tzet_complete_transcriptome.Trinity.fasta \
#	-l hemiptera_odb10 \
#	-m transcriptome \
#	-o busco_output \
#	--cpu 24

#conda deactivate
##-----------3)
## Blast the transcriptome on the uniprot database to do a full--length transcript analysis
## First download the data base
#cd /Xnfs/khila/abadiane/databases/uniprot_databases
#wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz

## Perform the full-length transcript analysis using blast
# Build the blastable database
#makeblastdb -in /Xnfs/khila/abadiane/databases/uniprot_databases/uniprot_sprot.fasta -dbtype prot

## Perform the blast search, reporting only the top aligments
#blastx -query /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/1_kallisto/Genome_ref/trinity_tzet_complete_transcriptome.Trinity.fasta \
#	-db /Xnfs/khila/abadiane/databases/uniprot_databases/uniprot_sprot.fasta \
#	-out /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/4_blastx/blastx.outfmt6 \
#    -evalue 1e-20 \
#	-num_threads 24 \
#	-max_target_seqs 1 \
#	-outfmt 6

## Examine the percent of the target being aligned to the best matching transcript
#/usr/lib/trinityrnaseq/util/analyze_blastPlus_topHit_coverage.pl \
#	/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/4_blastx/blastx.outfmt6 \
#	/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/1_kallisto/Genome_ref/trinity_tzet_complete_transcriptome.Trinity.fasta \
#	/Xnfs/khila/abadiane/databases/uniprot_databases/uniprot_sprot.fasta

## We can extend this process by grouping high scoring segment pairs per transcript and database hit instead of using the single best transcript
#/usr/lib/trinityrnaseq/util/misc/blast_outfmt6_group_segments.pl /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/4_blastx/blastx.outfmt6 \
#	/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/1_kallisto/Genome_ref/trinity_tzet_complete_transcriptome.Trinity.fasta \
#	/Xnfs/khila/abadiane/databases/uniprot_databases/uniprot_sprot.fasta > blastx.outfmt6.grouped
#/usr/lib/trinityrnaseq/util/misc/blast_outfmt6_group_segments.tophit_coverage.pl blastx.outfmt6.grouped > blastx.outfmt6.grouped.hist

##-----------4)
## Compute the Ex90N50 statistics
## First, we need to perform a transcript abundance estimation
#cd /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/4_blastx/
#mkdir abundance_output/

#PARAMS=(/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/real_complete_transcriptome/read_use_for_transcriptom/*fwd.fq)
#f1="${PARAMS[$SLURM_ARRAY_TASK_ID]}"
#f2=${PARAMS[$SLURM_ARRAY_TASK_ID]%%fwd.fq}"rev.fq"
#name=$(basename "${PARAMS[$SLURM_ARRAY_TASK_ID]}" "fwd.fq") 

#/usr/lib/trinityrnaseq/util/align_and_estimate_abundance.pl \
#	--transcripts /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/1_kallisto/Genome_ref/trinity_tzet_complete_transcriptome.Trinity.fasta \
#	--seqType fq \
#	--left $f1 \
#	--right $f2 \
#	--est_method kallisto \
#	--output_dir abundance_output/$name \
#	--thread_count 24 \
#	--gene_trans_map /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/1_kallisto/Genome_ref/trinity_tzet_complete_transcriptome.Trinity.fasta.gene_trans_map \
#	--prep_reference \
#	--coordsort_bam  

cd /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/4_blastx/
## Then, build transcript and gene expression matrices
/usr/lib/trinityrnaseq/util/abundance_estimates_to_matrix.pl \
	--est_method kallisto \
	--gene_trans_map /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/1_kallisto/Genome_ref/trinity_tzet_complete_transcriptome.Trinity.fasta.gene_trans_map \
	--name_sample_by_basedir abundance_output/*/*abundance.tsv





cd /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/4_blastx/
# Calculate the contig Ex90N50 statistics and Ex90 gene count
/usr/lib/trinityrnaseq/util/misc/contig_ExN50_statistic.pl \
	kallisto.isoform.counts.matrix \
	/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/1_kallisto/Genome_ref/trinity_tzet_complete_transcriptome.Trinity.fasta \
	transcript| tee ExN50.transcript.stats
/usr/lib/trinityrnaseq/util/misc/plot_ExN50_statistic.Rscript  ExN50.transcript.stats
xpdf ExN50.transcript.stats.plot.pdf

## Do the same at the gene level, which is more reliable
/usr/lib/trinityrnaseq/util/misc/contig_ExN50_statistic.pl \
         /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/4_blastx/kallisto.isoform.TMM.EXPR.matrix \
		 /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/2_Bowtie2/Genome_ref/trinity_tzet_complete_transcriptome.Trinity.fasta gene | tee ExN50.gene.stats

#followed by plotting:
/usr/lib/trinityrnaseq/util/misc/plot_ExN50_statistic.Rscript  ExN50.gene.stats
xpdf ExN50.gene.stats.plot.pdf






##Old transcriptome for comparaison
## Do the same at the gene level, which is more reliable
/usr/lib/trinityrnaseq/util/misc/contig_ExN50_statistic.pl \
         /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_novogene_illumina_X204SC23062088-Z01-F001/06_post-assembly_QC/abundance_output/abundance.tsv \
		 /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_novogene_illumina_X204SC23062088-Z01-F001/06_post-assembly_QC/tzet_transcriptome_trinity.fasta gene | tee ExN50.gene.stats

#followed by plotting:
/usr/lib/trinityrnaseq/util/misc/plot_ExN50_statistic.Rscript  ExN50.gene.stats
xpdf ExN50.gene.stats.plot.pdf
