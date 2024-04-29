#!/bin/bash
#SBATCH --array=0 #(0-(n-1) où n: nombre de job en parallèle) ###
#SBATCH --job-name=TrinityV2_kraken2_Tetraripis ###
#SBATCH --ntasks=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=16
#SBATCH --time=5-00:00:00 ###
#SBATCH --partition=Lake
#SBATCH --output=/home/tbessonn/stdout/%A_%a.out # standard output file format
#SBATCH --error=/home/tbessonn/stderr/%A_%a.err # error file format

cd /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/7_trinity/read_trinity_sortmerna_kraken2/

## De novo transcriptome assemble using Trinity.
Trinity \
    --seqType fq \
	--left /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/7_trinity/read_trinity_sortmerna_kraken2/Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7__sortmerna_non_rRNA_fwd.fq,/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/7_trinity/read_trinity_sortmerna_kraken2/Tzet_L2R1_EKRN230022263-1A_HFCCGDSX7__sortmerna_non_rRNA_fwd.fq,/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/7_trinity/read_trinity_sortmerna_kraken2/Tzet_L3R1_EKRN230022264-1A_HFCCGDSX7__sortmerna_non_rRNA_fwd.fq \
	--right /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/7_trinity/read_trinity_sortmerna_kraken2/Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7__sortmerna_non_rRNA_rev.fq,/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/7_trinity/read_trinity_sortmerna_kraken2/Tzet_L2R1_EKRN230022263-1A_HFCCGDSX7__sortmerna_non_rRNA_rev.fq,/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/7_trinity/read_trinity_sortmerna_kraken2/Tzet_L3R1_EKRN230022264-1A_HFCCGDSX7__sortmerna_non_rRNA_rev.fq \
	--CPU 16 \
	--max_memory 120G \
	--output /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/7_trinity/read_trinity_sortmerna_kraken2/tzet_transcriptome_trinityV2_kraken2 \
	--full_cleanup

## Create a file containing the stats of the assembly
/usr/lib/trinityrnaseq/util/TrinityStats.pl tzet_transcriptome_trinityV2_kraken2.fasta > tzet_trinity_statsV2_kraken2.log
