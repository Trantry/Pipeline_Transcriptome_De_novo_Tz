#!/bin/bash
#SBATCH --array=0 #(0-(n-1) où n: nombre de job en parallèle) ###
#SBATCH --job-name=Trinity_complete_transcriptome_Tetraripis ###
#SBATCH --ntasks=1
#SBATCH --mem=150G
#SBATCH --cpus-per-task=24
#SBATCH --time=8-00:00:00 ###
#SBATCH --partition=Lake
#SBATCH --output=/home/tbessonn/stdout/%A_%a.out # standard output file format
#SBATCH --error=/home/tbessonn/stderr/%A_%a.err # error file format
 
cd /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/6_trinity/complete_transcriptome/

## De novo transcriptome assemble using Trinity.
Trinity \
    --seqType fq \
	--left /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/6_trinity/complete_transcriptome/read_use_for_transcriptom/ARNTzetteli_EKRN230034977-1A_HGK23DSX7_L4_uncontaminated_rRNA-free_reads_1.fq,/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/6_trinity/complete_transcriptome/read_use_for_transcriptom/Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7__sortmerna_non_rRNA_fwd.fq,/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/6_trinity/complete_transcriptome/read_use_for_transcriptom/Tzet_L2R1_EKRN230022263-1A_HFCCGDSX7__sortmerna_non_rRNA_fwd.fq,/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/6_trinity/complete_transcriptome/read_use_for_transcriptom/Tzet_L3R1_EKRN230022264-1A_HFCCGDSX7__sortmerna_non_rRNA_fwd.fq \
	--right /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/6_trinity/complete_transcriptome/read_use_for_transcriptom/ARNTzetteli_EKRN230034977-1A_HGK23DSX7_L4_uncontaminated_rRNA-free_reads_2.fq,/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/6_trinity/complete_transcriptome/read_use_for_transcriptom/Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7__sortmerna_non_rRNA_rev.fq,/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/6_trinity/complete_transcriptome/read_use_for_transcriptom/Tzet_L2R1_EKRN230022263-1A_HFCCGDSX7__sortmerna_non_rRNA_rev.fq,/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/6_trinity/complete_transcriptome/read_use_for_transcriptom/Tzet_L3R1_EKRN230022264-1A_HFCCGDSX7__sortmerna_non_rRNA_rev.fq \
	--CPU 24 \
	--max_memory 150G \
	--output /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/6_trinity/complete_transcriptome/read_use_for_transcriptom/trinity_tzet_complete_transcriptome \
	--full_cleanup

## Create a file containing the stats of the assembly
/usr/lib/trinityrnaseq/util/TrinityStats.pl trinity_tzet_complete_transcriptome.fasta > trinity_tzet_complete_transcriptome.log

