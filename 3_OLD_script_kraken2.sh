cat "fichier1" "Fichier2" > "nouveau_nom". #concatenate

wc -l "nom_fichier" #compte les lignes/nombre de read

#!/bin/bash
#
#SBATCH --job-name=kraken2_Tetraripis_zetteli
#SBATCH --ntasks=24
#SBATCH --mem=100G
#SBATCH --time=0-24:00:00
#SBATCH --partition=E5
#SBATCH --output=%x.out
#SBATCH --error=%x.err

#concatnenate des version du L1R1-fwd
#cat Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_L1_1_fwd-paired.fq.gz \
#    Tzet_L1R1_EKRN230022262-1A_HFC7MDSX7_L3_1_fwd-paired.fq.gz \
#    Tzet_L1R1_EKRN230022262-1A_HFCCLDSX7_L2_1_fwd-paired.fq.gz > Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_concatenated_fwd.fq.gz #concatenate

#concatnenate des version du L1R1-rev
#cat Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_L1_2_rev-paired.fq.gz \
#    Tzet_L1R1_EKRN230022262-1A_HFC7MDSX7_L3_2_rev-paired.fq.gz \
#    Tzet_L1R1_EKRN230022262-1A_HFCCLDSX7_L2_2_rev-paired.fq.gz > Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_concatenated_rev.fq.gz #concatenate

cd //home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/kraken2/

## Run Kraken2 for leg1R1
kraken2 --db standard_database \
	--paired samples/Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_concatenated_fwd.fq.gz samples/Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_concatenated_rev.fq.gz \
	--threads 24 \
	--classified-out leg1_R1_classified_reads#.txt \
	--unclassified-out leg1_R1_unclassified_reads#.txt \
	--report leg1_R1_kraken2_report.txt \
	--gzip-compressed 
	
## Run Kraken2 for leg1R2
#kraken2 --db standard_database \
	--paired Tzet_L1R2_EKRN230022265-1A_HF3TCDSX7_L1_1_fwd-paired.fq.gz Tzet_L1R2_EKRN230022265-1A_HF3TCDSX7_L1_2_rev-paired.fq.gz \
	--threads 24 \
	--classified-out leg1_R2_classified_reads#.txt \
	--unclassified-out leg1_R2_unclassified_reads#.txt \
	--report leg1_R2_kraken2_report.txt \
	--gzip-compressed 

## Run Kraken2 for leg1R3
kraken2 --db standard_database \
	--paired Tzet_L1R3_EKRN230022268-1A_HF3TCDSX7_L1_1_fwd-paired.fq.gz Tzet_L1R3_EKRN230022268-1A_HF3TCDSX7_L1_2_rev-paired.fq.gz \
	--threads 24 \
	--classified-out leg1_R3_classified_reads#.txt \
	--unclassified-out leg1_R3_unclassified_reads#.txt \
	--report leg1_R3_kraken2_report.txt \
	--gzip-compressed 

#concatnenate des version du L2R1-fwd
#cat Tzet_L2R1_EKRN230022263-1A_HF3NCDSX7_L3_1_fwd-paired.fq.gz \
#    Tzet_L2R1_EKRN230022263-1A_HFCCGDSX7_L4_1_fwd-paired.fq.gz > Tzet_L2R1_EKRN230022263-1A_HFCCGDSX7_concatenated_fwd.fq.gz #concatenate

#concatnenate des version du L2R1-fwd
#cat Tzet_L2R1_EKRN230022263-1A_HF3NCDSX7_L3_2_rev-paired.fq.gz \
#    Tzet_L2R1_EKRN230022263-1A_HFCCGDSX7_L4_2_rev-paired.fq.gz > Tzet_L2R1_EKRN230022263-1A_HFCCGDSX7_concatenated_rev.fq.gz #concatenate

## Run Kraken2 for leg2R1
kraken2 --db standard_database \
	--paired Tzet_L2R1_EKRN230022263-1A_HFCCGDSX7_concatenated_fwd.fq.gz Tzet_L2R1_EKRN230022263-1A_HFCCGDSX7_concatenated_rev.fq.gz \
	--threads 24 \
	--classified-out leg2_R1_classified_reads#.txt \
	--unclassified-out leg2_R1_unclassified_reads#.txt \
	--report leg2_R1_kraken2_report.txt \
	--gzip-compressed

## Run Kraken2 for leg2R2
kraken2 --db standard_database \
	--paired Tzet_L2R2_EKRN230022266-1A_HF3NCDSX7_L3_1_fwd-paired.fq.gz Tzet_L2R2_EKRN230022266-1A_HF3NCDSX7_L3_2_rev-paired.fq.gz \
	--threads 24 \
	--classified-out leg2_R2_classified_reads#.txt \
	--unclassified-out leg2_R2_unclassified_reads#.txt \
	--report leg2_R2_kraken2_report.txt \
	--gzip-compressed

## Run Kraken2 for leg2R3
kraken2 --db standard_database \
	--paired Tzet_L2R3_EKRN230022269-1A_HF3NCDSX7_L3_1_fwd-paired.fq.gz Tzet_L2R3_EKRN230022269-1A_HF3NCDSX7_L3_2_rev-paired.fq.gz \
	--threads 24 \
	--classified-out leg2_R3_classified_reads#.txt \
	--unclassified-out leg2_R3_unclassified_reads#.txt \
	--report leg2_R3_kraken2_report.txt \
	--gzip-compressed

#concatnenate des version du L3R1-fwd
#cat Tzet_L3R1_EKRN230022264-1A_HF3NCDSX7_L3_1_fwd-paired.fq.gz \
#    Tzet_L3R1_EKRN230022264-1A_HFCCGDSX7_L4_1_fwd-paired.fq.gz > Tzet_L3R1_EKRN230022264-1A_HFCCGDSX7_concatenated_fwd.fq.gz #concatenate

#concatnenate des version du L3R1-rev
#cat Tzet_L3R1_EKRN230022264-1A_HF3NCDSX7_L3_2_rev-paired.fq.gz \
#    Tzet_L3R1_EKRN230022264-1A_HFCCGDSX7_L4_2_rev-paired.fq.gz > Tzet_L3R1_EKRN230022264-1A_HFCCGDSX7_concatenated_rev.fq.gz #concatenate

## Run Kraken2 for leg3R1
kraken2 --db standard_database \
	--paired Tzet_L3R1_EKRN230022264-1A_HFCCGDSX7_concatenated_fwd.fq.gz Tzet_L3R1_EKRN230022264-1A_HFCCGDSX7_concatenated_rev.fq.gz \
	--threads 24 \
	--classified-out leg3_R1_classified_reads#.txt \
	--unclassified-out leg3_R1_unclassified_reads#.txt \
	--report leg3_R1_kraken2_report.txt \
	--gzip-compressed

#concatnenate des version du L3R2-fwd
#cat Tzet_L3R2_EKRN230022267-1A_HF3NCDSX7_L3_1_fwd-paired.fq.gz \
#    Tzet_L3R2_EKRN230022267-1A_HFCCGDSX7_L4_1_fwd-paired.fq.gz > Tzet_L3R2_EKRN230022267-1A_HF3NCDSX7_concatenated_fwd.fq.gz #concatenate

#concatnenate des version du L3R2-rev
#cat Tzet_L3R2_EKRN230022267-1A_HF3NCDSX7_L3_2_rev-paired.fq.gz \
#    Tzet_L3R2_EKRN230022267-1A_HFCCGDSX7_L4_2_rev-paired.fq.gz > Tzet_L3R2_EKRN230022267-1A_HF3NCDSX7_concatenated_rev.fq.gz #concatenate

## Run Kraken2 for leg3R2
kraken2 --db standard_database \
	--paired Tzet_L3R2_EKRN230022267-1A_HF3NCDSX7_concatenated_fwd.fq.gz Tzet_L3R2_EKRN230022267-1A_HF3NCDSX7_concatenated_rev.fq.gz \
	--threads 24 \
	--classified-out leg3_R2_classified_reads#.txt \
	--unclassified-out leg3_R2_unclassified_reads#.txt \
	--report leg3_R2_kraken2_report.txt \
	--gzip-compressed

## Run Kraken2 for leg3R3
kraken2 --db standard_database \
	--paired Tzet_L3R3_EKRN230022270-1A_HF3NCDSX7_L3_1_fwd-paired.fq.gz Tzet_L3R3_EKRN230022270-1A_HF3NCDSX7_L3_2_rev-paired.fq.gz \
	--threads 24 \
	--classified-out leg3_R3_classified_reads#.txt \
	--unclassified-out leg3_R3_unclassified_reads#.txt \
	--report leg3_R3_kraken2_report.txt \
	--gzip-compressed
