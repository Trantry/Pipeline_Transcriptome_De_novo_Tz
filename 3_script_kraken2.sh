#!/bin/bash
#SBATCH --array=0-3 #(0-(n-1) où n: nombre de job en parallèle) ###
#SBATCH --job-name=Kraken2_Tetraripis ###
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=0-12:00:00 ###
#SBATCH --partition=E5
#SBATCH --output=/home/tbessonn/stdout/%A_%a.out # standard output file format
#SBATCH --error=/home/tbessonn/stderr/%A_%a.err # error file format

module use /applis/PSMN/debian11/E5/modules/all
 
cd /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/3_kraken2/

## Run Kraken2 for leg_concatenated
PARAMS=(/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/3_kraken2/samples/*_concatenated_fwd.fq.gz)
f1="${PARAMS[$SLURM_ARRAY_TASK_ID]}"
f2=${PARAMS[$SLURM_ARRAY_TASK_ID]%%_concatenated_fwd.fq.gz}"_concatenated_rev.fq.gz"
name=$(basename "${PARAMS[$SLURM_ARRAY_TASK_ID]}" "_concatenated_fwd.fq.gz") 

## Run Kraken2 for leg1R1
kraken2 --db standard_database \
	--paired $f1 $f2 \
	--threads 8 \
	--classified-out samples/classified/$name#"classified_reads.fq" \
	--unclassified-out samples/unclassified/$name#"classified_reads.fq" \
	--report samples/report/$name"_report.txt" \
	--gzip-compressed 

## Run Kraken2 for leg_non_concatenated
#PARAMS=(/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/3_kraken2/samples/*_1_fwd-paired.fq.gz)
#f1="${PARAMS[$SLURM_ARRAY_TASK_ID]}"
#f2=${PARAMS[$SLURM_ARRAY_TASK_ID]%%_1_fwd-paired.fq.gz}"_2_rev-paired.fq.gz"
#name=$(basename "${PARAMS[$SLURM_ARRAY_TASK_ID]}" "_1_fwd-paired.fq.gz") 

## Run Kraken2 for leg1R2
#kraken2 --db standard_database \
#	--paired $f1 $f2 \
#	--threads 8 \
#	--classified-out samples/classified/$name#"classified_reads.fq" \
#	--unclassified-out samples/unclassified/$name#"classified_reads.fq" \
#	--report samples/report/$name"_report.txt" \
#	--gzip-compressed 








#concatnenate des version du L1R1-fwd
#cat Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_L1_1_fwd-paired.fq.gz \
#    Tzet_L1R1_EKRN230022262-1A_HFC7MDSX7_L3_1_fwd-paired.fq.gz \
#    Tzet_L1R1_EKRN230022262-1A_HFCCLDSX7_L2_1_fwd-paired.fq.gz > Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_concatenated_fwd.fq.gz #concatenate

#concatnenate des version du L1R1-rev
#cat Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_L1_2_rev-paired.fq.gz \
#    Tzet_L1R1_EKRN230022262-1A_HFC7MDSX7_L3_2_rev-paired.fq.gz \
#    Tzet_L1R1_EKRN230022262-1A_HFCCLDSX7_L2_2_rev-paired.fq.gz > Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_concatenated_rev.fq.gz #concatenate

## Run Kraken2 for leg1R1

## Run Kraken2 for leg1R2

## Run Kraken2 for leg1R3

#concatnenate des version du L2R1-fwd
#cat Tzet_L2R1_EKRN230022263-1A_HF3NCDSX7_L3_1_fwd-paired.fq.gz \
#    Tzet_L2R1_EKRN230022263-1A_HFCCGDSX7_L4_1_fwd-paired.fq.gz > Tzet_L2R1_EKRN230022263-1A_HFCCGDSX7_concatenated_fwd.fq.gz #concatenate

#concatnenate des version du L2R1-fwd
#cat Tzet_L2R1_EKRN230022263-1A_HF3NCDSX7_L3_2_rev-paired.fq.gz \
#    Tzet_L2R1_EKRN230022263-1A_HFCCGDSX7_L4_2_rev-paired.fq.gz > Tzet_L2R1_EKRN230022263-1A_HFCCGDSX7_concatenated_rev.fq.gz #concatenate

## Run Kraken2 for leg2R1

## Run Kraken2 for leg2R2

## Run Kraken2 for leg2R3

#concatnenate des version du L3R1-fwd
#cat Tzet_L3R1_EKRN230022264-1A_HF3NCDSX7_L3_1_fwd-paired.fq.gz \
#    Tzet_L3R1_EKRN230022264-1A_HFCCGDSX7_L4_1_fwd-paired.fq.gz > Tzet_L3R1_EKRN230022264-1A_HFCCGDSX7_concatenated_fwd.fq.gz #concatenate

#concatnenate des version du L3R1-rev
#cat Tzet_L3R1_EKRN230022264-1A_HF3NCDSX7_L3_2_rev-paired.fq.gz \
#    Tzet_L3R1_EKRN230022264-1A_HFCCGDSX7_L4_2_rev-paired.fq.gz > Tzet_L3R1_EKRN230022264-1A_HFCCGDSX7_concatenated_rev.fq.gz #concatenate

## Run Kraken2 for leg3R1

#concatnenate des version du L3R2-fwd
#cat Tzet_L3R2_EKRN230022267-1A_HF3NCDSX7_L3_1_fwd-paired.fq.gz \
#    Tzet_L3R2_EKRN230022267-1A_HFCCGDSX7_L4_1_fwd-paired.fq.gz > Tzet_L3R2_EKRN230022267-1A_HF3NCDSX7_concatenated_fwd.fq.gz #concatenate

#concatnenate des version du L3R2-rev
#cat Tzet_L3R2_EKRN230022267-1A_HF3NCDSX7_L3_2_rev-paired.fq.gz \
#    Tzet_L3R2_EKRN230022267-1A_HFCCGDSX7_L4_2_rev-paired.fq.gz > Tzet_L3R2_EKRN230022267-1A_HF3NCDSX7_concatenated_rev.fq.gz #concatenate

## Run Kraken2 for leg3R2

## Run Kraken2 for leg3R3



## Copy the paired reads
#cp ../trinity/all_legs_trinity/*.fq.gz .
