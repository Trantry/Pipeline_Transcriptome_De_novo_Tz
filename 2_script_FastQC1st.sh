#!/bin/bash
#SBATCH --array=0-3 #(0-(n-1) où n: nombre de job en parallèle) ###
#SBATCH --job-name=Kraken2_Tetraripis ###
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=0-12:00:00 ###
#SBATCH --partition=E5
#SBATCH --output=/home/tbessonn/stdout/%A_%a.out # standard output file format
#SBATCH --error=/home/tbessonn/stderr/%A_%a.err # error file format

module use /applis/PSMN/debian11/E5/modules/all
 
cd /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/trimmomatic

PARAMS=(/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/3_kraken2/samples/*_concatenated_fwd.fq.gz) ####

f1="${PARAMS[$SLURM_ARRAY_TASK_ID]}"
f2=${PARAMS[$SLURM_ARRAY_TASK_ID]%%_concatenated_fwd.fq.gz}"_concatenated_rev.fq.gz"
name=$(basename "${PARAMS[$SLURM_ARRAY_TASK_ID]}" "_concatenated_fwd.fq.gz") 

fastqc \
  $f1 \
 	$f2 \
	-o fastqc2nd/

### Let's check data quality with fastqc

#L1R1
mkdir fastqc_2/
fastqc tzet_L1R1/v2/Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_L1_1_fwd-paired.fq.gz \
  tzet_L1R1/v2/Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_L1_2_rev-paired.fq.gz \
  tzet_L1R1/v2/Tzet_L1R1_EKRN230022262-1A_HFC7MDSX7_L3_1_fwd-paired.fq.gz \
  tzet_L1R1/v2/Tzet_L1R1_EKRN230022262-1A_HFC7MDSX7_L3_2_rev-paired.fq.gz \
  tzet_L1R1/v2/Tzet_L1R1_EKRN230022262-1A_HFCCLDSX7_L2_1_fwd-paired.fq.gz \
  tzet_L1R1/v2/Tzet_L1R1_EKRN230022262-1A_HFCCLDSX7_L2_2_rev-paired.fq.gz \
  -o fastqc_2/

#L1R2
fastqc tzet_L1R2/Tzet_L1R2_EKRN230022265-1A_HF3TCDSX7_L1_1_fwd-paired.fq.gz tzet_L1R2/Tzet_L1R2_EKRN230022265-1A_HF3TCDSX7_L1_2_rev-paired.fq.gz -o fastqc_2/

#L1R3
fastqc tzet_L1R3/Tzet_L1R3_EKRN230022268-1A_HF3TCDSX7_L1_1_fwd-paired.fq.gz tzet_L1R3/Tzet_L1R3_EKRN230022268-1A_HF3TCDSX7_L1_2_rev-paired.fq.gz -o fastqc_2/

#L2R1
fastqc tzet_L2R1/v2/Tzet_L2R1_EKRN230022263-1A_HF3NCDSX7_L3_1_fwd-paired.fq.gz tzet_L2R1/v2/Tzet_L2R1_EKRN230022263-1A_HF3NCDSX7_L3_2_rev-paired.fq.gz tzet_L2R1/v2/Tzet_L2R1_EKRN230022263-1A_HFCCGDSX7_L4_1_fwd-paired.fq.gz tzet_L2R1/v2/Tzet_L2R1_EKRN230022263-1A_HFCCGDSX7_L4_2_rev-paired.fq.gz -o fastqc_2/

#L2R2
fastqc tzet_L2R2/Tzet_L2R2_EKRN230022266-1A_HF3NCDSX7_L3_1_fwd-paired.fq.gz tzet_L2R2/Tzet_L2R2_EKRN230022266-1A_HF3NCDSX7_L3_2_rev-paired.fq.gz -o fastqc_2/

#L2R3
fastqc tzet_L2R3/Tzet_L2R3_EKRN230022269-1A_HF3NCDSX7_L3_1_fwd-paired.fq.gz tzet_L2R3/Tzet_L2R3_EKRN230022269-1A_HF3NCDSX7_L3_2_rev-paired.fq.gz -o fastqc_2/

#L3R1
fastqc tzet_L3R1/v2/Tzet_L3R1_EKRN230022264-1A_HF3NCDSX7_L3_1_fwd-paired.fq.gz tzet_L3R1/v2/Tzet_L3R1_EKRN230022264-1A_HF3NCDSX7_L3_2_rev-paired.fq.gz tzet_L3R1/v2/Tzet_L3R1_EKRN230022264-1A_HFCCGDSX7_L4_1_fwd-paired.fq.gz tzet_L3R1/v2/Tzet_L3R1_EKRN230022264-1A_HFCCGDSX7_L4_2_rev-paired.fq.gz -o fastqc_2/

#L3R2
fastqc tzet_L3R2/v2/Tzet_L3R2_EKRN230022267-1A_HF3NCDSX7_L3_1_fwd-paired.fq.gz tzet_L3R2/v2/Tzet_L3R2_EKRN230022267-1A_HF3NCDSX7_L3_2_rev-paired.fq.gz tzet_L3R2/v2/Tzet_L3R2_EKRN230022267-1A_HFCCGDSX7_L4_1_fwd-paired.fq.gz tzet_L3R2/v2/Tzet_L3R2_EKRN230022267-1A_HFCCGDSX7_L4_2_rev-paired.fq.gz -o fastqc_2/

#L3R3
fastqc tzet_L3R3/Tzet_L3R3_EKRN230022270-1A_HF3NCDSX7_L3_1_fwd-paired.fq.gz tzet_L3R3/Tzet_L3R3_EKRN230022270-1A_HF3NCDSX7_L3_2_rev-paired.fq.gz -o fastqc_2/
