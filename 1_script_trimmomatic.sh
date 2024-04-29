#!/bin/bash
#SBATCH --array=0-8#(0-(n-1) où n: nombre de job en parallèle) ###
#SBATCH --job-name=Trimmomatic_Tetraripis ###
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=0-12:00:00 ###
#SBATCH --partition=E5
#SBATCH --output=/home/tbessonn/stdout/%A_%a.out # standard output file format
#SBATCH --error=/home/tbessonn/stderr/%A_%a.err # error file format

module use /applis/PSMN/debian11/E5/modules/all

cd /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/trimmomatic


PARAMS=(/Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F003/X204SC23042056-Z01-F003/01.RawData/Tzet_L1R1/*_R1_001.fastq.gz)

f1="${PARAMS[$SLURM_ARRAY_TASK_ID]}"
f2=${PARAMS[$SLURM_ARRAY_TASK_ID]%%_R1_001.fastq.gz}"_R2_001.fastq.gz"
namef1=$(basename "${PARAMS[$SLURM_ARRAY_TASK_ID]}" "R1_001.fastq.gz" | sed -E 's/_S[0-9]+_L[0-9]+_//')

TrimmomaticPE \
  $f1 \
  $f2 \
  2_TrimmomaticPE/$namef1"_1_paired.fq.gz" \
  2_TrimmomaticPE/$namef1"_1_unpaired.fq.gz" 2_TrimmomaticPE/$namef1"_2_paired.fq.gz" \
  2_TrimmomaticPE/$namef1"_2_unpaired.fq.gz" \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

  
### This is a dataset containing rna-seq data for the 3 pairs of legs of Tetraripis zetteli, in 3 replicates each.
### Let's trim the rna-seq data with trimmomatic

# Trim adapters for L1R1            <input1>                                                                                                                                               <input 2>                                                                                                                                              <paired output 1>                                                      <unpaired output 1>                                                   <paired output 2>                                                    <unpaired output 2>                                                     
TrimmomaticPE -threads 4 -phred33 /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F002/X204SC23042056-Z01-F002_02/01.RawData/Tzet_L1R1/Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_L1_1.fq.gz /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F002/X204SC23042056-Z01-F002_02/01.RawData/Tzet_L1R1/Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_L1_2.fq.gz tzet_L1R1/Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_L1_1_fwd-paired.fq.gz tzet_L1R1/Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_L1_1_fwd-unpaired.fq.gz tzet_L1R1/Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_L1_2_rev-paired.fq.gz tzet_L1R1/Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_L1_2_rev-unpaired.fq.gz ILLUMINACLIP:/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/trimmomatic/TruSeq3-PE-2.fa:2:30:10 LEADING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>&1 | tee tzet_L1R1/trimmomatic_tzet_L1R1_test.out

# Trim adapters for L1R2
TrimmomaticPE -threads 4 -phred33 /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F002/X204SC23042056-Z01-F002_02/01.RawData/Tzet_L1R2/Tzet_L1R2_EKRN230022265-1A_HF3TCDSX7_L1_1.fq.gz /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F002/X204SC23042056-Z01-F002_02/01.RawData/Tzet_L1R2/Tzet_L1R2_EKRN230022265-1A_HF3TCDSX7_L1_2.fq.gz tzet_L1R2/Tzet_L1R2_EKRN230022265-1A_HF3TCDSX7_L1_1_fwd-paired.fq.gz tzet_L1R2/Tzet_L1R2_EKRN230022265-1A_HF3TCDSX7_L1_1_fwd-unpaired.fq.gz tzet_L1R2/Tzet_L1R2_EKRN230022265-1A_HF3TCDSX7_L1_2_rev-paired.fq.gz tzet_L1R2/Tzet_L1R2_EKRN230022265-1A_HF3TCDSX7_L1_2_rev-unpaired.fq.gz ILLUMINACLIP:/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/trimmomatic/TruSeq3-PE-2.fa:2:30:10 LEADING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>&1 | tee tzet_L1R2/trimmomatic_tzet_L1R2_test.out

# Trim adapters for L1R3
TrimmomaticPE -threads 4 -phred33 /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F002/X204SC23042056-Z01-F002_02/01.RawData/Tzet_L1R3/Tzet_L1R3_EKRN230022268-1A_HF3TCDSX7_L1_1.fq.gz /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F002/X204SC23042056-Z01-F002_02/01.RawData/Tzet_L1R3/Tzet_L1R3_EKRN230022268-1A_HF3TCDSX7_L1_2.fq.gz tzet_L1R3/Tzet_L1R3_EKRN230022268-1A_HF3TCDSX7_L1_1_fwd-paired.fq.gz tzet_L1R3/Tzet_L1R3_EKRN230022268-1A_HF3TCDSX7_L1_1_fwd-unpaired.fq.gz tzet_L1R3/Tzet_L1R3_EKRN230022268-1A_HF3TCDSX7_L1_rev-paired.fq.gz tzet_L1R3/Tzet_L1R3_EKRN230022268-1A_HF3TCDSX7_L1_2_rev-unpaired.fq.gz ILLUMINACLIP:/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/trimmomatic/TruSeq3-PE-2.fa:2:30:10 LEADING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>&1 | tee tzet_L1R3/trimmomatic_tzet_L1R3_test.out

# Trim adapters for L2R1
TrimmomaticPE -threads 4 -phred33 /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F002/X204SC23042056-Z01-F002_01/01.RawData/Tzet_L2R1/Tzet_L2R1_EKRN230022263-1A_HF3NCDSX7_L3_1.fq.gz /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F002/X204SC23042056-Z01-F002_01/01.RawData/Tzet_L2R1/Tzet_L2R1_EKRN230022263-1A_HF3NCDSX7_L3_2.fq.gz tzet_L2R1/Tzet_L2R1_EKRN230022263-1A_HF3NCDSX7_L3_1_fwd-paired.fq.gz tzet_L2R1/Tzet_L2R1_EKRN230022263-1A_HF3NCDSX7_L3_1_fwd-unpaired.fq.gz tzet_L2R1/Tzet_L2R1_EKRN230022263-1A_HF3NCDSX7_L3_2_rev-paired.fq.gz tzet_L2R1/Tzet_L2R1_EKRN230022263-1A_HF3NCDSX7_L3_2_rev-unpaired.fq.gz ILLUMINACLIP:/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/trimmomatic/TruSeq3-PE-2.fa:2:30:10 LEADING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>&1 | tee tzet_L2R1/trimmomatic_tzet_L2R1_test.out

# Trim adapters for L2R2
TrimmomaticPE -threads 4 -phred33 /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F002/X204SC23042056-Z01-F002_01/01.RawData/Tzet_L2R2/Tzet_L2R2_EKRN230022266-1A_HF3NCDSX7_L3_1.fq.gz /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F002/X204SC23042056-Z01-F002_01/01.RawData/Tzet_L2R2/Tzet_L2R2_EKRN230022266-1A_HF3NCDSX7_L3_2.fq.gz tzet_L2R2/Tzet_L2R2_EKRN230022266-1A_HF3NCDSX7_L3_1_fwd-paired.fq.gz tzet_L2R2/Tzet_L2R2_EKRN230022266-1A_HF3NCDSX7_L3_1_fwd-unpaired.fq.gz tzet_L2R2/Tzet_L2R2_EKRN230022266-1A_HF3NCDSX7_L3_2_rev-paired.fq.gz tzet_L2R2/Tzet_L2R2_EKRN230022266-1A_HF3NCDSX7_L3_2_rev-unpaired.fq.gz ILLUMINACLIP:/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/trimmomatic/TruSeq3-PE-2.fa:2:30:10 LEADING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>&1 | tee tzet_L2R2/trimmomatic_tzet_L2R2_test.out

# Trim adapters for L2R3
TrimmomaticPE -threads 4 -phred33 /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F002/X204SC23042056-Z01-F002_01/01.RawData/Tzet_L2R3/Tzet_L2R3_EKRN230022269-1A_HF3NCDSX7_L3_1.fq.gz /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F002/X204SC23042056-Z01-F002_01/01.RawData/Tzet_L2R3/Tzet_L2R3_EKRN230022269-1A_HF3NCDSX7_L3_2.fq.gz tzet_L2R3/Tzet_L2R3_EKRN230022269-1A_HF3NCDSX7_L3_1_fwd-paired.fq.gz tzet_L2R3/Tzet_L2R3_EKRN230022269-1A_HF3NCDSX7_L3_1_fwd-unpaired.fq.gz tzet_L2R3/Tzet_L2R3_EKRN230022269-1A_HF3NCDSX7_L3_2_rev-paired.fq.gz tzet_L2R3/Tzet_L2R3_EKRN230022269-1A_HF3NCDSX7_L3_2_rev-unpaired.fq.gz ILLUMINACLIP:/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/trimmomatic/TruSeq3-PE-2.fa:2:30:10 LEADING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>&1 | tee tzet_L2R3/trimmomatic_tzet_L2R3_test.out

# Trim adapters for L3R1
TrimmomaticPE -threads 4 -phred33 /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F002/X204SC23042056-Z01-F002_01/01.RawData/Tzet_L3R1/Tzet_L3R1_EKRN230022264-1A_HF3NCDSX7_L3_1.fq.gz /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F002/X204SC23042056-Z01-F002_01/01.RawData/Tzet_L3R1/Tzet_L3R1_EKRN230022264-1A_HF3NCDSX7_L3_2.fq.gz tzet_L3R1/Tzet_L3R1_EKRN230022264-1A_HF3NCDSX7_L3_1_fwd-paired.fq.gz tzet_L3R1/Tzet_L3R1_EKRN230022264-1A_HF3NCDSX7_L3_1_fwd-unpaired.fq.gz tzet_L3R1/Tzet_L3R1_EKRN230022264-1A_HF3NCDSX7_L3_2_rev-paired.fq.gz tzet_L3R1/Tzet_L3R1_EKRN230022264-1A_HF3NCDSX7_L3_2_rev-unpaired.fq.gz ILLUMINACLIP:/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/trimmomatic/TruSeq3-PE-2.fa:2:30:10 LEADING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>&1 | tee tzet_L3R1/trimmomatic_tzet_L3R1_test.out

# Trim adapters for L3R2
TrimmomaticPE -threads 4 -phred33 /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F002/X204SC23042056-Z01-F002_01/01.RawData/Tzet_L3R2/Tzet_L3R2_EKRN230022267-1A_HF3NCDSX7_L3_1.fq.gz /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F002/X204SC23042056-Z01-F002_01/01.RawData/Tzet_L3R2/Tzet_L3R2_EKRN230022267-1A_HF3NCDSX7_L3_2.fq.gz tzet_L3R2/Tzet_L3R2_EKRN230022267-1A_HF3NCDSX7_L3_1_fwd-paired.fq.gz tzet_L3R2/Tzet_L3R2_EKRN230022267-1A_HF3NCDSX7_L3_1_fwd-unpaired.fq.gz tzet_L3R2/Tzet_L3R2_EKRN230022267-1A_HF3NCDSX7_L3_2_rev-paired.fq.gz tzet_L3R2/Tzet_L3R2_EKRN230022267-1A_HF3NCDSX7_L3_2_rev-unpaired.fq.gz ILLUMINACLIP:/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/trimmomatic/TruSeq3-PE-2.fa:2:30:10 LEADING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>&1 | tee tzet_L3R2/trimmomatic_tzet_L3R2_test.out

# Trim adapters for L3R3
TrimmomaticPE -threads 4 -phred33 /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F002/X204SC23042056-Z01-F002_01/01.RawData/Tzet_L3R3/Tzet_L3R3_EKRN230022270-1A_HF3NCDSX7_L3_1.fq.gz /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F002/X204SC23042056-Z01-F002_01/01.RawData/Tzet_L3R3/Tzet_L3R3_EKRN230022270-1A_HF3NCDSX7_L3_2.fq.gz tzet_L3R3/Tzet_L3R3_EKRN230022270-1A_HF3NCDSX7_L3_1_fwd-paired.fq.gz tzet_L3R3/Tzet_L3R3_EKRN230022270-1A_HF3NCDSX7_L3_1_fwd-unpaired.fq.gz tzet_L3R3/Tzet_L3R3_EKRN230022270-1A_HF3NCDSX7_L3_2_rev-paired.fq.gz tzet_L3R3/Tzet_L3R3_EKRN230022270-1A_HF3NCDSX7_L3_2_rev-unpaired.fq.gz ILLUMINACLIP:/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/trimmomatic/TruSeq3-PE-2.fa:2:30:10 LEADING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>&1 | tee tzet_L3R3/trimmomatic_tzet_L3R3_test.out


### Some of the rna-seq data above were not optimal so novogene proceeded to a new sequencing round, which I called v2 in the data folder.

# Trim adapters for L1R1. We have 3 pairs of fq files so we need to do this 3 times
TrimmomaticPE -threads 4 -phred33 /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F003/X204SC23042056-Z01-F003/01.RawData/Tzet_L1R1/Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_L1_1.fq.gz /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F003/X204SC23042056-Z01-F003/01.RawData/Tzet_L1R1/Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_L1_2.fq.gz  tzet_L1R1/v2/Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_L1_1_fwd-paired.fq.gz tzet_L1R1/v2/Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_L1_1_fwd-unpaired.fq.gz tzet_L1R1/v2/Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_L1_2_rev-paired.fq.gz tzet_L1R1/v2/Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_L1_2_rev-unpaired.fq.gz ILLUMINACLIP:/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/trimmomatic/TruSeq3-PE-2.fa:2:30:10 LEADING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>&1 | tee tzet_L1R1/v2/trimmomatic_tzet_L1R1_v2_test.out

TrimmomaticPE -threads 4 -phred33 /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F003/X204SC23042056-Z01-F003/01.RawData/Tzet_L1R1/Tzet_L1R1_EKRN230022262-1A_HFC7MDSX7_L3_1.fq.gz /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F003/X204SC23042056-Z01-F003/01.RawData/Tzet_L1R1/Tzet_L1R1_EKRN230022262-1A_HFC7MDSX7_L3_2.fq.gz  tzet_L1R1/v2/Tzet_L1R1_EKRN230022262-1A_HFC7MDSX7_L3_1_fwd-paired.fq.gz tzet_L1R1/v2/Tzet_L1R1_EKRN230022262-1A_HFC7MDSX7_L3_1_fwd-unpaired.fq.gz tzet_L1R1/v2/Tzet_L1R1_EKRN230022262-1A_HFC7MDSX7_L3_2_rev-paired.fq.gz tzet_L1R1/v2/Tzet_L1R1_EKRN230022262-1A_HFC7MDSX7_L3_2_rev-unpaired.fq.gz ILLUMINACLIP:/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/trimmomatic/TruSeq3-PE-2.fa:2:30:10 LEADING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>&1 | tee tzet_L1R1/v2/trimmomatic_tzet_L1R1_v3_test.out

TrimmomaticPE -threads 4 -phred33 /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F003/X204SC23042056-Z01-F003/01.RawData/Tzet_L1R1/Tzet_L1R1_EKRN230022262-1A_HFCCLDSX7_L2_1.fq.gz /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F003/X204SC23042056-Z01-F003/01.RawData/Tzet_L1R1/Tzet_L1R1_EKRN230022262-1A_HFCCLDSX7_L2_2.fq.gz tzet_L1R1/v2/Tzet_L1R1_EKRN230022262-1A_HFCCLDSX7_L2_1_fwd-paired.fq.gz tzet_L1R1/v2/Tzet_L1R1_EKRN230022262-1A_HFCCLDSX7_L2_1_fwd-unpaired.fq.gz tzet_L1R1/v2/Tzet_L1R1_EKRN230022262-1A_HFCCLDSX7_L2_2_rev-paired.fq.gz tzet_L1R1/v2/Tzet_L1R1_EKRN230022262-1A_HFCCLDSX7_L2_2_rev-unpaired.fq.gz ILLUMINACLIP:/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/trimmomatic/TruSeq3-PE-2.fa:2:30:10 LEADING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>&1 | tee tzet_L1R1/v2/trimmomatic_tzet_L1R1_v4_test.out


# Trim adapters for L2R1. We have 2 pairs of fq files so we need to do this 2 times
TrimmomaticPE -threads 4 -phred33 /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F004/X204SC23042056-Z01-F004/01.RawData/Tzet_L2R1/Tzet_L2R1_EKRN230022263-1A_HF3NCDSX7_L3_1.fq.gz /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F004/X204SC23042056-Z01-F004/01.RawData/Tzet_L2R1/Tzet_L2R1_EKRN230022263-1A_HF3NCDSX7_L3_2.fq.gz tzet_L2R1/v2/Tzet_L2R1_EKRN230022263-1A_HF3NCDSX7_L3_1_fwd-paired.fq.gz tzet_L2R1/v2/Tzet_L2R1_EKRN230022263-1A_HF3NCDSX7_L3_1_fwd-unpaired.fq.gz tzet_L2R1/v2/Tzet_L2R1_EKRN230022263-1A_HF3NCDSX7_L3_2_rev-paired.fq.gz tzet_L2R1/v2/Tzet_L2R1_EKRN230022263-1A_HF3NCDSX7_L3_2_rev-unpaired.fq.gz ILLUMINACLIP:/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/trimmomatic/TruSeq3-PE-2.fa:2:30:10 LEADING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>&1 | tee tzet_L2R1/v2/trimmomatic_tzet_L2R1_v2_test.out

TrimmomaticPE -threads 4 -phred33 /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F004/X204SC23042056-Z01-F004/01.RawData/Tzet_L2R1/Tzet_L2R1_EKRN230022263-1A_HFCCGDSX7_L4_1.fq.gz /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F004/X204SC23042056-Z01-F004/01.RawData/Tzet_L2R1/Tzet_L2R1_EKRN230022263-1A_HFCCGDSX7_L4_2.fq.gz tzet_L2R1/v2/Tzet_L2R1_EKRN230022263-1A_HFCCGDSX7_L4_1_fwd-paired.fq.gz tzet_L2R1/v2/Tzet_L2R1_EKRN230022263-1A_HFCCGDSX7_L4_1_fwd-unpaired.fq.gz tzet_L2R1/v2/Tzet_L2R1_EKRN230022263-1A_HFCCGDSX7_L4_2_rev-paired.fq.gz tzet_L2R1/v2/Tzet_L2R1_EKRN230022263-1A_HFCCGDSX7_L4_2_rev-unpaired.fq.gz ILLUMINACLIP:/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/trimmomatic/TruSeq3-PE-2.fa:2:30:10 LEADING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>&1 | tee tzet_L2R1/v2/trimmomatic_tzet_L2R1_v3_test.out



# Trim adapters for L3R1. We have 2 pairs of fq files so we need to do this 2 times
TrimmomaticPE -threads 4 -phred33 /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F004/X204SC23042056-Z01-F004/01.RawData/Tzet_L3R1/Tzet_L3R1_EKRN230022264-1A_HF3NCDSX7_L3_1.fq.gz /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F004/X204SC23042056-Z01-F004/01.RawData/Tzet_L3R1/Tzet_L3R1_EKRN230022264-1A_HF3NCDSX7_L3_2.fq.gz tzet_L3R1/v2/Tzet_L3R1_EKRN230022264-1A_HF3NCDSX7_L3_1_fwd-paired.fq.gz tzet_L3R1/v2/Tzet_L3R1_EKRN230022264-1A_HF3NCDSX7_L3_1_fwd-unpaired.fq.gz tzet_L3R1/v2/Tzet_L3R1_EKRN230022264-1A_HF3NCDSX7_L3_2_rev-paired.fq.gz tzet_L3R1/v2/Tzet_L3R1_EKRN230022264-1A_HF3NCDSX7_L3_2_rev-unpaired.fq.gz ILLUMINACLIP:/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/trimmomatic/TruSeq3-PE-2.fa:2:30:10 LEADING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>&1 | tee tzet_L3R1/v2/trimmomatic_tzet_L3R1_v2_test.out

TrimmomaticPE -threads 4 -phred33 /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F004/X204SC23042056-Z01-F004/01.RawData/Tzet_L3R1/Tzet_L3R1_EKRN230022264-1A_HFCCGDSX7_L4_1.fq.gz /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F004/X204SC23042056-Z01-F004/01.RawData/Tzet_L3R1/Tzet_L3R1_EKRN230022264-1A_HFCCGDSX7_L4_2.fq.gz tzet_L3R1/v2/Tzet_L3R1_EKRN230022264-1A_HFCCGDSX7_L4_1_fwd-paired.fq.gz tzet_L3R1/v2/Tzet_L3R1_EKRN230022264-1A_HFCCGDSX7_L4_1_fwd-unpaired.fq.gz tzet_L3R1/v2/Tzet_L3R1_EKRN230022264-1A_HFCCGDSX7_L4_2_rev-paired.fq.gz tzet_L3R1/v2/Tzet_L3R1_EKRN230022264-1A_HFCCGDSX7_L4_2_rev-unpaired.fq.gz ILLUMINACLIP:/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/trimmomatic/TruSeq3-PE-2.fa:2:30:10 LEADING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>&1 | tee tzet_L3R1/v2/trimmomatic_tzet_L3R1_v3_test.out


# Trim adapters for L3R2. We have 2 pairs of fq files so we need to do this 2 times
TrimmomaticPE -threads 4 -phred33 /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F004/X204SC23042056-Z01-F004/01.RawData/Tzet_L3R2/Tzet_L3R2_EKRN230022267-1A_HF3NCDSX7_L3_1.fq.gz /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F004/X204SC23042056-Z01-F004/01.RawData/Tzet_L3R2/Tzet_L3R2_EKRN230022267-1A_HF3NCDSX7_L3_2.fq.gz tzet_L3R2/v2/Tzet_L3R2_EKRN230022267-1A_HF3NCDSX7_L3_1_fwd-paired.fq.gz tzet_L3R2/v2/Tzet_L3R2_EKRN230022267-1A_HF3NCDSX7_L3_1_fwd-unpaired.fq.gz tzet_L3R2/v2/Tzet_L3R2_EKRN230022267-1A_HF3NCDSX7_L3_2_rev-paired.fq.gz tzet_L3R2/v2/Tzet_L3R2_EKRN230022267-1A_HF3NCDSX7_L3_2_rev-unpaired.fq.gz ILLUMINACLIP:/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/trimmomatic/TruSeq3-PE-2.fa:2:30:10 LEADING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>&1 | tee tzet_L3R2/v2/trimmomatic_tzet_L3R2_v2_test.out

TrimmomaticPE -threads 4 -phred33 /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F004/X204SC23042056-Z01-F004/01.RawData/Tzet_L3R2/Tzet_L3R2_EKRN230022267-1A_HFCCGDSX7_L4_1.fq.gz /Xnfs/khila/omics_data/Novogene/X204SC23042056-Z01-F004/X204SC23042056-Z01-F004/01.RawData/Tzet_L3R2/Tzet_L3R2_EKRN230022267-1A_HFCCGDSX7_L4_2.fq.gz tzet_L3R2/v2/Tzet_L3R2_EKRN230022267-1A_HFCCGDSX7_L4_1_fwd-paired.fq.gz tzet_L3R2/v2/Tzet_L3R2_EKRN230022267-1A_HFCCGDSX7_L4_1_fwd-unpaired.fq.gz tzet_L3R2/v2/Tzet_L3R2_EKRN230022267-1A_HFCCGDSX7_L4_2_rev-paired.fq.gz tzet_L3R2/v2/Tzet_L3R2_EKRN230022267-1A_HFCCGDSX7_L4_2_rev-unpaired.fq.gz ILLUMINACLIP:/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/trimmomatic/TruSeq3-PE-2.fa:2:30:10 LEADING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>&1 | tee tzet_L3R2/v2/trimmomatic_tzet_L3R2_v3_test.out



### Let's check data quality with fastqc

#L1R1
mkdir fastqc_2/
fastqc tzet_L1R1/v2/Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_L1_1_fwd-paired.fq.gz tzet_L1R1/v2/Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7_L1_2_rev-paired.fq.gz tzet_L1R1/v2/Tzet_L1R1_EKRN230022262-1A_HFC7MDSX7_L3_1_fwd-paired.fq.gz tzet_L1R1/v2/Tzet_L1R1_EKRN230022262-1A_HFC7MDSX7_L3_2_rev-paired.fq.gz tzet_L1R1/v2/Tzet_L1R1_EKRN230022262-1A_HFCCLDSX7_L2_1_fwd-paired.fq.gz tzet_L1R1/v2/Tzet_L1R1_EKRN230022262-1A_HFCCLDSX7_L2_2_rev-paired.fq.gz -o fastqc_2/

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

