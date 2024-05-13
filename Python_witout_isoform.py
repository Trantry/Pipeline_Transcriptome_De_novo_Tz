#conda activate BUSCO

#conda install pandas
#python /home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/1_kallisto/kallisto_output/filtered_isoform.py

import pandas as pd
import os

def filter_longest_isoforms(input_file, output_file):
    # Read the TSV file into a pandas DataFrame
    df = pd.read_csv(input_file, sep='\t')

    # Split the 'target_id' into gene and isoform IDs
    df[['gene_id', 'isoform_id']] = df['target_id'].str.rsplit('_', n=1, expand=True)

    # Group by 'gene_id' and keep the row with the max 'length'
    df_max_length = df.loc[df.groupby('gene_id')['length'].idxmax()]

    # Drop the 'gene_id' and 'isoform_id' columns
    df_max_length = df_max_length.drop(columns=['gene_id', 'isoform_id'])

    # Write the DataFrame to a new TSV file
    df_max_length.to_csv(output_file, sep='\t', index=False)

def process_files(input_files):
    for input_file in input_files:
        # Generate output file name based on input file name
        base_name = os.path.basename(input_file)
        name_without_extension = os.path.splitext(base_name)[0]
        output_file = os.path.join(os.path.dirname(input_file), f"{name_without_extension}_without_isoforms.tsv")
        
        filter_longest_isoforms(input_file, output_file)

# Usage
input_files = [
    "/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/1_kallisto/kallisto_output/Tzet_L1R1_EKRN230022262-1A_HF3TCDSX7__sortmerna_non_rRNA_/abundance.tsv",
    "/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/1_kallisto/kallisto_output/Tzet_L1R2_EKRN230022265-1A_HF3TCDSX7_L1__sortmerna_non_rRNA_/abundance.tsv",
    "/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/1_kallisto/kallisto_output/Tzet_L1R3_EKRN230022268-1A_HF3TCDSX7_L1__sortmerna_non_rRNA_/abundance.tsv",
    "/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/1_kallisto/kallisto_output/Tzet_L2R1_EKRN230022263-1A_HFCCGDSX7__sortmerna_non_rRNA_/abundance.tsv",
    "/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/1_kallisto/kallisto_output/Tzet_L2R2_EKRN230022266-1A_HF3NCDSX7_L3__sortmerna_non_rRNA_/abundance.tsv",
    "/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/1_kallisto/kallisto_output/Tzet_L2R3_EKRN230022269-1A_HF3NCDSX7_L3__sortmerna_non_rRNA_/abundance.tsv",
    "/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/1_kallisto/kallisto_output/Tzet_L3R1_EKRN230022264-1A_HFCCGDSX7__sortmerna_non_rRNA_/abundance.tsv",
    "/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/1_kallisto/kallisto_output/Tzet_L3R2_EKRN230022267-1A_HF3NCDSX7__sortmerna_non_rRNA_/abundance.tsv",
    "/home/tbessonn/Tetraripis_zetteli/tetraripis_zetteli_leg_transcriptome_novogene/8_post_assembly_QC/1.5_complete-transcriptome/1_kallisto/kallisto_output/Tzet_L3R3_EKRN230022270-1A_HF3NCDSX7_L3__sortmerna_non_rRNA_/abundance.tsv"
]

process_files(input_files)
