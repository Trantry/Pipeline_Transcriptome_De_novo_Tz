# -*- coding: utf-8 -*-

def filter_sequences(file_path, output_file_path):
    with open(file_path, 'r') as file:
        transcriptome = file.read()
    filtered_transcriptome = ""
    sequences = transcriptome.split(">")
    for sequence in sequences:
        if sequence:
            header, seq = sequence.split("\n", 1)
            length = int(header.split("len=")[-1].split()[0])
            if length >= 300:
                filtered_transcriptome += ">" + sequence + "\n"
    with open(output_file_path, 'w') as output_file:
        output_file.write(filtered_transcriptome)

input_file_path = "/Users/trystan.bessonnier/Desktop/DESeq2_Transcriptome/trinity_tzet_complete_transcriptome.Trinity.fasta"
output_file_path = "/Users/trystan.bessonnier/Desktop/DESeq2_Transcriptome/filtered_transcriptome.fasta"

filter_sequences(input_file_path, output_file_path)
