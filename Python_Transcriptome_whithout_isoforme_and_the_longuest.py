#source path/to/venv/bin/activate

from Bio import SeqIO
import re

def filter_longest_isoforms(input_path, output_path):
    record_dict = SeqIO.to_dict(SeqIO.parse(input_path, "fasta"))
    longest_isoforms = {}

    for record_id, record in record_dict.items():
        gene_id = "_".join(record_id.split("_")[:-1])
        if gene_id not in longest_isoforms or len(record.seq) > len(longest_isoforms[gene_id].seq):
            longest_isoforms[gene_id] = record

    # Modifier les identifiants de séquence
    for record_id, record in longest_isoforms.items():
        record.id = re.sub(r"_i\d+$", "", record.id)
        record.description = record.id

        # Conserver les informations "len=numéro" et "path="
        len_info = "len=" + str(len(record.seq))
        path_info = "path=[0:0-" + str(len(record.seq)-1) + "]"
        record.description += " " + len_info + " " + path_info

    with open(output_path, "w") as output_file:
        SeqIO.write(longest_isoforms.values(), output_file, "fasta")

input_path = "/Users/trystan.bessonnier/Desktop/DESeq2_Transcriptome/trinity_tzet_complete_transcriptome.Trinity.fasta"
output_path = "/Users/trystan.bessonnier/Desktop/DESeq2_Transcriptome/Transcriptome_whithout_isoforme.fasta"

filter_longest_isoforms(input_path, output_path)
