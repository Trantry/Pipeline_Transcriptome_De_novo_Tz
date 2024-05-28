# Chemin du fichier d'entrée
# file_path = '/Users/trystan.bessonnier/Desktop/Blast/blast_output_sequence_overexpressed_shrink_LFC_gene_926_ids.txt'

# Initialisation du dictionnaire pour stocker les meilleurs hits
# best_hits = {}

# Lecture du fichier ligne par ligne
# with open(file_path, 'r') as file:
#     for line in file:
#         # Séparer les colonnes par tabulation
#         columns = line.strip().split('\t')
#         gene = columns[0]
#         refseq = columns[1]
#         aligned_nucleotides = int(columns[3])
        
#         # Si le gène n'est pas encore dans le dictionnaire, ou si le nombre de nucléotides alignés est plus grand que le précédent, on met à jour
#         if gene not in best_hits or aligned_nucleotides > best_hits[gene][1]:
#             best_hits[gene] = (refseq, aligned_nucleotides, line)

# Chemin du fichier de sortie
# output_file_path = '/Users/trystan.bessonnier/Desktop/Blast/best_hits.txt'

# Écriture des résultats dans un nouveau fichier
# with open(output_file_path, 'w') as output_file:
#     for gene, (refseq, aligned_nucleotides, line) in best_hits.items():
#         output_file.write(line)





#python3 -m venv myenv

#source myenv/bin/activate

#pip install --upgrade pip

#pip install biopython



import re

# Chemin du fichier FASTA d'entrée
fasta_path = '/Users/trystan.bessonnier/Desktop/Blast/retrieved_sequence_seqids_Protein_overexpressed_shrink_LFC_gene_926_ids.fasta'
# Chemin du fichier FASTA nettoyé
cleaned_fasta_path = '/Users/trystan.bessonnier/Desktop/Blast/cleaned_retrieved_sequence_seqids_Protein_overexpressed_shrink_LFC_gene_926_ids.fasta'

# Lire le fichier FASTA et nettoyer les en-têtes multiples
with open(fasta_path, 'r') as infile, open(cleaned_fasta_path, 'w') as outfile:
    for line in infile:
        if line.startswith('>'):
            # Séparer les en-têtes multiples et écrire chaque en-tête sur une nouvelle ligne
            headers = line.split('>')
            for header in headers:
                if header.strip():
                    outfile.write(f'>{header.strip()}\n')
        else:
            outfile.write(line)

print(f'Le fichier FASTA nettoyé a été créé : {cleaned_fasta_path}')


from Bio import SeqIO

# Chemin du fichier des meilleurs hits
best_hits_path = '/Users/trystan.bessonnier/Desktop/Blast/best_hits_926.txt'

# Chemin du fichier FASTA nettoyé
cleaned_fasta_path = '/Users/trystan.bessonnier/Desktop/Blast/cleaned_retrieved_sequence_seqids_Protein_overexpressed_shrink_LFC_gene_926_ids.fasta'

# Chemin du fichier FASTA de sortie
output_fasta_path = '/Users/trystan.bessonnier/Desktop/Blast/best_hits_sequences_926.fasta'

# Lire les meilleurs hits et stocker les RefSeq IDs
refseq_ids = set()
with open(best_hits_path, 'r') as file:
    for line in file:
        columns = line.strip().split('\t')
        refseq_id = columns[1]
        refseq_ids.add(refseq_id)

print(f'Nombre d\'IDs uniques dans best_hits_926.txt : {len(refseq_ids)}')

# Lire le fichier FASTA nettoyé et extraire les séquences correspondant aux RefSeq IDs
sequences_to_write = []
with open(cleaned_fasta_path, 'r') as fasta_file:
    for record in SeqIO.parse(fasta_file, 'fasta'):
        # Vérifier si l'entête contient des informations autres que ">"
        if any(char.isalnum() for char in record.id) and record.seq:
            if record.id.strip() in refseq_ids:
                sequences_to_write.append(record)

print(f'Nombre de séquences uniques trouvées dans le fichier FASTA : {len(sequences_to_write)}')

# Écrire les séquences sélectionnées dans un nouveau fichier FASTA
with open(output_fasta_path, 'w') as output_fasta_file:
    SeqIO.write(sequences_to_write, output_fasta_file, 'fasta')

print(f'Le fichier FASTA avec les meilleures séquences a été créé : {output_fasta_path}')


# Chemin du fichier FASTA
fasta_path = '/Users/trystan.bessonnier/Desktop/Blast/best_hits_sequences_926.fasta'

# Chemin du fichier de sortie pour les noms d'entête
output_header_file_path = '/Users/trystan.bessonnier/Desktop/Blast/headers_in_best_hits.fasta'

# Lire le fichier FASTA et extraire les noms d'entête
header_set = set()
with open(fasta_path, 'r') as fasta_file, open(output_header_file_path, 'w') as output_file:
    for line in fasta_file:
        if line.startswith('>'):
            header = line.strip()[1:]  # Supprimer le caractère '>'
            header_set.add(header)

    # Écrire les noms d'entête dans le fichier de sortie
    for header in header_set:
        output_file.write(">" + header + '\n')

print(f'Les noms d\'entête ont été écrits dans le fichier : {output_header_file_path}')





from collections import Counter

# Chemin du fichier FASTA
fasta_path = '/Users/trystan.bessonnier/Desktop/Blast/best_hits_sequences_926.fasta'

# Fonction pour extraire les identifiants répétés
def find_duplicate_ids(fasta_file):
    # Compteur des identifiants
    id_counter = Counter()
    # Liste pour stocker les identifiants répétés
    duplicate_ids = []
    # Parcourir les enregistrements FASTA
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                # Extraire l'identifiant en supprimant le caractère '>'
                seq_id = line.strip()[1:]
                # Ajouter l'identifiant à compter
                id_counter[seq_id] += 1
    # Filtrer les identifiants répétés
    for id, count in id_counter.items():
        if count > 1:
            duplicate_ids.append(id)
    return duplicate_ids

# Appel de la fonction pour trouver les identifiants répétés
duplicate_ids = find_duplicate_ids(fasta_path)

# Chemin du fichier de sortie
output_file_path = '/Users/trystan.bessonnier/Desktop/Blast/duplicate_ids_in_best_hits.fasta'

# Écrire les identifiants répétés dans le fichier de sortie
with open(output_file_path, 'w') as output_file:
    output_file.write("Identifiants répétés dans le fichier best_hits_sequences_926.fasta :\n")
    for id in duplicate_ids:
        output_file.write(">" + id + '\n')

print(f'Les identifiants répétés ont été écrits dans le fichier : {output_file_path}')



# from collections import Counter

# Chemin du fichier FASTA
#fasta_path = '/Users/trystan.bessonnier/Desktop/Blast/best_hits_sequences_926.fasta'

# Chemin du fichier de sortie pour les identifiants d'entête tripliqués
#output_triplicate_headers_file_path = '/Users/trystan.bessonnier/Desktop/Blast/triplicate_headers_in_best_hits.fasta'

# Lire le fichier FASTA et compter le nombre d'occurrences de chaque identifiant d'entête
# header_counter = Counter()
# with open(fasta_path, 'r') as fasta_file:
#     for line in fasta_file:
#         if line.startswith('>'):
#             header = line.strip()[1:]  # Supprimer le caractère '>'
#             header_counter[header] += 1

# Filtrer les identifiants d'entête qui apparaissent trois fois
# triplicate_headers = [header for header, count in header_counter.items() if count == 3]

# Écrire les identifiants d'entête tripliqués dans le fichier de sortie
# with open(output_triplicate_headers_file_path, 'w') as output_file:
#     for header in triplicate_headers:
#         output_file.write(">" + header + '\n')

# print(f'Les identifiants d\'entête tripliqués ont été écrits dans le fichier : {output_triplicate_headers_file_path}')


# from collections import Counter

# Chemin du fichier FASTA
#fasta_path = '/Users/trystan.bessonnier/Desktop/Blast/best_hits_sequences_926.fasta'

# Chemin du fichier de sortie pour les identifiants d'entête quadruplés
#output_quadruplicate_headers_file_path = '/Users/trystan.bessonnier/Desktop/Blast/quadruplicate_headers_in_best_hits.fasta'

# Lire le fichier FASTA et compter le nombre d'occurrences de chaque identifiant d'entête
# header_counter = Counter()
# with open(fasta_path, 'r') as fasta_file:
#     for line in fasta_file:
#         if line.startswith('>'):
#             header = line.strip()[1:]  # Supprimer le caractère '>'
#             header_counter[header] += 1

# Filtrer les identifiants d'entête qui apparaissent quatre fois
# quadruplicate_headers = [header for header, count in header_counter.items() if count == 4]

# Écrire les identifiants d'entête quadruplés dans le fichier de sortie
# with open(output_quadruplicate_headers_file_path, 'w') as output_file:
#     for header in quadruplicate_headers:
#         output_file.write(">" + header + '\n')

# print(f'Les identifiants d\'entête quadruplés ont été écrits dans le fichier : {output_quadruplicate_headers_file_path}')




from collections import defaultdict

# Chemin du fichier FASTA
fasta_path = '/Users/trystan.bessonnier/Desktop/Blast/best_hits_sequences_926.fasta'

# Chemin du fichier de sortie pour les identifiants d'entête avec leur nombre d'occurrences
output_duplicates_count_file_path = '/Users/trystan.bessonnier/Desktop/Blast/headers_with_duplicates_count.txt'

# Initialiser un dictionnaire pour stocker les entêtes avec leur nombre d'occurrences
header_counts = defaultdict(int)
total_proteins = 0  # Initialiser le compteur total de protéines

# Lire le fichier FASTA et compter le nombre d'occurrences de chaque identifiant d'entête
with open(fasta_path, 'r') as fasta_file:
    current_header = ""
    for line in fasta_file:
        if line.startswith('>'):
            current_header = line.strip()[1:]  # Supprimer le caractère '>'
            header_counts[current_header] += 1
            total_proteins += 1  # Incrémenter le compteur total de protéines

# Écrire les identifiants d'entête avec leur nombre d'occurrences dans le fichier de sortie
with open(output_duplicates_count_file_path, 'w') as output_file:
    for header, count in header_counts.items():
        output_file.write(f">{header} : {count}\n")

# Écrire le nombre total de protéines dans un message
print(f'Le nombre total de protéines dans le fichier FASTA est : {total_proteins}')

print(f'Les identifiants d\'entête avec leur nombre d\'occurrences et leurs entêtes ont été écrits dans le fichier : {output_duplicates_count_file_path}')
