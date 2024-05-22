# Installer les packages nécessaires si ce n'est pas déjà fait
install.packages("ggplot2")
library(ggplot2)

setwd("/Users/trystan.bessonnier/Desktop/Bioinfo/QC_transcriptome")

# Lire les fichiers CSV
new_data_transcript <- read.csv("new_data_transcript.csv")
old_data_transcript <- read.csv("old_data_transcript.csv")

# Ajouter une colonne pour différencier les deux ensembles de données
new_data_transcript$type <- "Nouveau"
old_data_transcript$type <- "Ancien"

# Combiner les deux ensembles de données
combined_data <- rbind(new_data_transcript, old_data_transcript)

# Tracer le graphique
ggplot(combined_data, aes(x = Ex, y = ExN50, color = type)) +
  geom_line() +
  geom_point() +
  labs(title = "Comparaison des EX90N50 transcriptomes",
       x = "Ex",
       y = "ExN50",
       color = "Type de transcriptome") +
  theme_minimal()

# Lire les fichiers CSV
new_data_genes <- read.csv("new_data_genes.csv")
old_data_genes <- read.csv("old_data_genes.csv")

# Ajouter une colonne pour différencier les deux ensembles de données
new_data_genes$type <- "Nouveau"
old_data_genes$type <- "Ancien"

# Combiner les deux ensembles de données
combined_data <- rbind(new_data_genes, old_data_genes)

# Tracer le graphique
ggplot(combined_data, aes(x = Ex, y = ExN50, color = type)) +
  geom_line() +
  geom_point() +
  labs(title = "Comparaison des EX90N50 transcriptomes",
       x = "Ex",
       y = "ExN50",
       color = "Type de transcriptome") +
  theme_minimal()



