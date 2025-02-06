#!/usr/bin/env python3
#-*- coding: utf-8 -*-
# ---------------------- #

##### SCRIPT ######

# ---------------------- #
import os
os.getcwd()
# dossier où sont les reads et la ref
os.chdir(path="/home/stagiaire/Documents/Projet_alignement_sequences/RecombCompet/data")

# Ajouter le chemin de minimap2 au PATH
os.environ["PATH"] += os.pathsep + "/home/stagiaire/minimap2"

# ---------------------- #

## library ##
import subprocess # pour gérér le mapping via minimap2
import pysam # pour convertir fichier sam en bam
import matplotlib.pyplot as plt # pour visualiser plot
# ---------------------- #


# Demander les préfixes des fichiers et le numéro à l'utilisateur.rice
ref_prefix = input("Entrez le préfixe du fichier de référence (ex. ref) : ")
read_prefix = input("Entrez le préfixe du fichier des reads (ex. Data) : ")
read_number = input("Entrez le numéro des reads (ex. 1-86) : ")

# Construire les noms de fichiers à partir des préfixes et du numéro
reference_genome = f"{ref_prefix}{read_number}.fa"
reads_file = f"{read_prefix}{read_number}.fastq"

# Créer un dossier pour stocker les résultats
result_dir = os.path.join("/home/stagiaire/Documents/Projet_alignement_sequences/RecombCompet", "res", f"res_{read_number}")
os.makedirs(result_dir, exist_ok=True)

# Chemins des fichiers de sortie
sam_file = os.path.join(result_dir, f"alignment_{read_number}.sam")
bam_file = os.path.join(result_dir, f"alignment_{read_number}.bam")
sorted_bam_file = os.path.join(result_dir, f"align_sorted_{read_number}.bam")


# Exécuter minimap2 pour effectuer le mapping avec le genome de reference et le fichier des reads
minimap2_command = [
    "minimap2", "-a", reference_genome, reads_file
]
with open(sam_file, "w") as sam_output:
    subprocess.run(minimap2_command, stdout=sam_output)



# Convertir un fichier SAM en BAM
with pysam.AlignmentFile(sam_file, "r") as sam: # ouverture en "read" pour lire le fichier
    with pysam.AlignmentFile(bam_file, "wb", header=sam.header) as bam:
        for read in sam:
            bam.write(read)
            

# Trier le fichier BAM
pysam.sort("-o", sorted_bam_file, bam_file)

# Indexer le fichier BAM trié
pysam.index(sorted_bam_file)

# Obtenir des statistiques sur le fichier BAM trié
stats_file_path = os.path.join(result_dir, f"align_sorted_{read_number}.bam_stats") # stockage fichier dans dossier res
flagstat_output = pysam.flagstat(sorted_bam_file)
with open(stats_file_path, "w") as stats_file:
    stats_file.write(flagstat_output)
    
# ------- --------------- #

# Demander à l'utilisateur.rice si les % de reads mappés de chaque catégorie
# doivent être calculés ? 
calc_percent = input("Calcul des pourcentages ? (oui/non) : ")
if calc_percent == "oui":

    # Seuils pour classer les tailles des reads
    thresholds = {
        "< 212 pb": 212,
        "< 474 pb": 474,
        "< 732 pb": 732,
        "< 982 pb": 982,
        "< 1232 pb": 1232,
        "< 1482 pb": 1482,
        ">= 1482 pb": float('inf')
    }

    # Initialiser un dictionnaire pour stocker les occurrences des tailles de reads
    read_size_counts = {bin_label: 0 for bin_label in thresholds.keys()}
     # les clefs (bin_label) sont les seuils définis précédemment
    # et toutes sont initialisées avec une valeur de comptage de 0

    # Ouvrir le fichier BAM trié et compter les tailles des reads selon les seuils définis
    with pysam.AlignmentFile(sorted_bam_file, "rb") as bam:
        for read in bam:
            if not read.is_unmapped:  # on considère pas ceux qui ne sont pas alignés
                read_length = read.query_length # longueur des reads
                for bin_label, upper_bound in thresholds.items(): # bin_label est la catégorie de taille (ex : "< 474 pb").
                # upper_bound est la valeur seuil correspondante (ex : 474 pb)
                    if read_length < upper_bound:
                        read_size_counts[bin_label] += 1 # si la longueur du read est < au seuil ; incrémente la valeur de comptage 
                        # du dico au seuil correspondant
                        break

    # Calculer le total des reads alignés (somme de tous)
    total_reads = sum(read_size_counts.values())

    # Calculer le nombre total de reads recombinants (tous sauf ceux dans le dernier seuil car ils n'ont pas recombiné)
    recombinant_reads = sum(count for bin_label, count in read_size_counts.items() if bin_label != ">= 1482 pb")

    # Calculer les pourcentages pour chaque catégorie (un dictionnaire)
    read_size_percentages = {bin_label: (count / total_reads) * 100 for bin_label, count in read_size_counts.items()}

    # Calculer les pourcentages pour les reads recombinants (sauf dernière catégorie) (un dictionnaire)
    recombinant_percentages = {
        bin_label: (count / recombinant_reads) * 100
        for bin_label, count in read_size_counts.items() if bin_label != ">= 1482 pb"
    }

    # Stocker les résultats dans un fichier texte
    output_file = os.path.join(result_dir, f"read_size_percentages_by_threshold_{read_number}.txt")
    with open(output_file, "w") as file: # ouverture en "write" pour écrire
        file.write("Taille du fragment\tNombre d'occurrences\tPourcentage total\tPourcentage recombinant\n")
        for size, count in read_size_counts.items():
            total_percentage = read_size_percentages[size]
            recombinant_percentage = recombinant_percentages.get(size, 0)
            file.write(f"{size}\t{count}\t{total_percentage:.2f}%\t{recombinant_percentage:.2f}%\n")

    # Afficher les statistiques dans la console
    print(f"Total reads recombinants : {recombinant_reads}")
    for size, count in read_size_counts.items():
        total_percentage = read_size_percentages[size]
        recombinant_percentage = recombinant_percentages.get(size, 0)
        print(f"Taille {size} : {count} occurrences ({total_percentage:.2f}% total, {recombinant_percentage:.2f}% recombinant)")

    # Afficher un histogramme
    plt.bar(read_size_percentages.keys(), read_size_percentages.values(), color='lightblue', edgecolor='black', label="Pourcentage total")
    plt.title("Distribution des tailles de reads alignés par seuils")
    plt.xlabel("Seuil de taille du fragment (pb)")
    plt.ylabel("Pourcentage de reads")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.legend()
    output_plot = os.path.join(result_dir, f"read_size_distribution_percentages_histogram_{read_number}.png") 
    plt.savefig(output_plot) # sauvegarder
    plt.show() # montrer le plot dans la console



