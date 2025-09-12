# RecombCompet
Building a genetic circuit for differentiation by recombination

Usage : 
```bash
python3 analyse_replicats.py -h
usage: analyse_replicats.py [-h] --ref REF --read READ [--rep REP] [--datadir DATADIR] [--resdir RESDIR]

Analyse attB: mapping, détection pics, % par intervalle, réplicats

options:
  -h, --help         show this help message and exit
  --ref REF          préfixe du fichier référence, sans extension (ex: ref2.12 -> ref2.12.fa)
  --read READ        préfixe des reads, sans extension (ex: data2.12 -> data2.12-1.fastq, data2.12-2.fastq...) Attention : Les réplicats
                     doivent suivre le format: prefix-1.fastq, prefix-2.fastq, etc.
  --rep REP          nombre de réplicats (default = 1)
  --datadir DATADIR  dossier contenant ref.fa et reads.fastq
  --resdir RESDIR    dossier pour sauvegarder les résultats
```
