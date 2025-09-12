#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import pysam
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from scipy.signal import find_peaks

def reads_by_intervals(boundaries, read_starts, total_reads, label):
    results = []
    intervals = []

    # avant le premier
    intervals.append((float("-inf"), boundaries[0]))
    # entre
    for i in range(len(boundaries)-1):
        intervals.append((boundaries[i], boundaries[i+1]))
    # après le dernier
    intervals.append((boundaries[-1], float("inf")))

    for i, (start, end) in enumerate(intervals):
        count = sum(1 for pos in read_starts if start <= pos < end)
        percent = count / total_reads * 100
        results.append((f"{label}{i}", start, end, count, percent))

    return results


def analyse_replicats(nomref, nomread, nbdereplicat=1,
                      data_dir="./RecombCompet/data",
                      res_dir="./RecombCompet/res"):
    """
    Analyse plusieurs réplicats et produit un barplot des % de recombinaison à chaque cassure, et si plusieurs réplicats, calcule l'écart-type
    """

    attB_seq = "gcccggatgatcctgacgacggagaccgccgtcgtcgacaagccggccga".upper()

    result_dir = os.path.join(res_dir, f"res_{nomref}")
    os.makedirs(result_dir, exist_ok=True)

    all_attB_percents = []

    for rep_idx in range(1, nbdereplicat+1):
        reads_file = os.path.join(data_dir, f"{nomread}-{rep_idx}.fastq")
        reference_genome = os.path.join(data_dir, f"{nomref}.fa")

        print(f"\n=== Réplica {rep_idx} : {reads_file} ===")

        sam_file = os.path.join(result_dir, f"alignment_{rep_idx}.sam")
        bam_file = os.path.join(result_dir, f"alignment_{rep_idx}.bam")
        sorted_bam_file = os.path.join(result_dir, f"align_sorted_{rep_idx}.bam")
        output_txt_file = os.path.join(result_dir, f"reads_attB_{rep_idx}.txt")
        hist_file = os.path.join(result_dir, f"hist_attB_{rep_idx}.png")
        barplot_file = os.path.join(result_dir, f"barplot_attB_{rep_idx}.png")

        # ----------------------
        # Détection attB
        attB_positions = []
        for record in SeqIO.parse(reference_genome, "fasta"):
            ref_seq = str(record.seq).upper()
            start = 0
            while True:
                idx = ref_seq.find(attB_seq, start)
                if idx == -1:
                    break
                attB_positions.append(idx)
                start = idx + 1
        attB_positions = sorted(attB_positions)
        print(f"Positions attB détectées : {attB_positions}")

        # ----------------------
        # Mapping avec minimap2
        print("Mapping reads...")
        with open(sam_file, "w") as sam_output:
            subprocess.run(["minimap2", "-a", reference_genome, reads_file], stdout=sam_output)

        # ----------------------
        # Conversion SAM -> BAM
        with pysam.AlignmentFile(sam_file, "r") as sam:
            with pysam.AlignmentFile(bam_file, "wb", header=sam.header) as bam:
                for read in sam:
                    bam.write(read)

        # Tri et indexation
        pysam.sort("-o", sorted_bam_file, bam_file)
        pysam.index(sorted_bam_file)

        # ----------------------
        # Extraction positions de départ
        read_starts = []
        total_reads = 0
        with pysam.AlignmentFile(sorted_bam_file, "rb") as bam:
            for read in bam.fetch():
                if read.is_unmapped:
                    continue
                total_reads += 1
                read_starts.append(read.reference_start)
        print(f"Total reads alignés : {total_reads}")

        # Détection pics
        counts, bin_edges = np.histogram(read_starts, bins=1000)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        peaks, _ = find_peaks(counts, prominence=10)
        peak_positions = sorted(bin_centers[peaks])
        print(f"Pics détectés : {peak_positions}")

        # ----------------------
        # Calcul % reads par intervalles
        results_attB = reads_by_intervals(attB_positions, read_starts, total_reads, "Seq")
        results_peaks = reads_by_intervals(peak_positions, read_starts, total_reads, "Pic")

        # Sauvegarde TXT
        with open(output_txt_file, "w") as ftxt:
            ftxt.write("Categorie\tIntervalle\tReads\tPourcentage\n")
            for res in results_attB + results_peaks:
                name, start, end, count, percent = res
                interval = f"{start:.0f}-{end:.0f}" if np.isfinite(start) and np.isfinite(end) else (f"<{end:.0f}" if start==-np.inf else f">{start:.0f}")
                ftxt.write(f"{name}\t{interval}\t{count}\t{percent:.2f}\n")

        all_attB_percents.append([res[4] for res in results_attB])

        # ----------------------
        # Histogramme
        plt.figure(figsize=(14,6))
        plt.hist(read_starts, bins=1000, color='lightgrey', edgecolor='black', alpha=0.6)
        colors_attB = plt.cm.viridis(np.linspace(0,1,len(results_attB)))
        for i, (name, start, end, _, percent) in enumerate(results_attB):
            if not np.isfinite(start):
                start = min(read_starts)
            if not np.isfinite(end):
                end = max(read_starts)
            plt.axvspan(start, end, color=colors_attB[i], alpha=0.3, label=f"{name}: {percent:.1f}%")
        for attB in attB_positions:
            plt.axvline(attB, color='green', linestyle=':', lw=2)
        for peak in peak_positions:
            plt.axvline(peak, color='red', linestyle='--', lw=1)
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.tight_layout()
        plt.savefig(hist_file, dpi=300)
        plt.close()

        # ----------------------
        # Barplot % reads attB
        means = all_attB_percents[-1]
        stds = [0]*len(means)
        plt.figure(figsize=(10,6))
        bars = plt.bar([res[0] for res in results_attB], means, yerr=stds, capsize=5, color=plt.cm.viridis(np.linspace(0,1,len(means))))
        plt.ylabel("% de reads")
        plt.xlabel("Catégories (Seq0-SeqN)")
        plt.title(f"Réplica {rep_idx}")
        plt.ylim(0, max(means)*1.2)
        for bar, mean in zip(bars, means):
            plt.text(bar.get_x()+bar.get_width()/2, bar.get_height()+1, f"{mean:.1f}%", ha='center', va='bottom')
        plt.tight_layout()
        plt.savefig(barplot_file, dpi=300)
        plt.close()

    # ----------------------
    # Barplot moyen ± écart-type sur réplicats
    all_attB_percents = np.array(all_attB_percents)
    mean_perc = np.mean(all_attB_percents, axis=0)
    std_perc = np.std(all_attB_percents, axis=0, ddof=1)
    categories = [res[0] for res in results_attB]

    plt.figure(figsize=(10,6))
    bars = plt.bar(categories, mean_perc, yerr=std_perc, capsize=5, color=plt.cm.viridis(np.linspace(0,1,len(categories))))
    plt.ylabel("% de reads")
    plt.xlabel("Catégories (Seq0-SeqN)")
    plt.title(f"Moyenne ± écart-type sur {nbdereplicat} réplicats")
    plt.ylim(0, max(mean_perc)*1.3)
    for bar, mean in zip(bars, mean_perc):
        plt.text(bar.get_x()+bar.get_width()/2, bar.get_height()+1, f"{mean:.1f}%", ha='center', va='bottom')
    plt.tight_layout()
    plt.savefig(os.path.join(result_dir, f"barplot_attB_mean_{nomref}.png"), dpi=300)
    plt.close()

    print(f"\n✅ Analyse terminée, résultats dans {result_dir}")

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyse attB: mapping, détection pics, % par intervalle, réplicats")
    parser.add_argument("--ref", required=True, help="préfixe du fichier référence, sans extension (ex: ref2.12 -> ref2.12.fa)")
    parser.add_argument("--read", required=True, help="préfixe des reads, sans extension (ex: data2.12 -> data2.12-1.fastq, data2.12-2.fastq...)\n"
             "Attention : Les réplicats doivent suivre le format: prefix-1.fastq, prefix-2.fastq, etc.")
    parser.add_argument("--rep", type=int, default=1, help="nombre de réplicats (default = 1)")
    parser.add_argument("--datadir", type=str, default="./RecombCompet/data",
                        help="dossier contenant ref.fa et reads.fastq")
    parser.add_argument("--resdir", type=str, default="./RecombCompet/res",
                        help="dossier pour sauvegarder les résultats")
    args = parser.parse_args()

    analyse_replicats(args.ref, args.read, nbdereplicat=args.rep,
                      data_dir=args.datadir, res_dir=args.resdir)