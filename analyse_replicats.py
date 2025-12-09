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
    """
    Calculate the number and percentage of reads in intervals defined by boundaries.
    
    Args:
        boundaries: List of positions defining interval boundaries (breakpoint positions)
        read_starts: List of read start positions
        total_reads: Total number of aligned reads
        label: Prefix label for interval names
    
    Returns:
        List of tuples (name, start, end, count, percentage)
    """
    results = []
    intervals = []

    # Before first breakpoint (non-recombinant)
    intervals.append((float("-inf"), boundaries[0]))
    # Between breakpoints (recombinant categories)
    for i in range(len(boundaries)-1):
        intervals.append((boundaries[i], boundaries[i+1]))
    # After last breakpoint
    intervals.append((boundaries[-1], float("inf")))

    for i, (start, end) in enumerate(intervals):
        count = sum(1 for pos in read_starts if start <= pos < end)
        percent = count / total_reads * 100
        results.append((f"{label}{i}", start, end, count, percent))

    return results


def analyse_replicats(nomref, nomread, nbdereplicat=1, 
                      data_dir="./RecombCompet/data",
                      res_dir="./RecombCompet/res", 
                      attb_file=None):
    """
    Analyze multiple replicates and produce barplots of recombination percentages at each breakpoint.
    If multiple replicates, calculates standard deviation.
    
    Breakpoints are automatically detected as attB_position + 21 nucleotides (cuts between "gg").
    
    Args:
        nomref: Reference genome file prefix (without .fa extension)
        nomread: Reads file prefix (without -N.fastq)
        nbdereplicat: Number of replicates
        data_dir: Directory containing reference and read files
        res_dir: Directory for saving results
        attb_file: attB sequence (file .txt)
    """

    # ----------------------
    # Load attB sequence (user file or default)
    # ----------------------
    if attb_file is None:
        default_attb_path = os.path.join(data_dir, "attB_default.txt")
        print(f"Using default attB file: {default_attb_path}")
        with open(default_attb_path) as f:
            attB_seq = f.read().strip().upper()
    else:
        print(f"Using user-provided attB file: {attb_file}")
        with open(attb_file) as f:
            attB_seq = f.read().strip().upper()

    breakpoint_offset = 21

    result_dir = os.path.join(res_dir, f"res_{nomref}")
    os.makedirs(result_dir, exist_ok=True)

    # Store percentages from each replicate for calculating mean and std
    all_breakpoint_percents = []

    for rep_idx in range(1, nbdereplicat + 1):
        reads_file = os.path.join(data_dir, f"{nomread}-{rep_idx}.fastq")
        reference_genome = os.path.join(data_dir, f"{nomref}.fa")

        print(f"\n=== Replicate {rep_idx}: {reads_file} ===")

        sam_file = os.path.join(result_dir, f"alignment_{rep_idx}.sam")
        bam_file = os.path.join(result_dir, f"alignment_{rep_idx}.bam")
        sorted_bam_file = os.path.join(result_dir, f"align_sorted_{rep_idx}.bam")
        output_txt_file = os.path.join(result_dir, f"reads_breakpoints_{rep_idx}.txt")
        hist_file = os.path.join(result_dir, f"hist_breakpoints_{rep_idx}.png")
        barplot_file = os.path.join(result_dir, f"barplot_breakpoints_{rep_idx}.png")

        # ----------------------
        # Detect attB positions and calculate breakpoint positions
        # ----------------------
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

        # Calculate breakpoint positions (attB + 21 nucleotides)
        breakpoint_positions = [pos + breakpoint_offset for pos in attB_positions]

        print(f"attB positions detected: {attB_positions}")
        print(f"Breakpoint positions (attB+{breakpoint_offset}): {breakpoint_positions}")

        # ----------------------
        # Map reads with minimap2
        # ----------------------
        print("Mapping reads...")
        with open(sam_file, "w") as sam_output:
            subprocess.run(["minimap2", "-a", reference_genome, reads_file], stdout=sam_output)

        # ----------------------
        # Convert SAM to BAM
        # ----------------------
        with pysam.AlignmentFile(sam_file, "r") as sam:
            with pysam.AlignmentFile(bam_file, "wb", header=sam.header) as bam:
                for read in sam:
                    bam.write(read)

        # Sort and index BAM file
        pysam.sort("-o", sorted_bam_file, bam_file)
        pysam.index(sorted_bam_file)

        # ----------------------
        # Extract read start positions
        # ----------------------
        read_starts = []
        total_reads = 0
        with pysam.AlignmentFile(sorted_bam_file, "rb") as bam:
            for read in bam.fetch():
                if read.is_unmapped:
                    continue
                total_reads += 1
                read_starts.append(read.reference_start)

        print(f"Total aligned reads: {total_reads}")

        # ----------------------
        # Detect peaks for visualization (optional)
        # ----------------------
        counts, bin_edges = np.histogram(read_starts, bins=1000)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        peaks, _ = find_peaks(counts, prominence=10)
        peak_positions = sorted(bin_centers[peaks])

        print(f"Peaks detected (for visualization): {peak_positions}")

        # ----------------------
        # Calculate % of reads per interval (based on BREAKPOINTS)
        # ----------------------
        results_breakpoints = reads_by_intervals(
            breakpoint_positions, read_starts, total_reads, "Break"
        )

        # Save results to text file
        with open(output_txt_file, "w") as ftxt:
            ftxt.write("# Breakpoint analysis: intervals based on attB+21 positions\n")
            ftxt.write(f"# attB positions: {attB_positions}\n")
            ftxt.write(f"# Breakpoint positions (attB+{breakpoint_offset}): {breakpoint_positions}\n")
            ftxt.write("Category\tInterval\tReads\tPercentage\n")
            for res in results_breakpoints:
                name, start, end, count, percent = res
                if np.isfinite(start) and np.isfinite(end):
                    interval = f"{start:.0f}-{end:.0f}"
                elif start == -np.inf:
                    interval = f"<{end:.0f}"
                else:
                    interval = f">{start:.0f}"
                ftxt.write(f"{name}\t{interval}\t{count}\t{percent:.2f}\n")

        # Store percentages for this replicate
        all_breakpoint_percents.append([res[4] for res in results_breakpoints])

        # ----------------------
        # Generate histogram with breakpoint intervals
        # ----------------------
        plt.figure(figsize=(14, 6))
        plt.hist(read_starts, bins=1000, color='lightgrey', edgecolor='black', alpha=0.6)

        colors_breaks = plt.cm.viridis(np.linspace(0, 1, len(results_breakpoints)))
        for i, (name, start, end, _, percent) in enumerate(results_breakpoints):
            s = min(read_starts) if not np.isfinite(start) else start
            e = max(read_starts) if not np.isfinite(end) else end
            plt.axvspan(s, e, color=colors_breaks[i], alpha=0.3, label=f"{name}: {percent:.1f}%")

        for attB in attB_positions:
            plt.axvline(attB, color='green', linestyle=':', lw=2, alpha=0.5)

        for bp in breakpoint_positions:
            plt.axvline(bp, color='red', linestyle='-', lw=2)

        for peak in peak_positions:
            plt.axvline(peak, color='orange', linestyle='--', lw=1, alpha=0.5)

        plt.xlabel("Position on reference genome")
        plt.ylabel("Number of reads")
        plt.title(f"Read distribution with breakpoints - Replicate {rep_idx}")
        plt.tight_layout()
        plt.savefig(hist_file, dpi=300)
        plt.close()

        # ----------------------
        # Generate barplot for this replicate
        # ----------------------
        means = all_breakpoint_percents[-1]
        stds = [0] * len(means)

        plt.figure(figsize=(10, 6))
        bars = plt.bar([res[0] for res in results_breakpoints], means, yerr=stds, capsize=5,
                       color=plt.cm.viridis(np.linspace(0, 1, len(means))))
        plt.ylabel("% of reads")
        plt.xlabel("Categories (Break0=non-recombinant, Break1-N=recombinant)")
        plt.title(f"Recombination at breakpoints - Replicate {rep_idx}")
        plt.ylim(0, max(means)*1.2 if max(means) > 0 else 10)

        for bar, mean in zip(bars, means):
            plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                     f"{mean:.1f}%", ha='center', va='bottom', fontsize=9)

        plt.tight_layout()
        plt.savefig(barplot_file, dpi=300)
        plt.close()

    # ----------------------
    # Mean barplot across replicates
    # ----------------------
    all_breakpoint_percents = np.array(all_breakpoint_percents)
    mean_perc = np.mean(all_breakpoint_percents, axis=0)
    std_perc = np.std(all_breakpoint_percents, axis=0, ddof=1) if nbdereplicat > 1 else np.zeros_like(mean_perc)
    categories = [res[0] for res in results_breakpoints]

    plt.figure(figsize=(10, 6))
    bars = plt.bar(categories, mean_perc, yerr=std_perc, capsize=5,
                   color=plt.cm.viridis(np.linspace(0, 1, len(categories))))
    plt.ylabel("% of reads")
    plt.xlabel("Categories (Break0=non-recombinant, Break1-N=recombinant)")
    plt.title(f"Mean ± standard deviation over {nbdereplicat} replicate(s)")
    plt.ylim(0, max(mean_perc)*1.3 if max(mean_perc) > 0 else 10)

    for bar, mean, std in zip(bars, mean_perc, std_perc):
        label = f"{mean:.1f}%" if std == 0 else f"{mean:.1f}±{std:.1f}%"
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                 label, ha='center', va='bottom', fontsize=9)

    plt.tight_layout()
    plt.savefig(os.path.join(result_dir, f"barplot_breakpoints_mean_{nomref}.png"), dpi=300)
    plt.close()

    # ----------------------
    # Save summary statistics
    # ----------------------
    summary_file = os.path.join(result_dir, f"summary_statistics_{nomref}.txt")
    with open(summary_file, "w") as f:
        f.write(f"Summary Statistics - {nbdereplicat} replicate(s)\n")
        f.write("=" * 60 + "\n")
        f.write(f"Breakpoint offset: attB + {breakpoint_offset} nucleotides\n")
        f.write(f"attB positions: {attB_positions}\n")
        f.write(f"Breakpoint positions: {breakpoint_positions}\n")
        f.write("=" * 60 + "\n\n")
        f.write("Category\tMean(%)\tStd(%)\n")
        for cat, mean, std in zip(categories, mean_perc, std_perc):
            f.write(f"{cat}\t{mean:.2f}\t{std:.2f}\n")
        f.write("\n")
        f.write("Note: Break0 = non-recombinant reads (before first breakpoint)\n")
        f.write("      Break1-N = recombinant reads (between breakpoints)\n")

    print(f" Analysis complete! Results saved in {result_dir}")
    print(f"   - Breakpoint positions (attB+{breakpoint_offset}): {breakpoint_positions}")
    print(f"   - Individual replicate plots: barplot_breakpoints_<N>.png")
    print(f"   - Mean plot with SD: barplot_breakpoints_mean_{nomref}.png")
    print(f"   - Summary statistics: summary_statistics_{nomref}.txt")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Recombination analysis: mapping, breakpoint detection (attB+21), percentage calculation per interval, replicates"
    )
    parser.add_argument("--ref", required=True, 
                       help="Reference file prefix without extension (e.g., ref2.12 -> ref2.12.fa)")
    parser.add_argument("--read", required=True, 
                       help="Reads file prefix without extension (e.g., data2.12 -> data2.12-1.fastq, data2.12-2.fastq...)\n"
                            "Note: Replicates must follow format: prefix-1.fastq, prefix-2.fastq, etc.")
    parser.add_argument("--rep", type=int, default=1, 
                       help="Number of replicates (default = 1)")
    parser.add_argument("--datadir", type=str, default="./RecombCompet/data",
                       help="Directory containing ref.fa and reads.fastq files")
    parser.add_argument("--resdir", type=str, default="./RecombCompet/res",
                       help="Directory for saving results")
    parser.add_argument("--attb", type=str, default=None,
                    help="Optional file containing attB sequence. If not provided, default attB_default.txt is used.")

    
    args = parser.parse_args()

analyse_replicats(args.ref, args.read, nbdereplicat=args.rep,
                  data_dir=args.datadir, res_dir=args.resdir,
                  attb_file=args.attb)
