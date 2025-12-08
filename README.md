# Recombination Analysis Tool

## Description

`analyse_replicats.py` is a Python script for analyzing recombination events in sequencing data by detecting breakpoints based on attB sites and calculating the percentage of recombinant vs non-recombinant reads.

The tool:
- Maps reads to a reference genome using minimap2
- Automatically detects attB sequences in the reference
- Identifies breakpoints at position attB+21 (cut site between "gg")
- Calculates the percentage of reads in each interval (non-recombinant vs recombinant categories)
- Supports multiple replicates with statistical analysis (mean ± standard deviation)
- Generates comprehensive visualizations and summary statistics

## Features

- Automatic breakpoint detection based on attB sequence  
- Multi-replicate support with standard deviation calculation  
- Visualization: histograms with colored intervals and barplots  
- Statistical summary exported to text files  
- Flexible directory structure

## Requirements

### Dependencies


# Core tools
minimap2
samtools

# Python packages
python >= 3.7
pysam
numpy
matplotlib
biopython
scipy


### Installation

```bash
# Install via conda (recommended)
conda create -n recomb_analysis python=3.9
conda activate recomb_analysis
conda install -c bioconda minimap2 samtools pysam biopython
conda install -c conda-forge numpy matplotlib scipy

# Or via pip
pip install pysam numpy matplotlib biopython scipy
```

## Usage

### Basic command

Usage : 
```bash
python3 analyse_replicats.py [-h] --ref REF --read READ [--rep REP] [--datadir DATADIR] [--resdir RESDIR]

Analyse attB: mapping, détection pics, % par intervalle, réplicats

options:
  -h, --help         show this help message and exit
  --ref REF          préfixe du fichier référence, sans extension (ex: ref2.12 -> ref2.12.fa)
  --read READ        préfixe des reads, sans extension (ex: data2.12 -> data2.12-1.fastq, data2.12-2.fastq...) Attention : Les réplicats
                     doivent suivre le format: prefix-1.fastq, prefix-2.fastq, etc.
  --rep REP          nombre de réplicats (default = 1)
  --datadir DATADIR  dossier contenant ref.fa et reads.fastq
  --resdir RESDIR    dossier pour sauvegarder les résultats
``

---

## Input Files Structure

### Directory organization

```
RecombCompet/
── data/
│   ── ref2.12.fa          # Reference genome
│   ── data2.12-1.fastq    # Replicate 1
│   ── data2.12-2.fastq    # Replicate 2
│   ── data2.12-3.fastq    # Replicate 3
└── res/                     # Results will be created here
```

### File naming convention

**Reference genome**: `<prefix>.fa`
- Example: `ref2.12.fa`

**Reads (replicates)**: `<prefix>-<N>.fastq`
- Example: `data2.12-1.fastq`, `data2.12-2.fastq`, `data2.12-3.fastq`
- N = replicate number (1, 2, 3, ...)

---

## Examples

### Single replicate analysis

```bash
python analyse_replicats.py \
    --ref ref2.12 \
    --read data2.12 \
    --rep 1
```

### Multiple replicates (3 replicates)

```bash
python analyse_replicats.py \
    --ref ref2.12 \
    --read data2.12 \
    --rep 3
```

### Custom directories

```bash
python analyse_replicats.py \
    --ref my_reference \
    --read my_reads \
    --rep 2 \
    --datadir /home/user/my_data \
    --resdir /home/user/my_results
```

---

## Output Files

All results are saved in `resdir/res_<reference_prefix>/`

### Per-replicate files

For each replicate N:

| File | Description |
|------|-------------|
| `alignment_N.sam` | SAM alignment file |
| `alignment_N.bam` | BAM alignment file |
| `align_sorted_N.bam` | Sorted and indexed BAM |
| `align_sorted_N.bam.bai` | BAM index |
| `reads_breakpoints_N.txt` | Detailed statistics per interval |
| `hist_breakpoints_N.png` | Histogram with colored intervals |
| `barplot_breakpoints_N.png` | Barplot of recombination percentages |

### Summary files (all replicates)

| File | Description |
|------|-------------|
| `barplot_breakpoints_mean_<ref>.png` | Mean percentages with error bars (SD) |
| `summary_statistics_<ref>.txt` | Statistical summary table |

---

## Output Interpretation

### Categories

- **Break0**: Non-recombinant reads (before first breakpoint)
- **Break1 to BreakN**: Recombinant reads (between successive breakpoints)
- **BreakN+1**: Reads after last breakpoint

### Example output (`reads_breakpoints_1.txt`)

```
# Breakpoint analysis: intervals based on attB+21 positions
# attB positions: [1000, 5000, 9000]
# Breakpoint positions (attB+21): [1021, 5021, 9021]
Category    Interval        Reads   Percentage
Break0      <1021          8500    85.00
Break1      1021-5021      800     8.00
Break2      5021-9021      500     5.00
Break3      >9021          200     2.00
```

### Visualizations

**Histogram** (`hist_breakpoints_N.png`):
- Grey bars: read distribution along reference
- Colored zones: intervals between breakpoints
- Green dotted lines: attB start positions
- Red solid lines: actual breakpoints (attB+21)
- Orange dashed lines: automatically detected peaks

**Barplot** (`barplot_breakpoints_N.png`):
- X-axis: recombination categories
- Y-axis: percentage of reads
- Error bars: standard deviation (if multiple replicates)

---

## Methodology

### Breakpoint Detection

1. **attB sequence identification**:
   ```
   attB = "gcccggatgatcctgacgacggagaccgccgtcgtcgacaagccggccga"
   ```
   The script searches for this 51 bp sequence in the reference genome.

2. **Breakpoint calculation**:
   ```
   Breakpoint position = attB_start + 21
   ```
   Position 21 corresponds to the cut site between "gg" nucleotides.

3. **Interval classification**:
   - Reads starting before the first breakpoint → non-recombinant
   - Reads starting between breakpoints → recombinant (different categories)

### Read Mapping

Uses **minimap2** with default parameters for accurate alignment:
```bash
minimap2 -a reference.fa reads.fastq > alignment.sam
```

### Statistics

For multiple replicates:
- **Mean**: average percentage across replicates
- **Standard deviation**: calculated with ddof=1 (Bessel's correction)



## Troubleshooting

### Common errors

**Error: "File not found"**
- Check that file names follow the correct convention
- Verify that `--datadir` points to the correct directory
- Ensure `.fa` and `.fastq` extensions are NOT included in `--ref` and `--read`

**Error: "No attB sequence found"**
- Verify that the attB sequence exists in your reference genome
- Check that the reference is in FASTA format

**Low alignment rate**
- Verify read quality
- Check that reads match the reference genome

### Performance

**Memory usage**: ~2-4 GB for typical datasets  
**Time**: ~2-10 minutes per replicate (depends on data size)

For large datasets:
- Consider downsampling reads
- Use more threads if minimap2 supports it (modify script if needed)

---

## License

This script is provided as-is for research purposes under the CClicence. 

---

## Contact

For questions or issues, please contact: 
- capucine.mayoud@umontpellier.fr
- anna-sophie.fiston-lavier@umontpellier.fr
---

## Version History

**v1.0** (2025)
- Initial release
- Automatic attB detection
- Multi-replicate support
- Statistical analysis and visualization
