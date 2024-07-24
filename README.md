

# Alleleome: Core- and Panalleleome Generation for Pan-Genomes of Species

## Introduction
"Alleleome" is a specialized package designed to explore and analyze natural sequence variations within the Open Reading Frames (ORFs) of alleles of core genes in a species' pan-genome, both at the amino acid and nucleotide levels. It identifies variants such as substitutions, insertions, and deletions through a series of steps:

1. Initial QCQA of sequences.
2. Building consensus for each gene's allele set.
3. Pairwise alignment of consensus sequences with individual alleles.
4. Identification and generation of amino acid variant datasets.
5. Analysis of synonymous and non-synonymous substitutions from codons and corresponding amino acid data.

The Alleleome workflow is specifically tailored to study the natural sequence variations in core genes of the pan-genome of a species, with an emphasis on variations at the amino acid and nucleotide level.

### Publication
This package is based on the Alleleome package created by Archana S. Harke ([anpache/Core-Alleleome](https://github.com/anpanche/Core-Alleleome)) For more detailed information, refer to their publication:
[Early Release on BioRxiv](https://www.biorxiv.org/content/biorxiv/early/2023/09/22/2023.09.22.558971.full.pdf)

## Table of Contents
- [Introduction](#introduction)
- [Getting Started](#getting-started)
- [Usage](#usage)
- [Features](#features)
- [Built With](#built-with)
- [License](#license)

## Getting Started

#### Prerequisites
- Alleleome is tested and confirmed for Linux systems with the Conda package manager.
- Requires Python version 3.10 or higher.
- For optimal performance, especially when processing a large dataset, such as 1400 core genes and their respective alleles across 3400 strains, a high RAM capacity is strongly recommended.

### Installation
#### Using GitHub Repository

1. Clone the repository (ensure Git LFS is set up as described above).

2. Navigate to the Alleleome directory:
   ```
   cd Alleleome
   ```
3. Activate the virtual environment as instructed above.

4. Install the package:
   ```
   pip install .
   ```

## Usage

### Modes
The package consists of 5 subprograms that should be run sequentially:

1. **prepare**: Collect information on the genomes and genes/loci, and perform QC/QA.
2. **fasta**: Use the collected information to create amino acid and nucleotide fasta files for alignment. At this steps it should be specified whether teh Core-alleleome should be anlyzed (default) or the Panalleleome (`--pan`).
3. **process**: Run alignments using blast and mafft.
4. **analyze**: Analyze and process the sequence alignments.
5. **preplot**: Generate files required for plotting the results on PanKB.

The computationally intensive **process** and **analyze** subprograms can (and should typically) be run in parallel by specifying the number of available cores using `-p <num>`.

### Running the Package
You can find the full usage and parameters of `alleleome` by using the `--help` function:
```bash
$ alleleome -h
usage: alleleome [-h] {prepare,fasta,process,analyze,preplot} ...

Alleleome - Explore and analyze natural sequence variations within the Open Reading Frames (ORFs) of alleles of genes in a species pan-genome.

positional arguments:
  {prepare,fasta,process,analyze,preplot}

options:
  -h, --help            show this help message and exit
```
The different modes should be run in the order: prepare > fasta > process > analyze > preplot. Info on the inputs and outputs for every mode can be found by invoking that mode with the --help argument. E.g.:
```bash
$ alleleome prepare -h
usage: alleleome prepare [-h] --gp_binary GP_BINARY --gp_locustag GP_LOCUSTAG --summary SUMMARY --summary_v2 SUMMARY_V2 --gbk_folder GBK_FOLDER --all_locustag ALL_LOCUSTAG --all_genes ALL_GENES --sel_locustag SEL_LOCUSTAG --sel_genes SEL_GENES

options:
  -h, --help            show this help message and exit
  --gp_binary GP_BINARY
                        Path to gene_presence_binary csv file.
  --gp_locustag GP_LOCUSTAG
                        Path to gene_presence_locustag csv file.
  --summary SUMMARY     Path to df_pangene_summary.csv file created by roary.
  --summary_v2 SUMMARY_V2
                        Path to the updated summary.
  --gbk_folder GBK_FOLDER
                        Folder containing GenBank files.
  --all_locustag ALL_LOCUSTAG
                        Path to all_locustags csv file.
  --all_genes ALL_GENES
                        Path to all_genes csv file.
  --sel_locustag SEL_LOCUSTAG
                        Path to sel_locustags csv file.
  --sel_genes SEL_GENES
                        Path to sel_genes csv file.
```

## Features
Alleleome introduces the concept of "ORF alleleome," encapsulating the gene alleles found across all strains of a species, thus providing a comprehensive view of genome-scale sequence variations. This analysis can be instrumental in understanding sequence diversity characteristics and natural selection processes across different species within a family. The study of the alleleome offers insights into the genetic basis of natural selection in a species.

Key features include:
- Analysis of sequence variants using the consensus sequence of ORFs.
- Identification of dominant amino acids and their variants at specific positions.
- Revealing natural sequence and structural variations compared with the consensus sequence and structural attributes.
- Identification of genome-scale synonymous and non-synonymous mutations through the analysis of codon changes and their corresponding amino acid changes."

## Built With
- Python and Biopython.
- Integrated with "BGCflow" workflow using SnakeMake.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
