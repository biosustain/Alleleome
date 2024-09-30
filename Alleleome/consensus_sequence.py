# Alleleome generation step I- Get the consensus of the allels (amino acid and nucleotide sequences) of each core gene
import logging
import os
from pathlib import Path
import gzip
import pandas as pd
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Align.Applications import MafftCommandline


def build_consensus(
    gene_list, out_dir, p=1
):
    logging.info("Starting build_consensus in build_consensus_sequence")
    out_dir = Path(out_dir)

    if p == 1:
        logging.info(f"Aligning sequences sequentially (p={p})")
        for gene_id in gene_list:
            build_single_gene_consensus(gene_id, out_dir)
    else:
        logging.info(f"Aligning sequences in parallel (p={p})")

        from multiprocessing import Pool
        from itertools import repeat

        chunksize = max(min(len(gene_list) // p, 500), 1)
        logging.info(f"Parallel chunksize = {chunksize}")
        with Pool(p) as pool:
            for _ in pool.imap_unordered(build_single_gene_consensus_parallel, zip(gene_list, repeat(out_dir)), chunksize=chunksize):
                pass
    
    logging.info("Completed build_consensus in build_consensus_sequence")

def build_single_gene_consensus_parallel(args):
    gene_id, out_dir = args
    build_single_gene_consensus(gene_id, out_dir)

def build_single_gene_consensus(gene_id, out_dir):
    out_dir = Path(out_dir)
    alignment_input_dir = out_dir / "input" / gene_id
    alignment_output_dir = out_dir / "output" / gene_id
    alignment_output_dir.mkdir(parents=True, exist_ok=True)

    for cons_type, ext in {"amino_acid": "faa", "nucleotide": "fna"}.items():
        generate_consensus(
            gene_id,
            cons_type,
            alignment_input_dir / f"pan_genes.{ext}",
            alignment_output_dir / f"mafft_{cons_type}_{gene_id}.fasta.gz",
            alignment_output_dir / f"{cons_type}_consensus_{gene_id}.{ext}.gz"
        )


def generate_consensus(gene_id, cons_type, input_path, mafft_output_path, consensus_output_path):
    if mafft_output_path.is_file() and consensus_output_path.is_file():
        logging.info(f"Outputs for {gene_id} already present, skipping.")
        return
    if os.stat(input_path).st_size == 0:
        logging.info(f"The input file is empty for {gene_id}")
        return

    mafft_cline = MafftCommandline(
            input=f"\"{input_path}\"",
            amino=(cons_type == "amino_acid"),
            nuc=(cons_type == "nucleotide"),
        )
    stdout, stderr = mafft_cline()

    with gzip.open(mafft_output_path, "wt") as handle:
        handle.write(stdout)

    with gzip.open(mafft_output_path, "rt") as handle:
        myalign = AlignIO.read(handle, "fasta")
    summary = AlignInfo.SummaryInfo(myalign)
    consensus = summary.dumb_consensus(threshold=0.5)
    seq = str(consensus).upper()

    with gzip.open(consensus_output_path, "wt") as consensus_file:
        consensus_seq = "".join([">" + gene_id + "\n" + seq])
        consensus_file.write(consensus_seq)
