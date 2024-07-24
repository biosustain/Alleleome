import pandas as pd
from Bio import AlignIO
from pathlib import Path
from collections import Counter
import logging
import gzip
from multiprocessing import Pool
from itertools import repeat


def find_dominant_aa(gene_list, out_dir, dominant_aa_path, p=1):
    logging.info("Starting: preplot: find_dominant_aa")
    out_dir = Path(out_dir)
    gene_list_len = len(gene_list)
    chunksize = max(min(gene_list_len // p, 50), 1)
    logging.info(f"Parallel chunksize = {chunksize}")
    counter = 0

    with open(dominant_aa_path, "w") as f:
        # Process each alignment file
        with Pool(p) as pool:
            for result in pool.imap_unordered(
                find_dominant_aa_h,
                zip(gene_list, repeat(out_dir)),
                chunksize=chunksize,
            ):
                logging.info(f"Processing FDA result of gene #{counter+1}/{gene_list_len}")
                if not result:
                    gene_list_len -= 1
                    continue
                df = pd.DataFrame(result)
                df.to_csv(f, header=(counter == 0), index=False)
                counter += 1

    # Convert the list of dictionaries to a DataFrame and save to CSV
    # consensus_dominant_aa_df = pd.DataFrame(consensus_data)
    # consensus_dominant_aa_df.to_csv(dominant_aa_path) #out_dir / 'final_core_consensus_dominant_aa_count_df.csv'
    logging.info("Finishing: preplot: find_dominant_aa")

def find_dominant_aa_h(args):
    gene, out_dir = args
    return find_dominant_aa_(gene, out_dir)

def find_dominant_aa_(gene, out_dir):
    out_dir = Path(out_dir)
    msa_file_name = out_dir / "output" / gene / f'mafft_amino_acid_{gene}.fasta.gz'
    # logging.info(msa_file_name)
    consensus_data = {"Gene": [], "AA_cons": [], "AA_pos": [], "AA_freq": [], "Sequence_type": []}
    with gzip.open(msa_file_name, "rt") as f:
        aln = AlignIO.read(f, 'fasta')
        for i in range(aln.get_alignment_length()):
            count = Counter(aln[:, i]).most_common(1)
            # logging.info(count)
            aa_pos = i + 1
            n = len(count)

            consensus_data["Gene"].extend([gene]*n)
            aa_cons, aa_freq = zip(*count)
            consensus_data["AA_cons"].extend(aa_cons)
            consensus_data["AA_freq"].extend(aa_freq)
            consensus_data["AA_pos"].extend([aa_pos]*n)
            consensus_data["Sequence_type"].extend(["Dominant"]*n)
    return consensus_data