# Alleleome generation step III- Parsing the results and generating the amino acid mutations
import logging
from . import amino_acid_variants_parallel, codon_mutations_parallel
from itertools import repeat
from multiprocessing import Pool
import pandas as pd
import gc
from pathlib import Path


def generate_amino_acid_vars(gene_list, out_dir, aa_vars_path, p=1):
    gene_list_len = len(gene_list)
    chunksize = max(min(gene_list_len // p, 50), 1)
    logging.info(f"Parallel chunksize = {chunksize}")
    counter = 0


    with open(aa_vars_path, "w") as f:
        with Pool(p) as pool:
            for result in pool.imap_unordered(
                amino_acid_variants_parallel.generate_amino_acid_vars,
                zip(gene_list, repeat(out_dir)),
                chunksize=chunksize,
            ):
                logging.info(f"Processing AAV result of gene #{counter+1}/{gene_list_len}")
                if not result:
                    gene_list_len -= 1
                    continue
                df = pd.DataFrame(result)
                df.to_csv(f, header=(counter == 0), index=False)
                counter += 1

def codon_mut(gene_list, out_dir, codon_mut_path, p=1):
    gene_list_len = len(gene_list)
    chunksize = max(min(gene_list_len // p, 5), 1)
    logging.info(f"Parallel chunksize = {chunksize}")
    counter = 0

    with open(codon_mut_path, "w") as f:
        with Pool(p) as pool:
            for result in pool.imap_unordered(
                codon_mutations_parallel.codon_mut,
                zip(gene_list, repeat(out_dir)),
                chunksize=chunksize,
            ):
                logging.info(f"Processing CM result of gene #{counter+1}/{gene_list_len}")
                if not result:
                    gene_list_len -= 1
                    continue
                df = pd.DataFrame(result)
                df.to_csv(f, header=(counter == 0), index=False)
                del result
                del df
                counter += 1
                if (counter // chunksize) == 0:
                    f.flush()
                    gc.collect()
