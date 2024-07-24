from pathlib import Path
import json
from Bio import SeqIO
from collections import defaultdict
import logging
import gzip

def calculate_msa_freq(gene_list, out_dir, aa_freq_dir):
    logging.info("Starting: preplot: calculate_msa_freq")
    out_dir = Path(out_dir)
    aa_freq_dir = Path(aa_freq_dir)

    # start generating required files
    for gene in gene_list:
        mafft_out_file = out_dir / "output" / gene / f"mafft_amino_acid_{gene}.fasta.gz"
        mafft_new_file = aa_freq_dir / gene / "MSA.fasta"

        # Calculate the AA frequency          
        # Define the standard set of amino acids
        amino_acids = "ACDEFGHIKLMNPQRSTVWY"

        # Read the fasta file
        with gzip.open(mafft_out_file, "rb") as f_in:
            with open(mafft_new_file, "wb") as f_out:
                data = f_in.read()
                f_out.write(data)
        with open(mafft_new_file, "r") as f:
            sequences = [str(record.seq) for record in SeqIO.parse(f, 'fasta')]

        # Assume the sequences are already aligned
        assert all(len(s) == len(sequences[0]) for s in sequences)

        # Use a list of defaultdicts to track counts of amino acids at each position
        counts = [defaultdict(int) for _ in range(len(sequences[0]))]

        # Update counts
        for sequence in sequences:
            for i, amino_acid in enumerate(sequence):
                counts[i][amino_acid] += 1

        # # Calculate percentage
        # percentages = [{amino_acid: count / len(sequences) * 100 for amino_acid, count in pos_counts.items()} for pos_counts in counts]
        n_sequences = len(sequences)
        # Prepare data for JSON
        data_for_json = {}
        for i, pos_counts in enumerate(counts):
            data_for_json[i+1] = {}
            for amino_acid, count in pos_counts.items():
                if count > 0:
                    data_for_json[i+1][amino_acid] = [count, round(100.0*count/n_sequences, 2)]

        aa_freq_path = aa_freq_dir / gene
        aa_freq_path.mkdir(parents=True, exist_ok=True)
        aa_freq_path = aa_freq_path / "AA_freq.json"
        # Write to JSON file
        with open(aa_freq_path, 'w') as f:
            json.dump(data_for_json, f)
    logging.info("Finishing: preplot: calculate_msa_freq")