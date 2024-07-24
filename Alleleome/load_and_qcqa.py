import argparse
import logging
from pathlib import Path

import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def remove_special_char(s):
    if "/" in str(s):
        return s.replace("/", "_")
    elif "'" in str(s):
        return s.replace("'", "_variant")
    elif "(" in str(s):
        s = s.replace("(", "_")
        return s.replace(")", "")
    else:
        return s

def load_roary_data(pangene_summary_path, gene_presence_binary_path, gene_presence_locustag_path):
    """
    Load data from Roary output.

    Parameters:
    roary_path (Path): Path to Roary output

    Returns:
    tuple: A tuple containing two pandas DataFrames
    """
    df_pangene_summary = pd.read_csv(
        pangene_summary_path,
        low_memory=False,
        index_col=0
    )
    df_gene_presence_binary = pd.read_csv(
        gene_presence_binary_path, index_col="Gene", low_memory=False
    )
    df_gene_presence_locustag = pd.read_csv(
        gene_presence_locustag_path, index_col="Gene", low_memory=False, dtype=str
    )
    df_gene_presence_locustag.replace('', np.nan, inplace=True)
    df_gene_presence_locustag.index = [
        remove_special_char(str(i)) for i in list(df_gene_presence_locustag.index)
    ]
    return df_pangene_summary, df_gene_presence_binary, df_gene_presence_locustag

def load(path):
    df = pd.read_csv(
        path, index_col=0, low_memory=False
    )
    return df

def save(df, path):
    return df.to_csv(path)

def parse_genbank_files(df_gene_presence_locustag, gbk_folder):
    """
    Parse GenBank files.

    Parameters:
    df_gene_presence_locustag (DataFrame): DataFrame with gene presence information
    gbk_folder (Path): Path to folder containing GenBank files

    Returns:
    DataFrame: DataFrame with parsed GenBank data
    """

    all_locustag_list = []
    genome_ids = list(df_gene_presence_locustag.columns)
    for i, genome_id in enumerate(genome_ids):
        logging.info(f"Writing genome #{i+1}/{len(genome_ids)} to fasta.")
        cur_df = df_gene_presence_locustag[genome_id]
        cur_df = cur_df[~cur_df.isna()]
        cur_df = pd.Series(cur_df.index.values, index=cur_df)
        genbank_file_path = Path(gbk_folder) / f"{genome_id}.gbk"
        for record in SeqIO.parse(genbank_file_path, "genbank"):
            for feature in record.features:
                genome_data_list = []
                tag = feature.qualifiers.get("locus_tag")
                if not tag:
                    continue
                try:
                    gene_id = cur_df[tag[0]]
                # if len(gene_id) != 1:
                except:
                    continue
                genome_data_list.append(tag[0])  # Locus tag
                genome_data_list.append(genome_id)  # Genome ID
                genome_data_list.append(gene_id)
                
                if "product" in feature.qualifiers.keys():
                    genome_data_list.append(
                        feature.qualifiers["product"][0]
                    )  # Prokka annotation
                else:
                    genome_data_list.append("")
                genome_data_list.append(
                    int(feature.location.start)
                )  # Start position
                genome_data_list.append(int(feature.location.end))  # End position
                genome_data_list.append(len(feature))
                all_locustag_list.append(genome_data_list)

    all_locustag_df = pd.DataFrame(
        all_locustag_list,
        columns=[
            "Locus_Tag",
            "Genome_ID",
            "Gene_ID",
            "Prokka_Annotation",
            "Start_Position",
            "End_Position",
            "Nucleotide_Len",
        ],
    )
    # Use genome + locus tag as index, since there may be duplicate locus tags accross genomes.
    all_locustag_df.index = all_locustag_df["Genome_ID"] + "@" + all_locustag_df["Locus_Tag"] 
    return all_locustag_df

def gene_class_def(freq, threshold_99, threshold_15):
    """
    Classify gene based on frequency.

    Parameters:
    freq (int): Frequency of the gene.
    threshold_99 (int): Threshold for core genes.
    threshold_15 (int): Threshold for rare genes.

    Returns:
    str: Gene class ('Core', 'Rare', or 'Accessory').
    """
    if freq >= threshold_99:
        return "Core"
    elif freq < threshold_15:
        return "Rare"
    else:
        return "Accessory"

def reassign_pangene_categories(df_pangene_summary, df_gene_presence_binary):
    n_genome = df_gene_presence_binary.shape[1]
    threshold_99 = round(n_genome * 0.99, 0)
    threshold_15 = round(n_genome * 0.15, 0)
    df_pangene_summary["pangenome_class_2"] = df_pangene_summary["No. isolates"].apply(
        gene_class_def, threshold_99=threshold_99, threshold_15=threshold_15
    )
    return df_pangene_summary

def prepare_qcqa(all_locustag_df, df_pangene_summary, df_gene_presence_binary, df_gene_presence_locustag):
    logging.info(f"all_locustag_df shape: {all_locustag_df.shape}")
    # logging.info("Calculating: Nucleotide lengths")
    # all_locustag_df["Nucleotide_Len"] = all_locustag_df["Nucleotide_Seq"].str.len()
    logging.info("Calculating: Gene_ID counts")
    all_genes_df = all_locustag_df["Gene_ID"].value_counts()
    logging.info("Calculating: Gene_ID counts: Renaming")
    all_genes_df.name = "Number_of_strains"
    logging.info("Calculating: Gene_ID counts: Transforming")
    all_genes_df = all_genes_df.to_frame()

    logging.info("Calculating: Gene lengths")
    gene_lengths = all_locustag_df.groupby("Gene_ID")["Nucleotide_Len"].agg(['sum', 'mean', 'std'])
    gene_lengths.rename(columns={"sum": "Gene_Len_Total", "mean": "Gene_Len_Mean", "std": "Gene_Len_Std"}, inplace=True)
    all_genes_df = all_genes_df.join(gene_lengths)
    all_genes_df["mean_2std"] = all_genes_df["Gene_Len_Mean"] - 2*all_genes_df["Gene_Len_Std"]

    logging.info("Calculating: QC1")
    all_locustag_df["less_than_2std_less_than_mean"] = all_locustag_df["Nucleotide_Len"] < all_locustag_df[["Gene_ID"]].join(all_genes_df[["mean_2std"]], on="Gene_ID")["mean_2std"]

    # Number of strains/genomes is equal to the number of columns in the presence df
    total_num_strains = df_gene_presence_locustag.shape[1]

    logging.info("Calculating: Number strains")
    num_5p_strains = total_num_strains * 0.05
    all_genes_df["5%_of_strains"] = num_5p_strains
    all_genes_df["greater_than_5%"] = (
            all_genes_df["Number_of_strains"] > num_5p_strains
        )
    
    logging.info("Calculating: Pangene categories")
    pangene_summary = reassign_pangene_categories(df_pangene_summary, df_gene_presence_binary)

    logging.info("Calculating: Joining results")
    all_genes_df = all_genes_df.join(pangene_summary["pangenome_class_2"])
    return pangene_summary, all_locustag_df, all_genes_df

def select_genes_and_alleles(all_locustag_df, all_genes_df):
    sel_genes_df = pd.DataFrame(index=all_genes_df.index)
    sel_genes_df["Core"] = all_genes_df["pangenome_class_2"] == "Core"
    sel_genes_df["Pan"] = all_genes_df["greater_than_5%"]

    sel_locustag_df = pd.DataFrame(index=all_locustag_df.index)
    sel_locustag_df["passed"] = ~all_locustag_df["less_than_2std_less_than_mean"]
    return sel_locustag_df, sel_genes_df


def gene_list(sel_genes_df, pan=False):
    return sel_genes_df.index[sel_genes_df["Pan" if pan else "Core"]]

def locustag_list(sel_locustag_df, pan=False):
    return sel_locustag_df.index[sel_locustag_df["passed"]]

def write_gene_list(gene_list, gene_list_path):
    with open(gene_list_path, "w") as f:
        f.writelines(f"{gene}\n" for gene in gene_list)

def load_gene_list(gene_list_path):
    with open(gene_list_path, "r") as f:
        gene_list = [gene.strip() for gene in f.readlines()]
    return gene_list