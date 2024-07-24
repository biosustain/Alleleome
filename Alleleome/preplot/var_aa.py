import pandas as pd
import numpy as np
from pathlib import Path
import json
import logging

def find_variable_aa(aa_vars_path, variable_aa_path):
    logging.info("Starting: preplot: find_variable_aa")

    df = pd.read_csv(aa_vars_path, na_filter=False, header=0, index_col=False)

    df.drop(index=df.loc[df.AA_mutation_type != "Substitution"].index, inplace=True)
    df.reset_index(inplace=True)

    df["Cons_aa"] = df["AA_cons_seq"].apply(list)
    df["Var_aa"] = df["AA_seq_change"].apply(list)
    s = pd.Series(
        [np.arange(x.AA_start_pos, x.AA_end_pos + 1, 1) for x in df.itertuples()],
        index=df.index,
    )
    df["AA_pos"] = pd.Series(data=s, index=df.index)
    df_f = df.explode(["AA_pos", "Cons_aa", "Var_aa"])

    df_f.to_csv(variable_aa_path)
    logging.info("Finishing: preplot: find_variable_aa")


def find_dominant_var_all(
    variable_aa_path,
    dominant_aa_path,
    dom_var_path,
    gaps_path,
    filt_norm_path,
    dom_var_out_dir,
    gene_list,
):
    logging.info("Starting: preplot: find_dominant_var_all")
    dom_var_out_dir = Path(dom_var_out_dir)

    # This is a file containing all separeted substitution positions along with variant details
    df = pd.read_csv(variable_aa_path, low_memory=False)

    df_var = (
        df[["AA_pos", "Var_aa", "Sequence_type"]]
        .groupby(df["Gene"])
        .value_counts()
        .reset_index(name="AA_freq")
    )

    df_dom = pd.read_csv(
        dominant_aa_path,
        low_memory=False,
    )

    df_dom.rename(
        columns={"AA_cons": "Amino_acid", "AA_freq": "Genome_count"}, inplace=True
    )
    df_var.rename(
        columns={"Var_aa": "Amino_acid", "AA_freq": "Genome_count"}, inplace=True
    )

    df_dom_var = pd.concat([df_var, df_dom]).sort_values("Gene")
    df_dom_var.reset_index(inplace=True)
    df_dom_var.drop(columns=["index"], inplace=True)

    df_dom_var.to_csv(dom_var_path)

    df_gap = df_dom_var.loc[(df_dom_var["Amino_acid"] == "-"), ["Gene", "AA_pos"]]
    df_gap.to_csv(gaps_path)

    df_filt = df_dom_var.loc[(df_dom_var["Amino_acid"] != "-"), :]

    cols = ["Genome_count"]
    g = df_filt.groupby("Gene")[cols]
    min1 = g.transform("min")
    max1 = g.transform("max")

    df_norm = df_filt.join(df_filt[cols].sub(min1).div(max1 - min1).add_suffix("_norm"))

    df_norm["Genome_count_norm"] = df_norm["Genome_count_norm"].apply(
        lambda x: round(x, 2)
    )
    df_norm.to_csv(filt_norm_path)

    # for getting the single gene details
    for gene, group in df_norm.groupby('Gene'):
        if not gene in gene_list:
            continue
        gene_path = dom_var_out_dir / gene
        gene_path.mkdir(parents=True, exist_ok=True)

        out_file = gene_path / f"{gene}_pan_aa_thresh_core_dom_var_pos.csv"
        if out_file.is_file():
            logging.info(f"Skipping {gene}_pan_aa_thresh_core_dom_var_pos.csv, because file already exists.")
            continue
        group.to_csv(out_file)
        # df_norm[df_norm.Gene == gene].to_csv(
        #     out_file
        # )
    logging.info("Finishing: preplot: find_dominant_var_all")

def dom_var_histogram(filt_norm_path, hist_path):
    logging.info("Starting: preplot: dom_var_histogram")
    # Read the data
    df = pd.read_csv(filt_norm_path)

    bins = np.arange(0, 1.01, 0.05)
    dominant_counts, dominant_bins = np.histogram(df.loc[df["Sequence_type"] == "Dominant", "Genome_count_norm"], bins=bins)
    variant_counts, variant_bins = np.histogram(df.loc[df["Sequence_type"] == "Variant", "Genome_count_norm"], bins=bins)

    # Preparing data for Highcharts
    dominant_data_highcharts = [{'x': float(x), 'y': int(y)} for x, y in zip(dominant_bins, dominant_counts)]
    variant_data_highcharts = [{'x': float(x), 'y': int(y)} for x, y in zip(variant_bins, variant_counts)]

    # For plotting repat last y values
    dominant_data_highcharts.append({'x': 1.0, 'y': dominant_data_highcharts[-1]['y']})
    variant_data_highcharts.append({'x': 1.0, 'y': variant_data_highcharts[-1]['y']})

    # Combine into one dictionary
    highcharts_data = {
        "dominant": dominant_data_highcharts,
        "variant": variant_data_highcharts
    }

    with open(hist_path, 'w') as f:
        json.dump(highcharts_data, f)
    logging.info("Finishing: preplot: dom_var_histogram")
