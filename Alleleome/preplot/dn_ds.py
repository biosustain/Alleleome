import pandas as pd
import numpy as np
import json
from pathlib import Path
import logging


def calculate_dn_ds(codon_muts_path, dn_ds_path, dn_ds_json_path):
    logging.info("Starting: preplot: calculate_dn_ds")

    logging.info("Executing: calcualte_dn_ds: read_csv")
    df = pd.read_csv(
        codon_muts_path,
        dtype=str,
        usecols=[
            "Gene",
            "Cons_codon",
            "Query_codon",
            "Codon_position",
            "AA_mutation_effect",
            "GCF_id",
        ],
    )
    logging.info("Executing: calcualte_dn_ds: groupby")
    df = (
        df.groupby(
            [
                "Gene",
                "Cons_codon",
                "Query_codon",
                "Codon_position",
                "AA_mutation_effect",
            ]
        )["GCF_id"]
        .count()
        .reset_index()
    )

    # Remove rows where 'Query_codon' contains 'N'
    logging.info("Executing: calcualte_dn_ds: select")
    df = df[~df["Query_codon"].str.contains("N")]

    logging.info("Executing: calcualte_dn_ds: groupby2")
    df_per_gene = (
        df.groupby(["Gene", "AA_mutation_effect"])["GCF_id"]
        .sum()
        .reset_index(name="Total_per_gene")
    )

    logging.info("Executing: calcualte_dn_ds: pivot")
    df_per_gene = df_per_gene.pivot(
        index=["Gene"], columns="AA_mutation_effect"
    )  # .fillna(0)
    logging.info("Executing: calcualte_dn_ds: map")
    df_per_gene.columns = df_per_gene.columns.to_flat_index().map(
        lambda x: f"{x[1].replace('-', '_')}_count"
    )
    logging.info("Executing: calcualte_dn_ds: rename")
    df_per_gene = df_per_gene.reset_index().rename_axis(None, axis=1)

    logging.info("Executing: calcualte_dn_ds: divide")
    df_per_gene["dN/dS"] = (
        df_per_gene["Non_synonymous_count"] / df_per_gene["Synonymous_count"]
    )
    logging.info("Executing: calcualte_dn_ds: np.inf")
    df_per_gene.replace([np.inf, -np.inf], np.nan, inplace=True)

    # Avoid comparing 'dN/dS' directly, since that involves float comparison
    logging.info("Executing: calcualte_dn_ds: group")
    df_per_gene["dN/dS Group"] = "dNdS_equal_to_one"
    logging.info("Executing: calcualte_dn_ds: group_greater")
    df_per_gene.loc[
        df_per_gene["Non_synonymous_count"] > df_per_gene["Synonymous_count"],
        "dN/dS Group",
    ] = "dNdS_greater_than_one"
    logging.info("Executing: calcualte_dn_ds: group_less")
    df_per_gene.loc[
        df_per_gene["Non_synonymous_count"] < df_per_gene["Synonymous_count"],
        "dN/dS Group",
    ] = "dNdS_less_than_one"

    logging.info("Executing: calcualte_dn_ds: dropna")
    df_per_gene.dropna(axis=0, inplace=True)
    logging.info("Executing: calcualte_dn_ds: to_csv")
    df_per_gene.to_csv(dn_ds_path)

    logging.info("Executing: calcualte_dn_ds: json_dict")
    json_dict = {
        x: df_per_gene.loc[
            df_per_gene["dN/dS Group"] == x,
            ["Synonymous_count", "Non_synonymous_count"],
        ]
        .to_numpy()
        .tolist()
        for x in ["dNdS_greater_than_one", "dNdS_equal_to_one", "dNdS_less_than_one"]
    }
    logging.info("Executing: calcualte_dn_ds: json_dump")
    with open(dn_ds_json_path, "w") as f:
        json.dump(json_dict, f)

    logging.info("Finishing: preplot: calculate_dn_ds")
