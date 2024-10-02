import logging
from pathlib import Path
import shutil

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def process_selected_genes(all_locustag_df, locustag_list, gene_list, out_dir, gbk_folder):
    out_dir = Path(out_dir)
    sel_locustag_df = all_locustag_df.loc[locustag_list]
    
    sel_gene_list = gene_list.to_list()

    for gene in gene_list:
        gene_dir = out_dir / "input" / gene
        gene_dir.mkdir(parents=True, exist_ok=True)
        fna_path = gene_dir / "pan_genes.fna"
        faa_path = gene_dir / "pan_genes.faa"
        if fna_path.is_file() and faa_path.is_file():
            sel_gene_list.remove(gene)
        # process_gene(gene, sel_locustag_df, fna_path, faa_path, tmp_folder)
    if not sel_gene_list:
        return
    genome_ids = sel_locustag_df["Genome_ID"].unique()
    for i, genome_id in enumerate(genome_ids):
        logging.info(f"Writing genome #{i+1}/{len(genome_ids)} to fasta.")
        genbank_file_path = Path(gbk_folder) / f"{genome_id}.gbk"
        for record in SeqIO.parse(genbank_file_path, "genbank"):
            for feature in record.features:
                tag = feature.qualifiers.get("locus_tag")
                if not tag:
                    continue
                locustag = tag[0]
                try:
                    locustag_info = sel_locustag_df.loc[f"{genome_id}@{locustag}", :]
                except:
                    continue
                gene_id = str(locustag_info["Gene_ID"])
                if not gene_id in sel_gene_list:
                    continue
                
                annotation = locustag_info["Prokka_Annotation"]
                if annotation == "":
                    annotation = "Unknown"

                gene_dir = out_dir / "input" / gene_id
                nuc_record = SeqRecord(
                    feature.extract(record.seq),
                    id=locustag, # locus tag
                    description=f"{annotation} | {genome_id}",
                )
                with open(gene_dir / f"pan_genes.fna", "a") as f:
                    SeqIO.write(nuc_record, f, "fasta")
                
                aa_record = SeqRecord(
                    feature.qualifiers["translation"][0] if "translation" in feature.qualifiers.keys() else "",
                    id=locustag, # locus tag
                    description=f"{annotation} | {genome_id}",
                )
                with open(gene_dir / f"pan_genes.faa", "a") as f:
                    SeqIO.write(aa_record, f, "fasta")

def process_all_genes(all_locustag_df, all_genes_df, locustag_list, gene_list, out_dir, gbk_folder):
    out_dir = Path(out_dir)
    sel_locustag_df = all_locustag_df.loc[locustag_list]
    
    sel_gene_list = gene_list.to_list()
    all_gene_list = all_genes_df.index.to_list()

    for gene in all_gene_list:
        gene_dir = out_dir / "input" / gene
        gene_dir.mkdir(parents=True, exist_ok=True)
        fna_path = gene_dir / "pan_genes.fna"
        other_fna_path = gene_dir / "others.fna"
        faa_path = gene_dir / "pan_genes.faa"
        other_faa_path = gene_dir / "others.fna"
        if fna_path.is_file() and faa_path.is_file() and other_fna_path.is_file() and other_faa_path.is_file():
            all_gene_list.remove(gene)
        # process_gene(gene, sel_locustag_df, fna_path, faa_path, tmp_folder)
    # if not sel_gene_list:
    #     return
    genome_ids = all_locustag_df["Genome_ID"].unique()
    for i, genome_id in enumerate(genome_ids):
        logging.info(f"Writing genome #{i+1}/{len(genome_ids)} to fasta.")
        genbank_file_path = Path(gbk_folder) / f"{genome_id}.gbk"
        for record in SeqIO.parse(genbank_file_path, "genbank"):
            for feature in record.features:
                tag = feature.qualifiers.get("locus_tag")
                if not tag:
                    continue
                locustag = tag[0]
                genome_locustag = f"{genome_id}@{locustag}"
                try:
                    locustag_info = all_locustag_df.loc[genome_locustag, :]
                except:
                    continue
                gene_id = str(locustag_info["Gene_ID"])
                if not gene_id in sel_gene_list:
                    continue
                
                annotation = locustag_info["Prokka_Annotation"]
                if annotation == "":
                    annotation = "Unknown"

                gene_dir = out_dir / "input" / gene_id
                nuc_record = SeqRecord(
                    feature.extract(record.seq),
                    id=locustag, # locus tag
                    description=f"{annotation} | {genome_id}",
                )
                if genome_locustag in locustag_list:
                    base_path = gene_dir / "pan_genes"
                else:
                    base_path = gene_dir / "others"
                with open(base_path.with_suffix(".fna"), "a") as f:
                    SeqIO.write(nuc_record, f, "fasta")
                
                aa_record = SeqRecord(
                    feature.qualifiers["translation"][0] if "translation" in feature.qualifiers.keys() else "",
                    id=locustag, # locus tag
                    description=f"{annotation} | {genome_id}",
                )
                with open(base_path.with_suffix(".faa"), "a") as f:
                    SeqIO.write(aa_record, f, "fasta")