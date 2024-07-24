import argparse
import datetime
import logging
import os

from . import (
    load_and_qcqa,
    write_fasta,
    consensus_sequence,
    sequence_alignment,
    mutations_parallel,
)
from .preplot import (
    dominant_aa,
    var_aa,
    dn_ds,
    msa_freq,
)

# log_directory = "./log"
# os.makedirs(log_directory, exist_ok=True)

# current_time = datetime.datetime.now()
# log_filename = os.path.join(
#     log_directory, f"alleleome_{current_time.strftime('%Y-%m-%d_%H:%M:%S')}.log"
# )

logging.basicConfig(
    # filename=log_filename,
    level=logging.DEBUG,
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def main_prepare(args):
    logging.info("Loading roary data")
    df_pangene_summary, df_gene_presence_binary, df_gene_presence_locustag = (
        load_and_qcqa.load_roary_data(args.summary, args.gp_binary, args.gp_locustag)
    )
    logging.info("Processing GBK files")
    all_locustag_df = load_and_qcqa.parse_genbank_files(
        df_gene_presence_locustag, args.gbk_folder,
    )
    logging.info("Preparing QCQA")
    pangene_summary, all_locustag_df, all_genes_df = load_and_qcqa.prepare_qcqa(
        all_locustag_df,
        df_pangene_summary,
        df_gene_presence_binary,
        df_gene_presence_locustag,
    )
    sel_locustag_df, sel_genes_df = load_and_qcqa.select_genes_and_alleles(
        all_locustag_df, all_genes_df
    )
    logging.info("Saving results")
    load_and_qcqa.save(pangene_summary, args.summary_v2)
    load_and_qcqa.save(all_locustag_df, args.all_locustag)
    load_and_qcqa.save(all_genes_df, args.all_genes)
    load_and_qcqa.save(sel_locustag_df, args.sel_locustag)
    load_and_qcqa.save(sel_genes_df, args.sel_genes)


def load_all_gene_data(args):
    all_genes_df = load_and_qcqa.load(args.all_genes)
    sel_locustag_df = load_and_qcqa.load(args.sel_locustag)
    sel_genes_df = load_and_qcqa.load(args.sel_genes)

    gene_list = load_and_qcqa.gene_list(sel_genes_df, pan=args.pan)
    locustag_list = load_and_qcqa.locustag_list(sel_locustag_df, pan=args.pan)
    return all_genes_df, sel_locustag_df, sel_genes_df, gene_list, locustag_list


def main_fasta(args):
    all_genes_df, sel_locustag_df, sel_genes_df, gene_list, locustag_list = (
        load_all_gene_data(args)
    )
    all_locustag_df = load_and_qcqa.load(args.all_locustag)
    write_fasta.process_selected_genes(
        all_locustag_df, locustag_list, gene_list, args.out_dir, args.gbk_folder
    )
    load_and_qcqa.write_gene_list(gene_list, args.gene_list)


def main_process(args):
    gene_list = load_and_qcqa.load_gene_list(args.gene_list)
    consensus_sequence.build_consensus(gene_list, args.out_dir, p=args.p)
    sequence_alignment.align_sequences(gene_list, args.out_dir, "amino_acid", p=args.p)
    sequence_alignment.align_sequences(gene_list, args.out_dir, "nucleotide", p=args.p)


def main_process_gene(args):
    consensus_sequence.build_single_gene_consensus(args.gene_id, args.out_dir)
    sequence_alignment.align_single_gene(args.gene_id, args.out_dir, "amino_acid")
    sequence_alignment.align_single_gene(args.gene_id, args.out_dir, "nucleotide")


def main_analyze(args):
    gene_list = load_and_qcqa.load_gene_list(args.gene_list)
    mutations_parallel.generate_amino_acid_vars(gene_list, args.out_dir, args.aa_vars, p=args.p)
    mutations_parallel.codon_mut(gene_list, args.out_dir, args.codon_muts, p=args.p)


def main_preplot(args):
    gene_list = load_and_qcqa.load_gene_list(args.gene_list)
    dn_ds.calculate_dn_ds(args.codon_muts, args.dn_ds, args.dn_ds_json)
    dominant_aa.find_dominant_aa(gene_list, args.out_dir, args.dominant_aa, p=args.p)
    var_aa.find_variable_aa(args.aa_vars, args.variable_aa)
    var_aa.find_dominant_var_all(
        args.variable_aa,
        args.dominant_aa,
        args.dom_var,
        args.gaps,
        args.filt_norm,
        args.dom_var_out_dir,
        gene_list,
    )
    var_aa.dom_var_histogram(args.filt_norm, args.hist)
    msa_freq.calculate_msa_freq(gene_list, args.out_dir, args.aa_freq_dir)

def ask_select_mode(args):
    logging.error("Please select a mode, see --help for more info.")


def main():
    logging.info("Application started")
    parser = argparse.ArgumentParser(
        description=(
            "Alleleome - Explore and analyze natural sequence "
            "variations within the Open Reading Frames (ORFs) of "
            "alleles of core genes in a species pan-genome."
        )
    )
    parser.set_defaults(func=ask_select_mode)
    subparsers = parser.add_subparsers()

    # parser.add_argument("mode", type=str, choices=["prepare", "fasta", "process", "process_gene", "analyze"])
    modes = {
        "prepare": main_prepare,
        "fasta": main_fasta,
        "process": main_process,
        "process_gene": main_process_gene,
        "analyze": main_analyze,
        "preplot": main_preplot,
    }
    parsers = {}
    for x, f in modes.items():
        parsers[x] = subparsers.add_parser(x)
        parsers[x].set_defaults(func=f)

    parsers["prepare"].add_argument(
        "--gp_binary",
        type=str,
        required=True,
        help="Path to gene_presence_binary csv file.",
    )
    parsers["prepare"].add_argument(
        "--gp_locustag",
        type=str,
        required=True,
        help="Path to gene_presence_locustag csv file.",
    )
    parsers["prepare"].add_argument(
        "--summary",
        required=True,
        help=("Path to df_pangene_summary.csv file created by roary."),
    )
    parsers["prepare"].add_argument(
        "--summary_v2",
        required=True,
        help=("Path to the updated summary."),
    )
    for x in ["prepare", "fasta"]:
        parsers[x].add_argument(
            "--gbk_folder",
            type=str,
            required=True,
            help="Folder containing GenBank files.",
        )
        parsers[x].add_argument(
            "--all_locustag",
            type=str,
            required=True,
            help="Path to all_locustags csv file.",
        )
        parsers[x].add_argument(
            "--all_genes", type=str, required=True, help="Path to all_genes csv file."
        )
        parsers[x].add_argument(
            "--sel_locustag",
            type=str,
            required=True,
            help="Path to sel_locustags csv file.",
        )
        parsers[x].add_argument(
            "--sel_genes", type=str, required=True, help="Path to sel_genes csv file."
        )
    for x in ["fasta", "process", "process_gene", "analyze", "preplot"]:
        parsers[x].add_argument(
            "--out_dir",
            type=str,
            required=True,
            help="Path to store alignment outputs.",
        )
    parsers["process_gene"].add_argument(
        "--gene_id",
        type=str,
        required=True,
        help="Gene ID when processing a single gene.",
    )
    for x in ["analyze", "preplot"]:
        parsers[x].add_argument(
            "--aa_vars",
            type=str,
            required=True,
            help="Path to pan_amino_acid_vars_df csv file.",
        )
    for x in ["analyze", "preplot"]:
        parsers[x].add_argument(
            "--codon_muts",
            type=str,
            required=True,
            help="Path to pan_gene_syno_non_syno_df csv file.",
        )
    parsers["fasta"].add_argument(
        "--pan",
        action=argparse.BooleanOptionalAction,
        help="Analyze the pangenome, instead of core.",
    )
    for x in ["fasta", "process", "analyze", "preplot"]:
        parsers[x].add_argument(
            "--gene_list",
            type=str,
            required=True,
            help="Path to gene_list.txt file.",
        )
    for x in ["process", "analyze", "preplot"]:
        parsers[x].add_argument(
            "-p",
            type=int,
            required=False,
            default=1,
            help="Number of parallel processes to spawn.",
        )
    parsers["preplot"].add_argument(
        "--dominant_aa",
        type=str,
        required=True,
        help="Path to final_core_consensus_dominant_aa_count_df csv file.",
    )
    parsers["preplot"].add_argument(
        "--variable_aa",
        type=str,
        required=True,
        help="Path to final_core_pan_aa_thresh_vars_all_substitutions_sep_df csv file.",
    )
    parsers["preplot"].add_argument(
        "--dom_var",
        type=str,
        required=True,
        help="Path to final_pan_aa_thresh_core_genes_dominant_variant_genome_count_pos csv file.",
    )
    parsers["preplot"].add_argument(
        "--gaps",
        type=str,
        required=True,
        help="Path to pan_aa_thresh_core_genes_aa_pos_with_gaps csv file.",
    )
    parsers["preplot"].add_argument(
        "--filt_norm",
        type=str,
        required=True,
        help="Path to final_pan_aa_thresh_core_genes_dom_var_genome_count_pos_normalized csv file.",
    )
    parsers["preplot"].add_argument(
        "--dom_var_out_dir",
        type=str,
        required=True,
        help="Path to store per gene variants.",
    )
    parsers["preplot"].add_argument(
        "--dn_ds",
        type=str,
        required=True,
        help="Path to final_dn_ds_count_per_gene csv file.",
    )
    parsers["preplot"].add_argument(
        "--dn_ds_json",
        type=str,
        required=True,
        help="Path to dn_ds.json file.",
    )
    parsers["preplot"].add_argument(
        "--hist",
        type=str,
        required=True,
        help="Path to the histogram data step_line.json.",
    )
    parsers["preplot"].add_argument(
        "--aa_freq_dir",
        type=str,
        required=True,
        help="Path to the directory to store AA_freq.json file per gene.",
    )
    # parser.add_argument(
    #     "--version", action="version", version="%(prog)s " + __version__
    # )

    args = parser.parse_args()
    args.func(args)
