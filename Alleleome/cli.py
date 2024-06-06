import argparse
import datetime
import logging
import os

from . import (
    QCQA_1,
    QCQA_2,
    QCQA_3,
    QCQA_4,
    __version__,
    amino_acid_sequence_alignment,
    build_consensus_sequence,
    codon_mutations,
    generate_amino_acid_variants,
    nucleotide_sequence_alignment,
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


def main():
    logging.info("Application started")
    try:
        parser = argparse.ArgumentParser(
            description=(
                "Alleleome - Explore and analyze natural sequence "
                "variations within the Open Reading Frames (ORFs) of "
                "alleles of core genes in a species pan-genome."
            )
        )
        parser.add_argument("type", type=str, choices=["Core", "Pan"])
        parser.add_argument("--path1", help="Path to pangenome_alignments directory")
        parser.add_argument(
            "--path2",
            help=(
                "Path to alleleome output directory containing "
                "pangene_summary_v2.csv file generated by Roary"
            ),
        )
        parser.add_argument(
            "--table",
            help=(
                "Path to a custom CSV pangene summary table. If not provided, "
                "pangene_summary_v2.csv file (generated by Roary) in the given "
                "Alleleome output directory a will be used."
            ),
            default=None,
        )
        parser.add_argument(
            "--log_to_terminal",
            action="store_true",
            help="Log message will be printed to the terminal instead of a file.",
        )
        parser.add_argument(
            "--version", action="version", version="%(prog)s " + __version__
        )

        args = parser.parse_args()

        if args.log_to_terminal:
            for handler in logging.root.handlers[:]:
                logging.root.removeHandler(handler)
            logging.basicConfig(
                level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
            )
        else:
            pass

        QCQA_1.process_nucleotide_sequences(
            args.path1, args.path2, pangene_summary_csv=args.table, pan_core=args.type
        )
        QCQA_2.analyze_gene_lengths(
            args.path1, args.path2, pangene_summary_csv=args.table, pan_core=args.type
        )
        QCQA_3.process_sequences(args.path1, args.path2)
        QCQA_4.process_genes(
            args.path1, args.path2, pangene_summary_csv=args.table, pan_core=args.type
        )
        build_consensus_sequence.build_consensus(
            args.path1, args.path2, pangene_summary_csv=args.table, pan_core=args.type
        )
        amino_acid_sequence_alignment.amino_acid_seq_align(
            args.path1, args.path2, pangene_summary_csv=args.table, pan_core=args.type
        )
        generate_amino_acid_variants.generate_amino_acid_vars(
            args.path1, args.path2, pangene_summary_csv=args.table, pan_core=args.type
        )
        nucleotide_sequence_alignment.nucleotide_seq_align(
            args.path1, args.path2, pangene_summary_csv=args.table, pan_core=args.type
        )
        codon_mutations.codon_mut(
            args.path1, args.path2, pangene_summary_csv=args.table, pan_core=args.type
        )
        logging.info("Application finished successfully")
    except Exception as e:
        logging.error(f"Application encountered an error: {e}", exc_info=True)
