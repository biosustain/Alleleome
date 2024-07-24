# Alleleome generation step III- Parsing the results and generating the amino acid mutations
import logging
from itertools import groupby
from operator import itemgetter
from pathlib import Path
import gzip
import pandas as pd
from Bio.Blast import NCBIXML


# String parsing helper functions
# fast find all chars in string: https://stackoverflow.com/questions/52452911/finding-all-positions-of-a-character-in-a-string
def _find_all_idx(string, character):
    idx = string.find(character)
    while idx != -1:
        yield idx
        idx = string.find(character, idx + 1)


def _get_char_idxs_in_str(s, c):
    char_idxs = list(_find_all_idx(s, c))
    char_idxs = [x for x in char_idxs]
    return char_idxs


CHAR_OF_INTEREST = "-"
assert _get_char_idxs_in_str("AR-D-G", CHAR_OF_INTEREST) == [2, 4]

# BLAST results glossary: https://www.ncbi.nlm.nih.gov/books/NBK62051/
# Description of special characters on the match string: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp
# The databases alignments are displayed as pairs of matches between query and subject sequence.
# A middle line between the query and subject sequence displays the status of a letter.
# For protein alignments (e.g, BLASTP/BLASTX/TBLASTN), identities present the letter,
# conservative substitutions present a "+", and nothing otherwise. This is the default view.
# https://www.ncbi.nlm.nih.gov/books/NBK62051/: "If the aligned residues have similar physico-chemical properties or have a
# positive score in the governing scoring matrix the substitution is said to be conservative."
# '-' in subject string indicate an insertion
# BLAST subject and query results string are always the same length.
# Won't find differences in the beginning of genes since doesn't compare query vs subject start and end match positions.
# This logic was removed since alleles can have enormous differences on either of their sequences.
# Best to keep the returned data seperate for finding
# consecutive mutations according to AA positions.

BLAST_INDEL_CHAR = "-"


def _get_ins_d(subject_str, query_str, subject_match_start_pos):
    ins_d = dict()
    for idx in _get_char_idxs_in_str(subject_str, BLAST_INDEL_CHAR):
        pos = idx + subject_match_start_pos
        ins_d[pos] = query_str[idx]
    return ins_d


# Won't find differences in the beginning of genes since doesn't compare query vs subject start and end match positions.
def _get_del_d(subject_str, query_str, subject_match_start_pos):
    del_d = dict()
    for idx in _get_char_idxs_in_str(query_str, BLAST_INDEL_CHAR):
        pos = idx + subject_match_start_pos
        del_d[pos] = subject_str[idx]
    return del_d


sbjct = "AR-D-G"  # insertion @ 3 and 5
qry = "-RBDHG"  # deletion @ 1
assert _get_ins_d(sbjct, qry, 1) == {3: "B", 5: "H"}
assert _get_ins_d(sbjct, qry, 2) == {
    4: "B",
    6: "H",
}  # test effects of subject_match_start_pos

assert _get_del_d(sbjct, qry, 1) == {1: "A"}
assert _get_del_d(sbjct, qry, 2) == {2: "A"}  # test effects of subject_match_start_pos


ALL_AMINO_ACIDS = {
    "A",
    "R",
    "N",
    "D",
    "B",
    "C",
    "E",
    "Q",
    "Z",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V",
}


def _get_sub_aa_d(
    subject_str, match_str, query_str, indel_aa_pos_l, subject_match_start_pos
):
    seq_chng_d = dict()

    for c_i in range(0, len(match_str)):
        if (
            match_str[c_i] not in ALL_AMINO_ACIDS
        ):  # finding BLAST characters indicating mismatch or gap
            pos = (
                subject_match_start_pos + c_i
            )  # subject_match_start_pos+c_i since subject_match_start_pos represents the first charater and the index c_i == 0
            if pos not in indel_aa_pos_l:
                seq_chng_d[pos] = {"s": subject_str[c_i], "q": query_str[c_i]}

    return seq_chng_d


assert _get_sub_aa_d("ARPGCQ", "+ PG Q", "MTPGBQ", [], 1) == {
    1: {"s": "A", "q": "M"},
    2: {"s": "R", "q": "T"},
    5: {"s": "C", "q": "B"},
}
assert _get_sub_aa_d("ARPGCQ", "+ PG Q", "MTPGBQ", [5], 1) == {
    1: {"s": "A", "q": "M"},
    2: {"s": "R", "q": "T"},
}  # testing use of indel_aa_pos_l
assert _get_sub_aa_d("ARPGCQ", "+ PG Q", "MTPGBQ", [], 2) == {
    2: {"s": "A", "q": "M"},
    3: {"s": "R", "q": "T"},
    6: {"s": "C", "q": "B"},
}  # test effects of subject_match_start_pos


def _get_consecutive_pos_l(pos_l):
    cnsc_pos_l = []
    for k, g in groupby(enumerate(pos_l), lambda ix: ix[0] - ix[1]):
        cnsc_pos_l.append(list(map(itemgetter(1), g)))
    return cnsc_pos_l


mut_pos_l = [1, 4, 5, 7, 8, 9, 22]
assert _get_consecutive_pos_l(mut_pos_l) == [[1], [4, 5], [7, 8, 9], [22]]

def generate_amino_acid_vars(args):
    gene, out_dir = args
    return generate_amino_acid_vars_(gene, out_dir)

def generate_amino_acid_vars_(
    gene, out_dir
):
    out_dir = Path(out_dir)
    blast_output_dir = (
        out_dir / "output" / gene
    )
    blast_output_file_path = blast_output_dir / (
        "amino_acid_blast_out_" + gene + ".xml.gz"
    )
    all_mutations = []

    with gzip.open(blast_output_file_path, "rt") as f:
        for record in NCBIXML.parse(f):
            if len(record.alignments) > 0:
                # Description of available members: https://biopython.org/docs/1.75/api/Bio.Blast.Record.html
                subject_match_start_pos = record.alignments[0].hsps[0].sbjct_start
                blast_subject_str = record.alignments[0].hsps[0].sbjct
                blast_match_str = record.alignments[0].hsps[0].match
                blast_query_str = record.alignments[0].hsps[0].query
                length_align = record.alignments[0].hsps[0].align_length
                align_score = record.alignments[0].hsps[0].score
                e_value = record.alignments[0].hsps[0].expect
                num_indentities = record.alignments[0].hsps[0].identities
                query_info = record.query
                query_info = query_info.split("|", 1)
                query_description = query_info[0]
                query_description = query_description.split(" ", 1)
                query_id = query_description[0]
                query_desc = query_description[1]
                GCF_id = query_info[1]
                # Won't find differences in the beginning of genes since doesn't compare query vs subject start and end match positions.
                ins_d = _get_ins_d(
                    blast_subject_str, blast_query_str, subject_match_start_pos
                )
                for cnsc_pos_l in _get_consecutive_pos_l(ins_d.keys()):
                    seq_chng_str = "".join([ins_d[p] for p in cnsc_pos_l])
                    mut_data = {
                        "Gene": gene,
                        "Mutation_source": "Pangenome variant",
                        "AA_start_pos": min(cnsc_pos_l),
                        "AA_end_pos": max(cnsc_pos_l),
                        "Mutation_size": len(cnsc_pos_l),
                        "AA_mutation_type": "Insertion",
                        "AA_cons_seq": "",
                        "AA_seq_change": seq_chng_str,
                        "Length_alignment": length_align,
                        "Score": align_score,
                        "E_value": e_value,
                        "Identity": num_indentities,
                        "Query_locus_tag": query_id,
                        "Query_description": query_desc,
                        "GCF_id": GCF_id,
                        "Sequence_type": "Variant",
                    }  # Can't use columns due to diff parsing of 'AA range' input
                    all_mutations.append(mut_data)
                # Won't find differences in the beginning of genes since doesn't compare query vs subject start and end match positions.
                del_d = _get_del_d(
                    blast_subject_str, blast_query_str, subject_match_start_pos
                )
                for cnsc_pos_l in _get_consecutive_pos_l(del_d.keys()):
                    seq_chng_str = "".join([del_d[p] for p in cnsc_pos_l])
                    mut_data = {
                        "Gene": gene,
                        "Mutation_source": "Pangenome variant",
                        "AA_start_pos": min(cnsc_pos_l),
                        "AA_end_pos": max(cnsc_pos_l),
                        "Mutation_size": len(cnsc_pos_l),
                        "AA_mutation_type": "Deletion",
                        "AA_cons_seq": seq_chng_str,
                        "AA_seq_change": "",
                        "Length_alignment": length_align,
                        "Score": align_score,
                        "E_value": e_value,
                        "Identity": num_indentities,
                        "Query_locus_tag": query_id,
                        "Query_description": query_desc,
                        "GCF_id": GCF_id,
                        "Sequence_type": "Variant",
                    }  # Can't use columns due to diff parsing of 'AA range' input
                    all_mutations.append(mut_data)

                sub_d = _get_sub_aa_d(
                    blast_subject_str,
                    blast_match_str,
                    blast_query_str,
                    list(ins_d.keys()) + list(del_d.keys()),
                    subject_match_start_pos,
                )
                for cnsc_pos_l in _get_consecutive_pos_l(sub_d.keys()):
                    seq_ref_str = "".join([sub_d[p]["s"] for p in cnsc_pos_l])
                    seq_chng_str = "".join([sub_d[p]["q"] for p in cnsc_pos_l])
                    mut_data = {
                        "Gene": gene,
                        "Mutation_source": "Pangenome variant",
                        "AA_start_pos": min(cnsc_pos_l),
                        "AA_end_pos": max(cnsc_pos_l),
                        "Mutation_size": len(cnsc_pos_l),
                        "AA_mutation_type": "Substitution",
                        "AA_cons_seq": seq_ref_str,
                        "AA_seq_change": seq_chng_str,
                        "Length_alignment": length_align,
                        "Score": align_score,
                        "E_value": e_value,
                        "Identity": num_indentities,
                        "Query_locus_tag": query_id,
                        "Query_description": query_desc,
                        "GCF_id": GCF_id,
                        "Sequence_type": "Variant",
                    }  # Can't use columns due to diff parsing of 'AA range' input
                    all_mutations.append(mut_data)
    return all_mutations