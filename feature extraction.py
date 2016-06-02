from cogent.parse.consan import sequence
from collections import Counter
from itertools import combinations
import re

NUCLEOTIDES = "ACUG"
def count_nucleotides(sequence):
    header = [char + "_count" for char in NUCLEOTIDES]
    counter = Counter(sequence)
    counts = [counter[char] + counter[char.lower] for char in NUCLEOTIDES]
    return header, counts

def count_denucleotides(sequence):
    pairs = [p[0]+p[1] for p in combinations(NUCLEOTIDES, 2)]
    header = [p + "_count" for p in pairs]
    counts = [sequence.count(p) + sequence.count(p.lower) for p in pairs]
    return header, counts
import numpy as np


def is_number(string):
    try:
        float(string)
    except ValueError:
        return False
    else:
        return True



POSITION = 0
MOTIF = 1
KMER  = 2
ZSCORE  = 3
PVALUE  = 4
workspace = "/home/jonathan/Documents/data"
all_RBPs = {}

def get_RBP_motifs_from_gene(gene_name):
    motif = "a - motif"
    occurances = "b - occurances"
    positions = "c - positions"
    k_mers = "d - k_mers"
    p_values = "e - p_values"
    average_pvalue = "f - average_pvalue"
    all_fields = ["motif", "occurances", "positions", "k_mers", "p_values", "average_pvalue"]


    f_name = workspace + "/rbp_data.tmp"
    f_name = workspace + "/rbp_motifs/prediction_example.txt"
    # get_motifs_from_rbpMap(f_name)

    prev_protein = ""
    protein = ""
    bNew_motif = False
    dMotif_counts = {}

    with open(f_name) as file:
        for line in file:
            line = line.split()
            if len(line) < 1: continue
            elif not is_number(line[0]):
                if line[0] == "Protein:":

                    # Not the first protein in the file
                    if protein:
                        # Add current gene to the list of genes that have a binding site for this RBP in their UTR
                        try:
                            all_RBPs[protein].append(gene_name)
                        except KeyError: # first gene for this RBP, create the list.
                             all_RBPs[protein] = [gene_name]

                        # calculate stats for previous protein
                        dMotif_counts[protein][average_pvalue] = np.average(np.array([float(x)
                                                                            for x in dMotif_counts[protein][p_values]]))

                    # For all including the first protein in the file
                    protein = line[1]
                    dMotif_counts[protein] = {  motif: "",
                                                occurances: 0,
                                                positions: [],
                                                k_mers: [],
                                                p_values: [] }

                continue # Not a data line, nothing more to do here

            # A data line, populate the dict with the data
            dMotif_counts[protein][motif] = line[MOTIF]
            dMotif_counts[protein][occurances] += 1
            dMotif_counts[protein][positions].append(line[POSITION])
            dMotif_counts[protein][k_mers].append(line[KMER])
            dMotif_counts[protein][p_values].append(line[PVALUE])

        # This is a list of lists (will hold all the lines in the file):
        #   - the outer list is created from dMotif_counts and sorted by the protein name.
        #   - the inner list is sorted by field name, and sets the order of the fields in the output file.

        def key_tuple_by_first(tup):
            return tup[0]

        def sort_by_first_in_tuple(list_of_tuples):
            return sorted(list_of_tuples, key=key_tuple_by_first)

        sorted_proteins_and_counts = sort_by_first_in_tuple(dMotif_counts.items())

        lines = [
                [gene_name, prot] +
                [value for key, value in sort_by_first_in_tuple(dFields.items())] for prot, dFields in sorted_proteins_and_counts
                ]

        header = ["gene_name", "protein_name", "protein_motif", "num_occurences", "AUG_postions", "k-mers", "p-values", "average_pvalue" ]

        print header
        for l in lines:
            print l
        print all_RBPs

def get_RBP_motifs_all_genes():
    f_name = workspace + "/rbp_motifs/prediction_example_multiple_genes.txt"
    collect_lines = False
    lines = []
    extrated_features = []
    gene_name = ""
    with open(f_name) as file:
        for line in file:
            line = line.split()
            if len(line) < 1 :
                continue

            # finish a gene section
            if re.match(r'\*+ ', line[0]):
                collect_lines = False
                extrated_features += get_RBP_motifs_from_gene(gene_name)
                continue

            # add line to be processed by get_RBP_motifs
            if collect_lines:
                lines.append(line)
                continue

            # this line is a gene name
            if re.match(r'ENSMUST.*', line[0]):
                gene_name = re.match(r'ENSMUST.*', line[0]).group(0)
                continue

            # start a new gene section
            if re.match(r'=+ ', line[0]):
                collect_lines = True
                continue


get_RBP_motifs_from_gene("my_gene")
