# from cogent.parse.consan import sequence
from collections import Counter
from itertools import combinations
import re
import numpy as np
import csv
import decimal
from blastdb_extract_sequences import fasta_output_file

# example to show how to write named matrices to file with pandas
# ----
# import numpy as np
# import pandas as pd
#
# A = np.random.randint(0, 10, size=36).reshape(6, 6)
# names = [_ for _ in 'abcdef']
# df = pd.DataFrame(A, index=names, columns=names)
# df.to_csv('df.csv', index=True, header=True, sep=' ')


# Constants
POSITION = 0
MOTIF = 1
KMER  = 2
ZSCORE  = 3
PVALUE  = 4
NUCLEOTIDES = "ACUG"

# Globals
all_RBPs = {}

#Paths

server = False
workspace = "/home/jonathan/Documents/data/"

if server:
    workspace = "/srv01/technion/jonathans/data/"

rbp_data_folder = workspace + "rbp_motif_analysis/"
rbp_motifs_per_gene_path = rbp_data_folder + "motifs_by_gene_name.tsv"
genes_per_protein_path = rbp_data_folder + "genes_by_protein_name.tsv"

 # Utility functions
def is_number(string):
    try:
        float(string)
    except ValueError:
        return False
    else:
        return True

def key_tuple_by_first(tup):
        return tup[0]

def sort_by_first_in_tuple(list_of_tuples):
    return sorted(list_of_tuples, key=key_tuple_by_first)

# Feature extraction functions
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

def get_RBP_motifs_from_gene(gene_name, input_lines):
    motif = "a - motif"
    occurances = "b - occurances"
    positions = "c - positions"
    k_mers = "d - k_mers"
    average_pvalue = "e - average_pvalue"
    p_values = "f - p_values"
    # all_fields = ["motif", "occurances", "positions", "k_mers", "p_values", "average_pvalue"]
    all_fields = ["motif", "p_values", "average_pvalue"]

    protein = ""
    dMotif_counts = {}

    for line in input_lines:
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
                    decimal.getcontext().prec = 2
                    avg = np.average(np.array([float(x) for x in dMotif_counts[protein][p_values]]))
                    avg = round(decimal.Decimal(avg),4 )
                    dMotif_counts[protein][average_pvalue] = avg

                # For all including the first protein in the file
                protein = line[1]
                dMotif_counts[protein] = {  motif: "",
                                            # occurances: 0,
                                            # positions: [],
                                            # k_mers: [],
                                            p_values: [] }

            continue # Not a data line, nothing more to do here

        # A data line, populate the dict with the data
        dMotif_counts[protein][motif] = line[MOTIF]
        # dMotif_counts[protein][occurances] += 1
        # dMotif_counts[protein][positions].append(line[POSITION])
        # dMotif_counts[protein][k_mers].append(line[KMER])
        dMotif_counts[protein][p_values].append(line[PVALUE])

    # Stringify the lists of values
    for protein in dMotif_counts.keys():
        # dMotif_counts[protein][positions] = ",".join(dMotif_counts[protein][positions])
        # dMotif_counts[protein][k_mers] = ",".join(dMotif_counts[protein][k_mers])
        dMotif_counts[protein][p_values] = ",".join(dMotif_counts[protein][p_values])

    # This is a list of lists (will hold all the lines in the file):
    #   - the outer list is created from dMotif_counts and sorted by the protein name.
    #   - the inner list is sorted by field name, and sets the order of the fields in the output file.


    sorted_proteins_and_counts = sort_by_first_in_tuple(dMotif_counts.items())

    lines = [
            [gene_name, prot] +
            [value for key, value in sort_by_first_in_tuple(dFields.items())[:-1]] for prot, dFields in sorted_proteins_and_counts
            ]

    # header = ["gene_name", "protein_name", "protein_motif", "num_occurences", "AUG_postions", "k-mers", "p-values", "average_pvalue" ]
    header = ["gene_name", "protein_name", "protein_motif","average_pvalue" ]

    return header, lines

def get_RBP_motifs_all_genes(rbp_output_file = rbp_data_folder + "All_Predictions.txt"):
    collect_lines = False
    lines = []
    extrated_features = []
    gene_name = ""
    header = ""
    with open(rbp_output_file) as file:
        for line in file:
            line = line.split()
            if len(line) < 1 :
                continue

            # finish a gene section
            if re.match(r'\*+', line[0]) and collect_lines:
                collect_lines = False
                header, values = get_RBP_motifs_from_gene(gene_name, lines)
                extrated_features += values
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
            if re.match(r'=+', line[0]):
                collect_lines = True
                continue

    with open(rbp_motifs_per_gene_path, 'w+') as file:
        writer = csv.writer(file,  delimiter="\t")
        writer.writerow(header)
        writer.writerows(extrated_features)

    header = ['protein_name', 'genes']
    with open(genes_per_protein_path, 'w+') as file:
        writer = csv.writer(file,  delimiter="\t")
        writer.writerow(header)
        all_proteins = [[protein] + genes for protein, genes in sort_by_first_in_tuple(all_RBPs.items())]
        writer.writerows(all_proteins)

        print  all_proteins

def update_kmer_dict(kmers_genes_dict, all_kmers_set, gene_name, sequence):
    K = 6
    kmers = [sequence[start:start+K] for start in range(len(sequence))]
    kmers_genes_dict[gene_name] = kmers
    [all_kmers_set.add(kmer) for kmer in kmers]

def get_kmers_all_genes(all_genes_fasta = fasta_output_file):
    collect_lines = []
    dKmers_per_gene = {}
    sAll_kmers = set()
    with open(all_genes_fasta) as file:
        for line in file:
            line = line.split()
            if len(line) < 1 :
                continue #skip empty lines

            # this line is a gene name
            if re.match(r'>ENSMUST.*', line[0]):
                gene_name = re.match(r'>(ENSMUST.*)', line[0]).group(1)
                continue

            else: # this is a sequence line
                update_kmer_dict(dKmers_per_gene, sAll_kmers, gene_name, line[0])
    for k, v in dKmers_per_gene.iteritems():
        print k, v
    print sAll_kmers

if __name__ == '__main__':

    # get_RBP_motifs_all_genes()
    get_kmers_all_genes()