#!/usr/bin/python

import csv
import pandas as pd
import sys


# exon-starts are 0-based, exon-ends are 1-based
# my_input_file = sys.argv[1]
# my_output_file = sys.argv[2]

my_input_file = "/home/jonathan/Documents/data/mm9_prev_version/mm9_ensGene_eric.gpe"
my_output_file = "/home/jonathan/Documents/data/extracted_utrs.bed"
N = 3

NAME_FIELD = 1
NAME2_FIELD = 12
CHROM_FIELD = 2
STRAND_FIELD = 3
CDS_START_FIELD = 6
CDS_END_FIELD = 7
EXON_COUNT_FIELD = 8
EXON_STARTS_FIELD = 9
EXON_ENDS_FIELD = 10

all_lines = list()
print "Proccessing ", my_input_file, "..."
# df = pd.read_csv(my_input_file, sep="\t")
# for index, row in df.iterrows():
#     # print "\n\n"
#     exonStarts = [int(x) for x in row['exonStarts'].split(",")[:-1]]
#     exonEnds = [int(x) for x in row['exonEnds'].split(",")[:-1]]
#     if row["strand"] == '+':
#         exons_to_keep = [index for index, exonStart in enumerate(exonStarts) if exonStart < int(row["cdsStart"])]
#         # print "proccessing + strand"
#         # print "exonStarts: ", exonStarts
#         # print "exonEnds: ", exonEnds
#         # print "pairwise diffs, End - Start: ", [end - start for start, end in zip(exonStarts,exonEnds)]
#         # print "cdsStart: ", int(row["cdsStart"])
#         # print "exons_to_keep: ", exons_to_keep
#         # print "diffs to edge, CDS_start - Start : ", [int(row["cdsStart"]) - v for  v in exonStarts]
#     else:
#         exons_to_keep = [index for index, exonEnd in enumerate(exonEnds) if exonEnd > int(row["cdsEnd"])] #for "-" strands, translation direction is: 3' ---(CDS) UAC ---- start <---end--- 5'
#         # print "proccessing - strand"
#         # print "exonStarts: ", exonStarts
#         # print "exonEnds: ", exonEnds
#         # print "pairwise diffs, End - Start: ", [end - start for start, end in zip(exonStarts,exonEnds)]
#         # print "cdsEnd: ", int(row["cdsEnd"])
#         # print "exons_to_keep: ", exons_to_keep
#         # print "diffs to edge, End - CDS_end: ", [v - int(row["cdsEnd"]) for  v in exonEnds]
#
#     new_exonStarts = [exonStarts[i] for i in exons_to_keep]
#     new_exonEnds = [exonEnds[i] for i in exons_to_keep]
#
#     new_lines = [[row["chrom"], start, end, row["name"] + "_" + row["name2"], row["strand"]] for start,end in zip(new_exonStarts, new_exonEnds)]
#
#     all_lines += new_lines
#
#
# print "writing to ", my_output_file, "..."
# output_header = ["chrom", "start", "end", "name", "strand"]
# with open(my_output_file, 'w') as tsv:
#     writer = csv.writer(tsv,  delimiter="\t")
#     writer.writerow(output_header)
#     writer.writerows(all_lines)
#
#
#
# exit()

###########################################################################

count = 0
with open(my_input_file, 'r') as tsv:
    reader = csv.reader(tsv,  delimiter="\t")
    reader.next()
    for line in reader:
        strand = line[STRAND_FIELD]

        name = line[NAME_FIELD]
        name2 = line[NAME2_FIELD]
        chrom = line[CHROM_FIELD]
        exon_starts =  [int(x) for x in line[EXON_STARTS_FIELD].split(",")[:-1]] # -1 for removing the empty member at the end of the list
        exon_ends =  [int(x) for x in line[EXON_ENDS_FIELD].split(",")[:-1]] # -1 for removing the empty member at the end of the list

        # print exon_ends
        # print "cds_start:", cds_start
        # print "cds_end:", cds_end
        # print "strand", strand

        if strand == '+':
            cds_start = int(line[CDS_START_FIELD])
            exons_to_keep = [index for index, value in enumerate(exon_starts) if int(value) < cds_start]

            print "cdsStart: ", cds_start
            print "exons_to_keep: ", exons_to_keep
            starts_to_keep = [ex for i, ex in enumerate(exon_starts) if i in exons_to_keep]
            print "diffs to edge, CDS_start - start : ", [cds_start - v for  v in exon_starts]
            print "diffs to edge, chosen exons, CDS_start - start : ", [cds_start - v for  v in starts_to_keep]

        else:
            cds_end = int(line[CDS_END_FIELD])
            exons_to_keep = [index for index, value in enumerate(exon_ends) if int(value) > cds_end]

        new_exon_starts = [exon_starts[i] for i in exons_to_keep]
        new_exon_ends = [exon_ends[i] for i in exons_to_keep]

        new_lines = [[chrom, start, end, name + "_" + name2, strand] for start,end in zip(exon_starts, exon_ends)]
        all_lines += new_lines

        # print "adding lines:"
        # for l in  new_lines:
        #     print l

        count += 1

print "writing to ", my_output_file, "..."
output_header = ["chrom", "start", "end", "name", "strand"]
with open(my_output_file, 'w') as tsv:
    writer = csv.writer(tsv,  delimiter="\t")
    writer.writerow(output_header)
    writer.writerows(all_lines)
