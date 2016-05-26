#!/Local/Anaconda-2.0.1-Linux-x86_64/bin/python

import csv
from subprocess import check_output

# exon-starts are 0-based, exon-ends are 1-based
# my_input_file = sys.argv[1]
# my_output_file = sys.argv[2]

my_input_file = "../data/annotations/mm9_prev_version/mm9_ensGene_eric.gpe"
output_file = "../data/extracted_sequences/extracted_utrs_blastdb.fa"
blastb_path = "/storage/md_reut/footprint/mm9/blastdb/mm9"
entries_file = "../data/entries.txt"

MAX_LINES = 100

START = 0
END = 1000000000

NAME_FIELD = 1
CHROM_FIELD = 2
STRAND_FIELD = 3
TX_START_FIELD = 4
TX_END_FIELD = 5
CDS_START_FIELD = 6
CDS_END_FIELD = 7
EXON_COUNT_FIELD = 8
EXON_STARTS_FIELD = 9
EXON_ENDS_FIELD = 10
NAME2_FIELD = 12
RGB = 111 #some random value
SCORE = 500

def run_blastdbcmd(entries):
    # write entries to temporary file
    with open(entries_file, 'w') as file:
        # print "writing to", entries_file, "..."
        file.write("\n".join(entries))

    blastdb_str = "blastdbcmd -db {db} -entry_batch {entries}".format(db = blastb_path, entries = entries_file)
    grep_str = 'grep -v ">"'
    command_str = blastdb_str + "|" + grep_str

    # execute the command
    command_result = check_output(command_str, shell=True)
    return command_result




def proccess_line(line):
     # extract info from line
    strand = line[STRAND_FIELD]
    name = line[NAME_FIELD]
    name2 = line[NAME2_FIELD]
    chrom = line[CHROM_FIELD]
    tx_start = line[TX_START_FIELD]
    tx_end = line[TX_END_FIELD]

    # A note on 0/1 based representation: the basic data file we're processing here is a bigGenePred, which is almost
    # the internal representation of the UCSC genome browser. Start coordinates are 0 based and end coordinates are
    # 1 based. We need coordinates for blastdbcmd, which requires consistent 1 based representation.

    # break the list fields into lists
    exon_starts =  [int(x) +1 for x in line[EXON_STARTS_FIELD].split(",")[:-1]] # -1 for removing the empty member at the end of the list
                                                                                # +1 because the bed file start indices are zero-based
    exon_ends =  [int(x) for x in line[EXON_ENDS_FIELD].split(",")[:-1]] # -1 for removing the empty member at the end of the list

    entries = ""

    if strand == '+':
        cds_start = int(line[CDS_START_FIELD])

        # 3' ----------------------------------------  5'
        # 5' -----ex_start-----ex_end-----cds_start--> 3'

        # 3' ----------------------------------------  5'
        # 5' -----ex_start----cds_start-----ex_end---> 3'

        # list of tuples of coordinates, where the the start of the exon is before the start of the CDS
        exons_to_keep = [(int(exon_start),int(exon_ends[i]))
                         for i, exon_start in enumerate(exon_starts)
                         if int(exon_start) < cds_start]
        entries = ["{chr} {start}-{end} plus".format(chr=chrom, start=start, end=end) for start, end in exons_to_keep]

    elif strand == '-':

        # 3' <---cds_start---ex_start-----ex_end----- 5'
        # 5' ---------------------------------------- 3'

        # 3' <-----ex_start----cds_start-----ex_end-- 5'
        # 5' ---------------------------------------- 3'

        cds_end = int(line[CDS_END_FIELD])

        # list of tuples of coordinates, where the the start of the exon is before the start of the CDS
        exons_to_keep = [(int(exon_starts[i]),int(exon_end))
                         for i, exon_end in enumerate(exon_ends)
                         if int(exon_end) > cds_end]
        entries = ["{chr} {start}-{end} minus".format(chr=chrom, start=start, end=end) for start, end in exons_to_keep]
        entries.reverse() # since this will give the reverse compliment, in order to get a correct concatenation we want to reverse the order of the exons.


    if len(entries) > 1:
        first_kept_start = exons_to_keep[0][0]
        first_kept_end = exons_to_keep[0][1]
        first_kept_size = first_kept_end - first_kept_start
        res = run_blastdbcmd(entries)
        utr_sequence = res.replace('\n', '')
        return  [name + "_" + name2, chrom, first_kept_start, first_kept_size, strand, utr_sequence]
    return []

if __name__ == '__main__':

    print "Proccessing", my_input_file, "..."
    output_lines = []
    count = 0
    with open(my_input_file, 'r') as tsv:
        reader = csv.reader(tsv,  delimiter="\t")
        reader.next()
        for line in reader:
            if count > MAX_LINES:
                break
            count += 1

            if line[CDS_START_FIELD] == line[CDS_END_FIELD]:
                continue #skeep non-coding RNAs

            proccessed_line = proccess_line(line)
            if len(proccessed_line) > 0:
                output_lines.append(proccessed_line)

        output_header = ["name", "chr", "first_exon_start", "first_exon_size", "strand", "5'_UTR"]
    with open(output_file, 'w') as tsv:
            writer = csv.writer(tsv,  delimiter="\t")
            print "writing to", output_file, "..."
            writer.writerow(output_header)
            writer.writerows(output_lines)

