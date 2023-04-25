'''
Annotate sequence context for Bismark coverage files.

Bismark's coverage file contain information for all sequence contexts.

This script annotate the .cov file with the sequence context plus/minus one base around
the called methylated site i.e. pos-1 and pos+1.

This makes it possible to filter and merge adjacent methylation calls.

Merging methylation is important to obtain full coverage for the backward strand sites.
If a CpG methylation call is made, the methylation status of G on the backward strand
is listed separately in the coverage file. By identifiying these cases it is possible to merge the
coverage information for these sites.
'''

import argparse
import os
import re
import numpy as np
import pandas as pd
from Bio import SeqIO

pd.options.mode.chained_assignment = None  # default='warn'


COORDINATE_BASE = 1  # 1 for Bismark coverage file, 0 for BED file
# set to the following to get upstream and downstream information separately ["+1", "-1", "-1+1"]:
DIRECTIONS = ["-1+1"]

# define command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta_file", help="path to the input FASTA file")
parser.add_argument("-b", "--bed_file", help="path to the input BED file")
parser.add_argument("-o", "--output_bed", help="path to the output BED file")
parser.add_argument("-c", "--chunk_size", help="chunk size for reading FASTA file", type=int)

args = parser.parse_args()

# for testing
#
#args = parser.parse_args("-f test/input.fa -b test/input.bed -c 10 -o test/output".split(" "))

# read the input FASTA file using BioPython
genome = SeqIO.to_dict(SeqIO.parse(args.fasta_file, "fasta"))

# read the input BED file using pandas
bed = pd.read_table(args.bed_file, names=["seqname", "start", "end", "perc_mCpG", "numCs", "numTs"])

# define a function to extract DNA sequences
def get_sequence(bed_chunk: np.array, direction: str = "+1") -> list:
    seqnames = bed_chunk[:, 0]

    if direction == "+1":
        starts = map(lambda x: x - COORDINATE_BASE, bed_chunk[:, 1])
        ends = map(lambda x: x - COORDINATE_BASE + 2, bed_chunk[:, 2])
    elif direction == "-1":
        starts = map(lambda x: x - COORDINATE_BASE - 1, bed_chunk[:, 1])
        ends = map(lambda x: x, bed_chunk[:, 2])
    elif direction == "-1+1":
        starts = map(lambda x: x - COORDINATE_BASE - 1, bed_chunk[:, 1])
        ends = map(lambda x: x - COORDINATE_BASE + 2, bed_chunk[:, 2])
    else:
        print(f"Unknown direction {direction}. Quitting.")
        quit()

    sequence = map(lambda _seq, _start, _end: genome[_seq][_start:_end], seqnames, starts, ends)
    return list(sequence)


def process_genome(bed_chunks: list) -> None:
    n_chunk = len(bed_chunks)
    print("start processing")

    for chunk_no, chunk in enumerate(bed_chunks):
        print(f"processing chunk no {chunk_no+1} / {n_chunk}", end="\r")
        for _direction in DIRECTIONS:
            chunk_seqrecords = get_sequence(chunk[["seqname", "start", "end"]].values, direction=_direction)
            chunk_sequences = list(map(lambda x: str(x.seq), chunk_seqrecords))
            chunk.loc[:, f"REF_{_direction}"] = chunk_sequences

        chunk.to_csv(args.output_bed, mode='a',
                     header=not os.path.exists(args.output_bed), sep="\t", index=False,
                     )
    print(f"finished. wrote annotated file to {args.output_bed}")


# split the BED dataframe into chunks of specified size
chunk_size = args.chunk_size
bed_chunks = [bed[i:i+chunk_size] for i in range(0, len(bed), chunk_size)]

# extract DNA sequences for each BED chunk
if os.path.exists(args.output_bed):
    print("Error: Outputfile exists!")
    quit()

process_genome(bed_chunks)
# done


bed = pd.read_table("get_CpG_from_genome/test/output_tri.bed.gz")
df[df.isna().any(axis=1)]
bed.dropna(inplace=True)
# define a function to extract DNA sequences

def merge_reverse_strand_calls(df, context_column_label="REF_-1+1"):

    if context_column_label not in df.columns:
        print(f"Error. Column {context_column_label} not found in DataFrame.")

    dfs_split = {
        "forward": df[df["REF_-1+1"].str.endswith("CG")],
        "reverse": df[df["REF_-1+1"].str.startswith("CG")]
    }


    map(lambda df: df.drop("perc_mCpG", axis=1, inplace=True), dfs_split.values]
# filter and subset for testing
df_cpg_fw, df_cpg_rev = [df[df["seqname"] == "chr25"] for df in [df_cpg_fw, df_cpg_rev]]
df_cpg_fw, df_cpg_rev = [df.head() for df in [df_cpg_fw, df_cpg_rev]]

df_cpg_rev["start-1"] = df_cpg_rev["start"] -1

df_cpg_merged = df_cpg_fw.merge(df_cpg_rev,
                                how="left",
                                left_on=["seqname", "start"],
                                right_on=["seqname", "start-1"],
                                suffixes=["_fw", "_rev"])

df_cpg_merged["numCs"] = df_cpg_merged["numCs_fw"] +df_cpg_merged["numCs_rev"]
df_cpg_merged["numTs"] = df_cpg_merged["numTs_fw"] +df_cpg_merged["numTs_rev"]

df_cpg_merged["perc_mCpG"] = \
    df_cpg_merged["numCs"] / (df_cpg_merged["numCs"] + df_cpg_merged["numTs"]) * 100

column_filter = re.compile(".*_fw$|.*_rev$")
column_filter_list = list(filter(column_filter.search, df_cpg_merged.columns))
df_cpg_merged.drop(column_filter_list, axis=1, inplace = True)

del df_cpg_merged

for idx, row in bed.iterrows():
    seqname = row["seqname"]
    pos = row["start"]
    context = row["REF_-1+1"]

    # is a forward CpG site? then nothing to do
    row_cpg = df_cpg[(df_cpg["seqname"] == seqname) & (df_cpg["start"] == pos)]

    # is a reverse CpG site
    if context.startswith("CG"):
        # and previous site is
        row_cpg_previous = df_cpg[(df_cpg["seqname"] == seqname) & (df_cpg["start"] == pos-1)]
        row_loc = (df_cpg["seqname"] == seqname) & (df_cpg["start"] == pos - 1)
        if not row_cpg_previous.empty:
            df_cpg.loc[row_loc, ["numCs", "numTs"]] += row["numCs"], row["numTs"]
            df_cpg.loc[row_loc, "perc_mCpG"] = \
                df_cpg.loc[row_loc, "numCs"] / (df_cpg.loc[row_loc, "numCs"] + df_cpg.loc[row_loc, "numTs"]) * 100
        else:
            df_cpg[row_loc] = row


df_cpg.to_csv(args.output_bed.replace(".bed.", ".CpG_forward.bed."), sep="\t", header=False, index=False, float_format="%.5f")
