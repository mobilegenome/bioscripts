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
import numpy as np
import pandas as pd
from Bio import SeqIO
from pathlib import Path

pd.options.mode.chained_assignment = None  # default='warn'


COORDINATE_BASE = 1  # 1 for Bismark coverage file, 0 for BED file
# set to the following to get upstream and downstream information separately ["+1", "-1", "-1+1"]:
DIRECTIONS = ["-1+1"]


# for testing
#
#args = parser.parse_args("-f test/input.fa -b test/input.bismark.cov -c 10 -o test/output".split(" "))


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


def process_genome(bed_chunks: list, fpath_output: [Path | str]) -> None:
    n_chunk = len(bed_chunks)
    print("start processing")

    for chunk_no, chunk in enumerate(bed_chunks):
        print(f"processing chunk no {chunk_no+1} / {n_chunk}", end="\r")
        for _direction in DIRECTIONS:
            chunk_seqrecords = get_sequence(chunk[["seqname", "start", "end"]].values, direction=_direction)
            chunk_sequences = list(map(lambda x: str(x.seq), chunk_seqrecords))
            chunk.loc[:, f"REF_{_direction}"] = chunk_sequences

        chunk.to_csv(fpath_output, mode='a',
                     header=not os.path.exists(fpath_output), sep="\t", index=False,
                     )
    print(f"finished. wrote annotated file to {fpath_output}")


def get_CpGs(df_bed: pd.DataFrame) -> None:
    # split the BED dataframe into chunks of specified size
    chunk_size = args.chunk_size
    bed_chunks = [df_bed[i:i+chunk_size] for i in range(0, len(df_bed), chunk_size)]

    # extract DNA sequences for each BED chunk
    if os.path.exists(args.output_bed):
        print("Error: Outputfile exists!")
        quit()

    process_genome(bed_chunks, args.output_bed)


if __name__ == "__main__":
    # define command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta_file", help="path to the input FASTA file", required=True)
    parser.add_argument("-b", "--bed_file", help="path to the input BED file", required=True)
    parser.add_argument("-o", "--output_bed", help="path to the output BED file", required=True)
    parser.add_argument("-c", "--chunk_size", help="chunk size for reading FASTA file", type=int)

    args = parser.parse_args()

    # read the input FASTA file using BioPython
    genome = SeqIO.to_dict(SeqIO.parse(args.fasta_file, "fasta"))

    # read the input BED file using pandas
    df_bed = pd.read_table(args.bed_file, names=["seqname", "start", "end", "perc_mCpG", "numCs", "numTs"],
                           dtype={"seqname": str,
                                  "start": int,
                                  "end": int,
                                  "perc_mCpG": float,
                                  "numCs": int,
                                  "numTs": int})

    get_CpGs(df_bed)

