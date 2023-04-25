'''
Merge bookended forward and reverse methylation calls.

This script reads the output of "get_CpG_from_genome.py", separate strand
and merge calls.
'''


import argparse
import re
import numpy as np
import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'

# define command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-b", "--bed_file", help="path to the input BED file")
parser.add_argument("-o", "--output_bed", help="path to the output BED file")

args = parser.parse_args()



def merge_reverse_strand_calls(df: pd.DataFrame,
                               context_column_label: str = "REF_-1+1") -> pd.DataFrame:

    # find and report any NaNs in data
    _df_nan = df[df.isna().any(axis=1)]
    if not _df_nan.empty:
        print("Warning: the following elements of the DataFrame contain NaNs")
        print(_df_nan)
        df = df.dropna()

    if context_column_label not in df.columns:
        print(f"Error. Column {context_column_label} not found in DataFrame.")

    print("prepare dataframs")

    dfs_split = {
        "forward": df[df["REF_-1+1"].str.endswith("CG")],
        "reverse": df[df["REF_-1+1"].str.startswith("CG")]
    }

    # drop unused perc_mCpG column
    dfs_split = {_label: _df.drop("perc_mCpG", axis=1) for _label, _df in dfs_split.items()}

    dfs_split["reverse"]["start-1"] = dfs_split["reverse"]["start"] - 1

    print("joining dataframes")
    df_merged = dfs_split["forward"].merge(dfs_split["reverse"],
                                                how="outer",
                                                left_on=["seqname", "start"],
                                                right_on=["seqname", "start-1"],
                                                suffixes=["_fw", "_rev"])

    # replace unmatched sites, that are only occurring on the reverse strand, ie. where
    # no matching site is found on the forward strand
    row_indexer = np.isnan(df_merged.start_fw) & df_merged.start_rev
    df_merged.loc[row_indexer, ["numCs_fw", "numTs_fw"]] = 0
    df_merged.loc[row_indexer, "start"] = df_merged[row_indexer].start_rev


    # replace unmatched sites on the reverse side, where the corresponding side on the reverse strand
    # is missing
    row_indexer = np.isnan(df_merged.start_rev) & df_merged.start_fw
    df_merged.loc[row_indexer, "numCs_rev"] = 0
    df_merged.loc[row_indexer, "numTs_rev"] = 0
    df_merged.loc[row_indexer, "start"] = df_merged[row_indexer].start_fw


    print("adding read counts")
    df_merged["numCs"] = df_merged.numCs_fw + df_merged.numCs_rev
    df_merged["numTs"] = df_merged.numTs_fw + df_merged.numTs_rev



    print("calculate percentage mCpG")
    df_merged["perc_mCpG"] = df_merged.numCs / (df_merged.numCs + df_merged.numTs) * 100


    # fill left-over NaN "start",
    # i.e. matching forward and reverse sites with the forward start position
    row_indexer = np.isnan(df_merged.start) & df_merged.start_fw & df_merged.start_rev
    df_merged.loc[row_indexer, "start"] = df_merged.start_fw

    # check that no start value is missing:
    if not df_merged.loc[np.isnan(df_merged.start)].empty:
        print("Error. Found missing sites")

    # create a new "end" column that's equal to start
    df_merged["end"] = df_merged["start"]

    # convert count and coordinate columns to integers
    integer_columns = ["start", "end", "numCs", "numTs"]
    df_merged = df_merged.astype({col_label: "int" for col_label in integer_columns})

    print("drop unused columns")
    _column_filter = re.compile(".*_fw$|.*_rev$")
    _column_filter_list = list(filter(_column_filter.search,
                                      df_merged.columns)) + ["start-1"]

    df_merged.drop(_column_filter_list, axis=1, inplace=True)
    df_merged = df_merged[["seqname", "start", "end", "perc_mCpG", "numCs", "numTs"]]

    return df_merged


#args = parser.parse_args("-b get_CpG_from_genome/test/output_tri.bed.gz -o merge_reverse_strand_calls/test/output.bed".split(" "))

df = pd.read_table(args.bed_file)
df2 = merge_reverse_strand_calls(df)
print(f"save dataframe with dimensions {df2.shape} to file")
df2.to_csv(args.output_bed, sep="\t", index=False, header=False, float_format='%.3f')
