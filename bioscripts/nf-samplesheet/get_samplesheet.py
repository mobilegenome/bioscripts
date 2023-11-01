#!/usr/env python

import pandas as pd
import argparse
import re
import sys

from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("samples",
                    type=str,
                    nargs='+',
                    help="one or more sample name")
parser.add_argument("--single",
                    action="store_true",
                    default=False,
                    help="single-end reads?")
parser.add_argument("-d", "--dir",
                    help="directory with FASTQ files")

args = parser.parse_args()

in_dir = Path(args.dir)

FILE_SUFFIX = ".fastq.gz"
FILE_PATTERN_REGEX = re.compile("^(?P<sample>.+)_(?P<lane>L[1-4])_(?P<read>R[12])[._](?P<rest>.*$)")




files = list(map(lambda path: path.absolute(), in_dir.glob(f"*{FILE_SUFFIX}")))

df = pd.DataFrame(data={"fpath": files})
df.loc[:, "fname"] = df.fpath.apply(lambda fpath: fpath.name)

df[["sample", "lane", "read", "_rest"]] = df.fname.str.extract(FILE_PATTERN_REGEX, expand=True)

df = df[df["sample"].isin(args.samples)]

df = df.pivot(index=["sample", "lane"], columns="read", values="fpath")

df = df.rename(columns={"R1": "fastq_1",
                   "R2": "fastq_2"})

if not args.single:
    df["single"] = "false"

df["genome"] = ""
df.reset_index(inplace=True)

df.drop(columns=["lane"], inplace=True)
df = df[["sample", "single", "fastq_1", "fastq_2", "genome"]]
df.to_csv(sys.stdout, index=None)