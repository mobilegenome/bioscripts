from Bio import SeqIO
import argparse
import re
import os
from pathlib import Path


def split_fasta(input_file, pattern, output_dir):
    output_dir = Path(output_dir)
    regex = re.compile(pattern)
    for record in SeqIO.parse(input_file, "fasta"):
        if regex.match(record.id):
            # Use the sequence ID as the filename
            output_filename = Path.joinpath(output_dir, "{}.fasta".format(record.id))
            if not output_dir.exists():
                print(f"Creating {output_dir}")
                output_dir.mkdir(parents=True, exist_ok=True)

            record.id = "{}:0:{}\n".format(record.id, len(record))
            record.description = ""
            SeqIO.write(record,
                        output_filename,
                        "fasta")



if __name__ == "__main__":
    # Set up command line arguments
    parser = argparse.ArgumentParser(description="Split a large FASTA file into individual files.")
    parser.add_argument("-i", "--input_file", metavar="INPUT_FILE", type=str,
                        help="The path to the input FASTA file.")
    parser.add_argument("-p", "--pattern", metavar="PATTERN", type=str,
                        help="The regex pattern to match sequence names.")
    parser.add_argument("-o", "--output_dir", metavar="OUTPUT_DIR", type=str,
                        help="The directory where the output files will be written.")
    args = parser.parse_args()

    # Call the split_fasta function with the input file, pattern, and output directory
    split_fasta(args.input_file, args.pattern, args.output_dir)