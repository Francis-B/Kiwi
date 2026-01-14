"""
Some functions to be used for the digestion
"""

from collections import defaultdict


def read_fasta(fasta_file, header_extractor=None):
    """Convert fasta file to python dictionary"""
    fasta = defaultdict(str)
    with open(fasta_file, "r", encoding="UTF-8") as file:
        header = ""
        for line in file:
            # Get header
            if line.startswith(">"):
                if header_extractor:
                    header = header_extractor(line)
                else:
                    header = line.strip(">\n")

            # Add sequence
            else:
                assert header != "", "Trying to assign sequence to an empty header"

                fasta[header] += line.strip("\n")

    return fasta
