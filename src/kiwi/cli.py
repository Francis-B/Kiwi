"""
Function used to parse argument when digestool is run with CLI
"""

import argparse

from kiwi.digestion import Experiment


def parse_args():
    """Parse command line arguments to run a cleave protein of the given
    fasta file"""

    parser = argparse.ArgumentParser(
        description="Digest proteins of a given fasta file"
    )
    parser.add_argument("fasta_file", metavar="</path/to/file.fasta>")
    parser.add_argument(
        "-l",
        "--min_length",
        metavar="",
        default=7,
        help="minimal length of peptide sequences (default: 7)",
    )
    parser.add_argument(
        "-L",
        "--max_length",
        metavar="",
        default=None,
        help="maximal length of peptide sequences (default: None)",
    )
    parser.add_argument(
        "-M",
        "--miscleavages",
        metavar="",
        default=2,
        help="maximum of miscleavages allowed (default: 2)",
    )
    parser.add_argument(
        "-m",
        "--mass",
        metavar="",
        default=4600,
        help="maximal molecular mass of peptide sequences (default: 4600 dalton)",
    )
    parser.add_argument(
        "-c",
        "--clip_nterm_m",
        action="store_true",
        help="Do not consider sequences without it N-term methionine (default: False)",
    )
    parser.add_argument(
        "-u ",
        "--unique",
        action="store_true",
        help="add flag to inform if the peptides are unique or not",
    )
    parser.add_argument(
        "-e",
        "--enzyme",
        metavar="",
        default="trypsin",
        help="enzyme used to perform digestion (default: trypsin)",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="",
        help="output directory (default: /path/to/file.fasta_digestedPeptides.csv)",
    )

    args = parser.parse_args()
    return args


def main():
    # Get command line args
    args = parse_args()

    # Initalize Experiment
    experiment = Experiment(args.fasta_file)

    # Set params
    experiment.set_min_length(args.min_length)

    if args.max_length is not None:
        experiment.set_max_length(args.max_length)

    experiment.set_max_miscleavages(args.miscleavages)
    experiment.set_max_mass(args.mass)
    experiment.set_clip_nterm_m(args.clip_nterm_m)
    experiment.set_enzyme(args.enzyme)

    # Set output dir
    if args.output is not None:
        experiment.set_outdir(args.output)

    print(
        "Performing digestion with:\n\n"
        + f"min_length: {experiment.params['min_length']}\n"
        + f"max_length: {experiment.params['max_length']}\n"
        + f"enzyme: {experiment.params['enzyme']}\n"
        + f"max_mass: {experiment.params['max_mass']}\n"
        + f"max_miscleavages: {experiment.params['max_miscleavages']}\n"
    )

    # Cleave proteins
    experiment.cleave_proteins()

    # Check if peptides are unique
    if args.unique is True:
        experiment.check_peptide_uniqueness()
    experiment.write()
    print(f"List of peptide sequences saved as {experiment.outdir}")


if __name__ == "__main__":
    main()
