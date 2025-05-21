#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
"""
This script define the class Digestion which perform the in silico digestion and
check the uniqueness of peptide sequences.

It contains also the lines for command line usage.
"""

# ------------------------------------------------------------------------------
from textwrap import dedent

from enzyme import Enzyme
from peptide import Peptide
import reader
import cli
# ------------------------------------------------------------------------------


class Digestion:
    """Instantiate this class to perform in silico digestion and find unique
    peptide"""

    # Default parameter. Can be change with self.set<param>(value)
    _enzyme = "trypsin"  # See enzyme.py for more enzymes choice
    _min_length = 7
    _max_length = None
    _max_miscleavages = 1
    _max_mass = 4600  # in dalton

    def __init__(self, fasta):
        self.fasta = fasta
        self.params = {
            "enzyme": self.__class__._enzyme,
            "min_length": self.__class__._min_length,
            "max_length": self.__class__._max_length,
            "max_miscleavages": self.__class__._max_miscleavages,
            "max_mass": self.__class__._max_mass,
            "clip_nterm_m": False,
        }
        self.proteins = reader.fasta_to_dict(self.fasta)
        self.peptides = []
        self.outdir = f"{fasta}_digestedPeptides.csv"
        self.is_unique = {}

    def __str__(self):
        return dedent(f"""\
                        Parameters:

                        Enzyme: {self.params["enzyme"]}
                        Minimum peptide sequence length: {self.params["min_length"]}
                        Maximum peptide sequence length: {self.params["max_length"]}
                        Maximum of miscleavages: {self.params["max_miscleavages"]}
                        Maximum peptide's mass: {self.params["max_mass"]}
                        """)

    def set_enzyme(self, enzyme_name):
        self.params["enzyme"] = enzyme_name.lower()

    def set_min_length(self, minimum):
        self.params["min_length"] = int(minimum)

    def set_max_length(self, maximum):
        self.params["max_length"] = int(maximum)

    def set_max_miscleavages(self, miscleavages):
        self.params["max_miscleavages"] = int(miscleavages)

    def set_max_mass(self, _mass):
        self.params["max_mass"] = float(_mass)

    def set_clip_nterm_m(self, bool_):
        self.params["clip_nterm_m"] = bool_

    def set_outdir(self, dir_):
        self.outdir = dir_

    def _is_cleaved(self):
        """
        Check if proteins were digested.

        Returns
        -------
        bool
            Wheter or not proteins were digested.
        """
        return len(self.peptides) != 0

    def _is_valid_peptide(self, peptide):
        """
        Determine if peptide respect the given criteria.

        Parameters
        ----------
        peptide : Peptide
            Peptide to check.
        """
        # Check length
        if self.params["min_length"] is not None and self.params["min_length"] > len(
            peptide
        ):
            return False
        # Check length
        if self.params["max_length"] is not None and self.params["max_length"] < len(
            peptide
        ):
            return False
        # Check mass
        if (
            self.params["max_mass"] is not None
            and self.params["max_mass"] < peptide.get_mass()
        ):
            return False

        return True

    def _get_valid_peptides(self, peptides):
        """
        Get the peptides that respect the given criteria.

        Parameters
        ----------
        peptides : list of Peptides
            All peptides to check.

        Returns
        -------
        list of Peptide
            All peptides that respect to previously-defined of defaut criteria.
        """
        return [peptide for peptide in peptides if self._is_valid_peptide(peptide)]

    def _format_line(self, peptide, include_mass):
        """
        From the given peptide, get the information needed for the file.

        Parameters
        ----------
        peptide : Peptide
            Peptide object from which the information will be extracted.
        include_mass : bool
            Wheter or not to include mass in the file.

        Returns
        -------
        str
            Line to be written in the tsv file.
        """
        parent_proteins = ",".join(peptide.parent_proteins)
        if include_mass:
            return "\t".join([peptide.sequence, parent_proteins, peptide.mass]) + "\n"
        else:
            return "\t".join([peptide.sequence, parent_proteins]) + "\n"

    def write(self, outdir=None, include_mass=False):
        """
        Write results to a file.

        Parameters
        ----------
        outdir : str, optional
            Output directory where to save the peptides list and attributes
        include_mass : bool, default=False
            Wheter or not to include peptide mass in the tsv file
        """
        if not self._is_cleaved():
            self.cleave_proteins()
        if outdir is not None:
            self.set_outdir(outdir)

        with open(self.outdir, "w", encoding="UTF-8") as out_file:
            for peptide in self.peptides:
                out_file.write(self._format_line(peptide, include_mass))

    def cleave_proteins(self):
        """
        Perform an in silico digestion of proteins with previously selected
        or default parameters

        Returns
        -------
        None
            Update the self.peptides list
        """
        enzyme = Enzyme(self.params["enzyme"])
        max_miscleavages = self.params["max_miscleavages"]
        clip_nterm_m = self.params["clip_nterm_m"]
        all_peptides = {}

        for id_, protein_sequence in self.proteins.items():
            # Digestion protein sequence and get peptides
            peptides = [
                Peptide(sequence, id_)
                for sequence in enzyme.digest(
                    protein_sequence, max_miscleavages, clip_nterm_m
                )
            ]
            valid_peptides = self._get_valid_peptides(peptides)

            # Update peptides' dictionary
            for peptide in valid_peptides:
                sequence = peptide.sequence
                # Update parent protein list if peptide already exist
                if sequence in all_peptides:
                    all_peptides[sequence].add_parent_protein(id_)
                # If peptide not yet in dictionary, add it
                else:
                    all_peptides[sequence] = peptide

        self.peptides = list(all_peptides.values())

    def get_uniques_peptides(self):
        """
        Get unique peptides.

        Return
        ------
        list of str
            List of unique peptide sequences.
        """
        return [peptide.sequence for peptide in self.peptides if peptide.is_unique()]

    def get_shared_peptides(self):
        """
        Get shared peptides (i.e. peptide with multiple parent proteins)

        Return
        ------
        list of str
            List of shared peptide sequences.
        """
        return [
            peptide.sequence for peptide in self.peptides if not peptide.is_unique()
        ]

    def get_peptide_degeneracy(self):
        """
        Get peptide degeneracy dict.

        Return
        ------
        dict
            A dictionary with peptide sequences as keys and boolean as key (True
            if unique, False if not).
        """
        return {peptide.sequence: peptide.is_unique() for peptide in self.peptides}


# ------------------------------------------------------------------------------
# For command line usage

if __name__ == "__main__":
    args = cli.parse_args()

    run = Digestion(args.fasta_file)
    run.set_min_length(args.min_length)
    if args.max_length is not None:
        run.set_max_length(args.max_length)
    run.set_max_miscleavages(args.miscleavages)
    run.set_max_mass(args.mass)
    run.set_clip_nterm_m(args.clip_nterm_m)
    run.set_enzyme(args.enzyme)
    if args.output is not None:
        run.set_outdir(args.output)

    print(
        dedent(f"""
                Performing digestion with: 
                
                min_length = {run.params["min_length"]}
                max_length = {run.params["max_length"]}
                enzyme = '{run.params["enzyme"]}'
                max_mass = {run.params["max_mass"]}
                max_miscleavages = {run.params["max_miscleavages"]}
                """)
    )

    run.cleave_proteins()
    if args.unique is True:
        run.check_peptide_uniqueness()
    run.write()
    print(f"List of peptide sequences saved as {run.outdir}")
