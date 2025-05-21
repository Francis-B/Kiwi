#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
"""
This module contains the Enzyme class that is used by the Digestion class to
perform the in silico digestion of the protein sequences.
"""

# ------------------------------------------------------------------------------
import re


class Enzyme:
    # Regular expression used to find the cleavage site
    cleavage_site = {
        "trypsin-p": "(?<=[RK](?!P))",  # Cleaves after arginine or lysine not followed by proline
        "trypsin": "(?<=[RK])",  # Cleaves after arginine or lysine
        "lysc-p": "(?<=K(?!P))",  # Cleaves after lysine not followed by proline
        "lysc": "(?<=K)",  # Cleaves after lysine
        "lysn": "(?=K)",  # Cleaves before lysine
        "argc": "(?<=R)",  # Cleaves after arginine
        "aspc": "(?<=D)",  # Cleaves after aspartic acid
        "aspn": "(?=D)",  # Cleaves before aspartic acid
        "gluc": "(?<=E)",  # Cleaves after glutamic acid
        "glun": "(?=E)",  # Cleaves before glutamic acid
        "chymotrypsin+": "(?<=[YWFLM])",  # Cleaves after tyrosine, tryptophane, phenylalanine, leucine, methionine
        "chymotrypsin": "(?<=[YWF])",  # Cleaves after tyrosine, tryptophane, phenylalanine
    }

    def __init__(self, enzyme):
        self.digestion_site = self.__class__.cleavage_site[enzyme]

    @classmethod
    def get_available_enzyme(cls):
        """
        Get all available enzyme.

        Returns
        -------
        list
            Name of each enzyme that is available for the in silico digestion.
        """
        print(list(cls.cleavage_site.keys()))

    def _join_sequences(self, base_sequences, miss):
        """
        Join adjacent peptide sequences to simulate miscleavages.

        Parameters
        ----------
        base_sequences : list
            Sequence of fully digested peptides.
        miscleavages : int
            Number of miscleavages. A value of 1 would mean that 2 base_sequences
            will be joined; 2 would mean that 3 will be joined; etc.
        """
        num_sequences = len(base_sequences) - miss  # Number of sequences to loop on
        return ["".join(base_sequences[n : n + miss + 1]) for n in range(num_sequences)]

    def digest(self, sequence, miscleavages, clip_nterm_m):
        """
        Cleave protein sequences.

        Parameters
        ----------
        sequence : str
                Protein sequence to digest.
        miscleavages : int
                Number of maximum miscleavages allowed

        Returns
        -------
        list
            All peptide sequences obtained with the digestion.
        """

        base_sequences = re.split(self.digestion_site, sequence)

        # Initiate peptides list with all the fully digested peptide
        peptide_sequences = base_sequences

        miss = 1
        while miss <= miscleavages:
            joined_sequences = self._join_sequences(base_sequences, miscleavages)

            # If there is no joined sequences, break out of the loop
            if joined_sequences == []:
                break
            # Also consider N-terminal peptide without its methionine
            if not clip_nterm_m:
                first_pep = joined_sequences[0]
                joined_sequences.append(re.sub(r"^[Mm]", "", first_pep))

            peptide_sequences.extend(joined_sequences)

            miss += 1

        return peptide_sequences
