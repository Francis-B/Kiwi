import matplotlib.pyplot as plt


class Protein:
    def __init__(self, id_, sequence):
        self.id = id_
        self.sequence = sequence
        self.peptides = []

    def get_id(self):
        return self.id

    def add_peptide(self, peptide):
        self.peptides.append(peptide)

    def add_peptides(self, peptides):
        self.peptides.extend(peptides)

    def get_peptides(self):
        """
        Get all valid peptides of protein.

        Returns
        -------
        List of Peptide
            List of all valid peptides of the proteins for the given digestion
            parameters.
        """
        return self.peptides

    def get_unique_peptides(self):
        """
        Get all valid unique peptides of protein.

        Returns
        -------
        List of Peptide
            List of all valid and unique peptides of the proteins.
        """
        return [peptide for peptide in self.peptides if peptide.is_unique()]

    def get_shared_peptides(self):
        """
        Get all valid non-unique peptides of the protein.

        Returns
        -------
        List of Peptide
            List of all valid and non-unique peptides of the proteins.
        """
        return [peptide for peptide in self.peptides if not peptide.is_unique()]

    def get_positions(self, peptide_type):
        """
        Get the position of all peptides. Position of adjacent peptides will be
        merged.

        Parameters
        ----------
        peptide_type : str
            Type of peptide for which to position are to be retrieved. Options
            are "all", "unique" and "shared".
        """
        if peptide_type == "all":
            peptides = self.get_peptides()
        elif peptide_type == "unique":
            peptides = self.get_unique_peptides()
        elif peptide_type == "shared":
            peptides = self.get_shared_peptides()
        else:
            print(peptide_type, " is not a valid option.")
            return

        positions = []
        for peptide in peptides:
            start = self.sequence.fint(peptide.sequence)
            end = start + len(peptide)
            positions.append((start, end))

        # TODO: Concatenate all adjacent positions
        return positions

    def plot_peptides_map(self):
        """
        Plot a map of the valid peptides of this proteins.

        Returns
        -------
        Matplotlib.Figure
            Map of the valid peptides of this proteins.
        """
        pass
