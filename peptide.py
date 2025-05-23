from mass import amino_acid, atom


class Peptide:
    def __init__(self, sequence, id_):
        self.sequence = sequence
        self.mass = self.get_mass()
        self.parent_proteins = [id_]

    def __len__(self):
        return len(self.sequence)

    def __repr__(self):
        return f"Peptide({self.sequence})"

    def get_mass(self):
        """
        Sum the mass of all residues and add the mass of H (N-terminal) and OH
        (C-terminus) to obtain the sequence mass
        """
        _mass = (
            sum((amino_acid[aa] for aa in self.sequence)) + atom["H"] * 2 + atom["O"]
        )
        return round(_mass, 4)

    def add_parent_protein(self, id_):
        """
        Add protein id to the parent proteins list.

        Parameters
        ----------
        id_ : str
            Protein id to be added.

        Returns
        -------
        None
            Update self.parent_proteins.
        """
        self.parent_proteins.append(id_)

    def is_unique(self):
        """
        Check if peptide is unique (found in only one protein from database)

        Returns
        -------
        bool
            Wheter or not peptide is unique.
        """
        return len(self.parent_proteins) == 1
