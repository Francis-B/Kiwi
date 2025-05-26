from digestion import Digestion

digestion = Digestion("data/uniprot_SwissProt_Human_1_11_2017.fasta")

digestion.set_max_miscleavages(2)
digestion.cleave_proteins()

protein = list(digestion.proteins.values())[0]

protein.get_peptides()

len(protein)
protein.get_detectable_loc("all")
