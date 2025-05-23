from digestion import Digestion

digestion = Digestion("data/uniprot_SwissProt_Human_1_11_2017.fasta")

digestion.set_max_miscleavages(2)
digestion.cleave_proteins()

protein = list(digestion.proteins.values())[0]

protein.get_peptides()

"abcdefgh".find("cd")

positions = []
for peptide in protein.get_peptides():
    start = protein.sequence.find(peptide.sequence)
    end = start + len(peptide)
    positions.append((start, end))

print(positions)
