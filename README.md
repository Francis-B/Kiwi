# Kiwi: a fast regex-based protein digestion tool

## Description

Kiwi is a python tool that performs a proteolytic cleavage on proteins from a given fasta and that can determine which obtained peptides are unique.

The digestion can be tuned with the following parameters:

- Enzyme used to perform the digestion.
- Minimal and maximal length of the peptide sequences;
- Maximal molecular mass of the peptide sequences;
- Number of miscleavages allowed;

## Instructions

To use it, clone or download this repository.

**To run in python script:**

From a script in the Kiwi folder:

```python
from digestion import Digestion
#-------------------------------------------------------------------------------

fastaFile = '/path/to/file.fasta'
digestion = Digestion(fastaFile) # Create an instance with default parameters and the
                                 # protein sequences dictionary

# Change parameters
digestion.set_min_length(<int>)
digestion.set_max_length(<int>)
digestion.set_max_mass(<float>)
digestion.set_max_miscleavages(<int>)
digestion.set_outdir(<filepath>)
digestion.set_enzyme(<str>)  # Implemented enzymes can be found in enzyme.py

# Print parameters
print(digestion)

# Running
digestion.cleave_proteins()
digestion.check_sequences_uniqueness()
digestion.write()  # write the result into a file

# Access the list of sequences
digestion.peptides
```

By default, the peptide sequences returned are at least 7 amino acids long, have a maximum of 1 miscleavage, have a molecular mass under 4600 dalton and result from a tryptic digestion.

If no directory is passed in digestion.write(), the file is
automatically saved as /path/to/file.fasta_digestedPeptides.csv.

**To run from command line:**

Linux users can run digestool from the terminal by typing:

```bash
/path/to/digestool.py /path/to/file.fasta
```

To show the help message:

```bash
$ /path/to/digestion.py -h

  usage: digestion.py [-h] [-l] [-M] [-m] [-a] [-e] [-o] </path/to/file.fasta>

  Digest proteins of a given fasta file

  positional arguments:
    </path/to/file.fasta>

  optional arguments:
    -h, --help            show this help message and exit
    -l , --length         minimal length of peptide sequences (default: 7)
    -M , --miscleavages   maximum of miscleavages allowed (default: 1)
    -m , --mass           maximal molecular mass of peptide sequences (default: 4600 dalton)
    -u , --unique         the list returned will contain a flag to inform if the peptides are unique or not
    -e , --enzyme         enzyme used to perform digestion (default: trypsin)
    -o , --output         output directory (default: /path/to/file.fasta_digestedPeptides.csv)
```

## Note on unique peptides

This script was mainly developped has a tool for mass spectrometry database preparation. In this context, we consider a peptide sequence as unique if no other proteins yield the **exact** same sequence after enzymatic digestion.

So, for example, if two given proteins yield respectively the following peptides after a tryptic digestion:

  1. EIQILLR
  2. APELDFGEIQILLR

Even if the first peptide can be found in the second one, it is still consider has unique.

## Requirements

python >= 3.6
