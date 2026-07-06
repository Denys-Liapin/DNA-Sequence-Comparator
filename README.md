# DNA-Sequence-Comparator

## Overview
A versatile Python command-line interface tool designed to compare DNA sequences, identify single nucleotide polymorphisms, and visualize mutational contexts. 

## Scientific Context
This tool is built as a universal sequence comparator. Whether analyzing viral strains, bacterial adaptations, or genetic variations in plants and mammals, this script allows researchers to quickly spot structural differences and potential functional mutations without relying on resource-heavy, genome-wide aligners. It is designed to be sequence-agnostic, providing fast, localized visualization of SNPs in any provided FASTA sequence.

## Usage
Run the script from your terminal by passing your reference and sample `.fasta` files directly:

```bash
git clone https://github.com/Denys-Liapin/DNA-Sequence-Comparator.git
cd DNA-Sequence-Comparator
python3 dna_comparator.py <reference_file.fasta> <sample_file.fasta>
```

## Key Features
- **FASTA Parsing:** Efficiently reads standard biological sequence files using `Biopython`.
- **Contextual Visualization:** Generates clean, color-coded terminal outputs to highlight mutations within their immediate genetic neighborhood.
- **Variant Counting:** Automatically calculates the total number of sequence discrepancies.

## Tech Stack
- Python 3
- Biopython

<details>
<summary>

work example

</summary>

<img width="977" height="483" alt="image" src="https://github.com/user-attachments/assets/b4c6e664-47de-4646-b7de-42444bf72345" />

</details>
