from Bio import SeqIO
from Bio.Seq import Seq

def load_and_check_fasta(reference_file, sample_file):
    print("--- Running DNA Analyzer ---")
    
    try:
        # Read files containing a single biological sequence each
        ref_record = SeqIO.read(reference_file, "fasta")
        sample_record = SeqIO.read(sample_file, "fasta")
        
        # Extract sequences and convert to uppercase for reliability
        ref_seq = ref_record.seq.upper()
        sample_seq = sample_record.seq.upper()
        
        print(f"Reference loaded. Length: {len(ref_seq)} bp (base pairs)")
        print(f"Sample loaded. Length: {len(sample_seq)} bp")
        
        # Sequence length validation
        if len(ref_seq) != len(sample_seq):
            print("Warning: Sequences have different lengths! Direct comparison is impossible.")
            return None, None
        else:
            print("Lengths match, proceeding with the analysis.\n")
            return ref_seq, sample_seq
            
    except Exception as error:
        print(f"Error reading FASTA files: {error}")
        return None, None

def analyze_dna_sequences(ref_dna, sample_dna):
    if ref_dna is None or sample_dna is None:
        return False

    if ref_dna == sample_dna:
        print("DNA sequences are absolutely identical.")
        return False
    else:
        print("Discrepancies detected in DNA strands. Starting point mutation search:")
        find_and_display_mutations(ref_dna, sample_dna)
        return True

def find_and_display_mutations(ref_dna, sample_dna):
    mutations_count = 0
    
    for index, (nuc_ref, nuc_sample) in enumerate(zip(ref_dna, sample_dna)):
        if nuc_ref != nuc_sample:
            mutations_count += 1
            show_mutation_context(ref_dna, sample_dna, mutation_pos=index)
            check_codon_mutation(ref_dna, sample_dna, mutation_pos=index)
            
    print("\n" + "="*40)
    print(f"Total point mutations found: {mutations_count}")

def show_mutation_context(ref_dna, sample_dna, mutation_pos, window=10):
    # Define the boundaries of the sliding window around the mutation
    start = max(0, mutation_pos - window)
    end = min(len(ref_dna), mutation_pos + window + 1)
    
    ref_chunk = ref_dna[start:end]
    sample_chunk = sample_dna[start:end]
    
    # Build the graphical alignment match line
    match_line = ""
    for nuc_r, nuc_s in zip(ref_chunk, sample_chunk):
        if nuc_r == nuc_s:
            match_line += "|"
        else:
            match_line += "*" # Mismatch indicator
            
    # ANSI escape codes for terminal color formatting
    RED_COLOR = '\033[91m'
    RESET_COLOR = '\033[0m'
    
    # Highlight the mutated nucleotide in red
    chunk_index = mutation_pos - start
    colored_sample_chunk = (
        sample_chunk[:chunk_index] + 
        RED_COLOR + sample_chunk[chunk_index] + RESET_COLOR + 
        sample_chunk[chunk_index+1:]
    )

    print(f"\nMutation Context (DNA position {mutation_pos + 1}):")
    print(f"Reference : {ref_chunk}")
    print(f"Match     : {match_line}")
    print(f"Sample    : {colored_sample_chunk}")

def check_codon_mutation(ref_dna, sample_dna, mutation_pos):
    # Determine the index and 1-based position of the amino acid
    amino_acid_index = mutation_pos // 3
    amino_acid_number = amino_acid_index + 1
    
    # Find the boundaries of the corresponding triplet (codon)
    codon_start = amino_acid_index * 3
    codon_end = codon_start + 3
    
    ref_codon = ref_dna[codon_start:codon_end]
    sample_codon = sample_dna[codon_start:codon_end]
    
    # Ignore incomplete codons at the end of the sequence
    if len(ref_codon) < 3 or len(sample_codon) < 3:
        return

    # Translate the specific triplet into an amino acid
    amino_acid_ref = Seq(ref_codon).translate()
    amino_acid_sample = Seq(sample_codon).translate()
    
    # Check if the DNA mutation causes a missense/nonsense change in the protein structure
    if amino_acid_ref != amino_acid_sample:
        codon_nucleotide_pos = codon_start + 1
        print(f"-> Consequence: Missense mutation in amino acid #{amino_acid_number} "
              f"(Codon starts at nucleotide position {codon_nucleotide_pos}): "
              f"Reference: {amino_acid_ref} -> Sample: {amino_acid_sample}")


# Main execution block
reference_sequence, sample_sequence = load_and_check_fasta("tea_etalon.fasta", "tea_reducted.fasta")
is_mutated = analyze_dna_sequences(reference_sequence, sample_sequence)