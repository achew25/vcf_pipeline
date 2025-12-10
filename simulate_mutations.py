# simulate_mutations.py

import os
import argparse
import sys
# Import functions
from functions.parse_fasta import parse_fasta
from functions.make_mutations import make_mutations 
# NOTE: The ground truth writing function is now integrated within make_mutations.py

def simulate_mutations(args):
    """
    Parses FASTA, introduces mutations, and writes the mutated sequence 
    to a file, and outputs the ground truth variants to a separate file.
    """

    # Define intermediary directory and file paths
    intermediary_dir = "check_intermediary_file"
    intermediary_seq_path = os.path.join(intermediary_dir, args.output_seq)
    
    # NEW: Define path for the Ground Truth file
    intermediary_gt_path = os.path.join(intermediary_dir, args.output_gt)
    
    # Ensure the intermediary directory exists
    print(f"Ensuring intermediary directory '{intermediary_dir}' exists...", file=sys.stderr)
    os.makedirs(intermediary_dir, exist_ok=True)


    # 1. Parse FASTA Input
    print(f"Parsing FASTA: {args.fasta_file}...", file=sys.stderr)
    file_path = os.path.join("data", "raw", args.fasta_file)
    try:
        sequence_to_mutate_dict = parse_fasta(file_path)
    except FileNotFoundError:
        print(f"Error: FASTA file not found at {file_path}", file=sys.stderr)
        sys.exit(1)

    # Use the specified sequence in the dictionary
    try:
        fasta_id, sequence_to_mutate = list(sequence_to_mutate_dict.items())[args.which_fasta]
    except IndexError:
        print(f"Error: Index {args.which_fasta} out of range for sequences in FASTA file.", file=sys.stderr)
        sys.exit(1)
        
    database_tag = fasta_id.strip().split(" ", maxsplit=1)
    print(f"Target Sequence: {database_tag[0]}", file=sys.stderr)

    # 2. Make Mutations
    # NOTE: The make_mutations function is responsible for calling 
    # write_ground_truth with the path provided below.
    print(f"Introducing {args.snps} SNPs and {args.indels} INDELs...", file=sys.stderr)
    
    # To properly integrate the ground truth writing into the file system,
    # the make_mutations function *must* be updated to accept the output path 
    # for the ground truth file. Assuming you have updated `make_mutations` 
    # to accept `output_gt_path` as a fifth argument:
    
    mutated_sequence = make_mutations(
        sequence_to_mutate, 
        args.snps, 
        args.indels, 
        args.max_indel_len,
        intermediary_gt_path # NEW: Pass the ground truth file path
    )
    # If you updated `make_mutations` to use a fixed filename like "ground_truth_variants.txt" 
    # and moved the file handling logic (like `write_ground_truth`) into it, you don't need this fifth argument. 
    # The VCF-like formatting is crucial for evaluation. 

    # 3. Output Mutated Sequence to TXT file
    print(f"Writing mutated sequence to: {intermediary_seq_path}", file=sys.stderr)
    try:
        with open(intermediary_seq_path, 'w') as f:
            f.write(mutated_sequence + '\n')
    except IOError as e:
        print(f"Error writing to output file {intermediary_seq_path}: {e}", file=sys.stderr)
        sys.exit(1)

    print("Mutation simulation complete! Check both sequence and ground truth files.", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="Part 1: Simulate mutations and write mutated sequence to a .txt file inside 'check_intermediary_file'.")
    parser.add_argument("fasta_file", help="Input FASTA file name (e.g., EcoliK12-MG1655.fasta)")
    parser.add_argument("--output-seq", required=True, help="Output file name for the mutated sequence (e.g., mutated_sequence.txt).")
    
    # NEW ARGUMENT: Output file for the ground truth
    parser.add_argument("--output-gt", required=True, help="Output file name for the ground truth variants (e.g., ground_truth.txt).")
    
    parser.add_argument("--snps", type=int, default=300, help="Number of SNPs to introduce. Default: 300")
    parser.add_argument("--indels", type=int, default=10, help="Number of INDELs to introduce. Default: 10")
    parser.add_argument("--max_indel_len", type=int, default=10, help="Maximum length of an INDEL. Default: 10")
    parser.add_argument("--which_fasta", type=int, default=0, help="Index of the sequence in the FASTA file to use. Default: 0")
    
    args = parser.parse_args()
    simulate_mutations(args)

if __name__ == "__main__":
    main()