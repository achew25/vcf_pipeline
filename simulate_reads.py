# simulate_reads.py

import os
import argparse
import sys
# Import functions
from functions.get_depth import get_depth
from functions.create_fastq import create_fastq

def simulate_reads(args):
    """Reads mutated sequence from a file in 'check_intermediary_file', generates reads, and creates FASTQ output."""

    # Define intermediary directory and file path
    intermediary_dir = "check_intermediary_file"
    intermediary_path = os.path.join(intermediary_dir, args.input_seq)

    # 1. Read Mutated Sequence from TXT file
    print(f"Reading mutated sequence from: {intermediary_path}", file=sys.stderr)
    
    try:
        with open(intermediary_path, 'r') as f:
            mutated_sequence = f.readline().strip()
    except FileNotFoundError:
        print(f"Error: Input sequence file not found at {intermediary_path}. Did Part 1 run successfully?", file=sys.stderr)
        sys.exit(1)
    except IOError as e:
        print(f"Error reading input file {intermediary_path}: {e}", file=sys.stderr)
        sys.exit(1)

    if not mutated_sequence:
        print("Error: Input sequence file is empty or missing data.", file=sys.stderr)
        sys.exit(1)
        
    # Use the provided tag for output file naming
    database_tag = [args.fasta_tag] 
    print(f"Processing sequence tag: {database_tag[0]}", file=sys.stderr)

    # 2. Get Fragmented Reads (Simulate Sequencing Depth)
    print(f"Generating fragmented reads with avg length {args.read_len} at {args.depth} depth...", file=sys.stderr)
    mutated_fragments = get_depth(mutated_sequence, args.depth, args.read_len)

    # 3. Create FASTQ Output
    fastq_filename = database_tag[0] + ".fastq"
    os.makedirs(args.outdir, exist_ok=True) 
    file_path_and_name = os.path.join(args.outdir, fastq_filename)
    
    print(f"Writing FASTQ file to: {file_path_and_name}", file=sys.stderr)

    create_fastq(
        file_path_and_name, database_tag, args.which_fasta, mutated_fragments, args.quality_ascii
    )
    print("FASTQ simulation complete!", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="Part 2: Read mutated sequence from 'check_intermediary_file', simulate reads, and output FASTQ file.")
    parser.add_argument("fasta_tag", help="The FASTA ID tag for the sequence being processed (e.g., NC_000913.3).")
    parser.add_argument("--input-seq", required=True, help="Input file name for the mutated sequence (e.g., mutated_sequence.txt).")
    parser.add_argument("--outdir", default="data/processed", help="Output directory for FASTQ file. Default: data/processed")
    parser.add_argument("--depth", type=int, default=30, help="Depth or coverage of reads on reference genome. Default: 30")
    parser.add_argument("--read_len", type=int, default=100, help="Average fragment/read length. Default: 100")
    parser.add_argument("--which_fasta", type=int, default=0, help="Index of the sequence in the FASTA file used in Part 1. Default: 0")
    parser.add_argument("--quality_ascii", type=str, default="~", help="Writes ascii quality score as ~ for perfect. Default: ~")
    
    args = parser.parse_args()
    simulate_reads(args)

if __name__ == "__main__":
    main()