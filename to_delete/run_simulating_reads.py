# run_pipeline.py

import os
import argparse
# Import functions
from functions.parse_fasta import parse_fasta
from functions.make_mutations import make_mutations
from functions.get_depth import get_depth
from functions.create_fastq import create_fastq

def run_simulating_reads_pipeline(args):
    """Executes the pipeline."""

    # 1. Parse FASTA Input
    print(f"Parsing FASTA: {args.fasta_file}...")
    file_path = os.path.join("data", "raw", args.fasta_file)
    sequence_to_mutate_dict = parse_fasta(file_path)

    # Use the first sequence in the dictionary for mutation (as per your original script)
    fasta_id, sequence_to_mutate = list(sequence_to_mutate_dict.items())[args.which_fasta]
    database_tag = fasta_id.strip().split(" ", maxsplit=1)
    print(f"Target Sequence: {database_tag[0]}")

    # 2. Make Mutations
    print(f"Introducing {args.snps} SNPs and {args.indels} INDELs...")
    mutated_sequence = make_mutations(
        sequence_to_mutate, args.snps, args.indels, args.max_indel_len
    )

    # 3. Get Fragmented Reads (Simulate Sequencing Depth)
    print(f"Generating fragmented reads with avg length {args.read_len} at {args.depth} depth...")
    mutated_fragments = get_depth(mutated_sequence, args.depth, args.read_len)

    # 4. Create FASTQ Output
    fastq_filename = database_tag[0] + ".fastq"
    os.makedirs(args.outdir, exist_ok=True) # Ensure output directory exists
    file_path_and_name = os.path.join(args.outdir, fastq_filename)
    
    print(f"Writing FASTQ file to: {file_path_and_name}")
    # Assume perfect reads ascii = "~"

    create_fastq(
        file_path_and_name, database_tag, args.which_fasta, mutated_fragments, args.quality_ascii
    )
    print("SNP and INDEL simulation complete!")


def main():
    parser = argparse.ArgumentParser(description="Simulate single strand illumina reads from a reference genome")
    parser.add_argument("fasta_file", help="Input FASTA file name (e.g., EcoliK12-MG1655.fasta)")
    parser.add_argument("--outdir", default="data/processed", help="Output directory for FASTQ file. Default: data/processed")
    parser.add_argument("--snps", type=int, default=300, help="Number of SNPs to introduce. Default: 300")
    parser.add_argument("--indels", type=int, default=10, help="Number of INDELs to introduce. Default: 10")
    parser.add_argument("--depth", type=int, default=30, help="Depth or coverage of reads on reference genome. Default: 30")
    parser.add_argument("--max_indel_len", type=int, default=10, help="Maximum length of an INDEL. Default: 10")
    parser.add_argument("--read_len", type=int, default=100, help="Average fragment/read length. Default: 100")
    parser.add_argument("--which_fasta", type=int, default=0, help="Index of the sequence in the FASTA file to use. Default: 0 (first entry)")
    parser.add_argument("--quality_ascii", type=str, default="~", help="Writes ascii quality score as ~ for perfect. Default: ~")
    
    args = parser.parse_args()
    run_simulating_reads_pipeline(args)

if __name__ == "__main__":
    main()