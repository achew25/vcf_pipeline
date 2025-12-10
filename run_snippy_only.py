import argparse
import subprocess
import os
import sys

# Configuration
# Paths to the external tools.
BWA_TOOL = "bwa"
SAMTOOLS_TOOL = "samtools"
FREEBAYES_TOOL = "freebayes"

# Directory Constants
INPUT_DIR = "data/raw"
OUTPUT_DIR = "data/snippy_only"


def run_command(command, step_name):
    """Helper function to execute a command and handle errors."""
    print(f"\n>>> Running Step: {step_name}")
    print(f"Command: {' '.join(command)}")
    try:
        # Check=True will raise an exception for non-zero exit codes
        # We redirect stdout/stderr to the console for real-time feedback
        subprocess.run(command, check=True, stdout=sys.stdout, stderr=sys.stderr)
    except subprocess.CalledProcessError as e:
        print(f"\n!!! ERROR in {step_name} !!!")
        print(f"Command failed with exit code {e.returncode}")
        sys.exit(1)
    except FileNotFoundError:
        print(f"\n!!! ERROR: Required tool not found for {step_name} !!!")
        print("Please ensure BWA, SAMTOOLS, and FREEBAYES are installed and in your PATH.")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="A simple Bash-based pipeline for reference mapping and variant calling.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument(
        "-r", "--reference", required=True, 
        help=f"Filename of the reference genome FASTA file (.fa/.fasta) within the '{INPUT_DIR}' directory."
    )
    parser.add_argument(
        "-1", "--read1", required=True, 
        help=f"Filename of the R1 (forward) sequencing reads FASTQ file (.fq/.fastq) within the '{INPUT_DIR}' directory."
    )
    parser.add_argument(
        "-2", "--read2", required=False, default=None,
        help=f"Filename of the R2 (reverse) sequencing reads FASTQ file (.fq/.fastq) within the '{INPUT_DIR}' directory. (Optional)"
    )
    parser.add_argument(
        "-o", "--output_prefix", required=True, 
        help="Prefix for all output files (e.g., 'sample_A'). Output files go to the '{OUTPUT_DIR}' directory."
    )
    
    args = parser.parse_args()

    # Setup Directories
    if not os.path.isdir(INPUT_DIR):
        print(f"!!! WARNING: Input directory '{INPUT_DIR}' not found. Please create it and place your files inside.")
    
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Output directory '{OUTPUT_DIR}' ensured to exist.")


    # Input Path Resolution
    REF_FASTA = os.path.join(INPUT_DIR, args.reference)
    R1_FASTQ = os.path.join(INPUT_DIR, args.read1)
    
    if args.read2:
        R2_FASTQ = os.path.join(INPUT_DIR, args.read2)
        print("Running in PAIR-END mode.")
    else:
        R2_FASTQ = None
        print("Running in SINGLE-END mode (No --read2 provided).")

    #Output Path Resolution 
    BAM_UNSORTED = os.path.join(OUTPUT_DIR, f"{args.output_prefix}.unsorted.bam")
    BAM_SORTED = os.path.join(OUTPUT_DIR, f"{args.output_prefix}.sorted.bam")
    VCF_RAW = os.path.join(OUTPUT_DIR, f"{args.output_prefix}.raw.vcf")
    
    # 1. Index the Reference Genome (BWA)
    run_command([BWA_TOOL, "index", REF_FASTA], "1. BWA Indexing")

    # 2. Align Reads to Reference (BWA MEM)
    print("\n>>> Running Step: 2. BWA Alignment (Piped)")
    
    bwa_cmd = [BWA_TOOL, "mem", "-t", "4", REF_FASTA, R1_FASTQ] # -t for threads
    if R2_FASTQ:
        bwa_cmd.append(R2_FASTQ)
        
    # Open the output file for the unsorted BAM
    with open(BAM_UNSORTED, 'wb') as bam_out:
        # Start BWA MEM process
        bwa_process = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, stderr=sys.stderr)
        
        # Start SAMTOOLS view process (to convert SAM to BAM)
        samtools_view_cmd = [SAMTOOLS_TOOL, "view", "-b", "-"] # -b for BAM, - for stdin
        samtools_view_process = subprocess.Popen(samtools_view_cmd, stdin=bwa_process.stdout, stdout=bam_out, stderr=sys.stderr)
        
        # Close the pipe from bwa_process to allow it to finish
        bwa_process.stdout.close()
        
        # Wait for both processes to finish
        samtools_view_process.wait()
        bwa_process.wait()
        
        # Check for errors
        if bwa_process.returncode != 0 or samtools_view_process.returncode != 0:
            print("\n!!! ERROR in BWA/SAMTOOLS Piped Alignment !!!")
            sys.exit(1)
            
    print("Alignment and SAM to BAM conversion complete.")

    # 3. Sort and Index the BAM file (SAMTOOLS)
    run_command([SAMTOOLS_TOOL, "sort", "-o", BAM_SORTED, BAM_UNSORTED], "3. SAMTOOLS Sort")
    run_command([SAMTOOLS_TOOL, "index", BAM_SORTED], "4. SAMTOOLS Index")
    
    # 5. Variant Calling (FREEBAYES)
    run_command(
        [
            FREEBAYES_TOOL,
            "-f", REF_FASTA,
            BAM_SORTED,
            "-v", VCF_RAW,
        ], 
        "5. Freebayes Variant Calling"
    )

    print("\n--- Pipeline Complete ---")
    print(f"Final BAM file: {BAM_SORTED}")
    print(f"Final VCF file: {VCF_RAW}")

if __name__ == "__main__":
    main()