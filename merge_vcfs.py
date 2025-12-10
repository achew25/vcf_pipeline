import argparse
import subprocess
import os
import sys

def run_command(command):
    """Executes a shell command and checks for errors."""
    print(f"Executing: {' '.join(command)}")
    try:
        # Run the command, capture stdout/stderr, and check the return code
        result = subprocess.run(
            command,
            check=True,  # Raise a CalledProcessError if the exit code is non-zero
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True  # Decode stdout/stderr as text
        )
        # print(f"Output: {result.stdout}")
    except subprocess.CalledProcessError as e:
        print(f"\n--- ERROR executing command: {' '.join(command)} ---", file=sys.stderr)
        print(f"STDOUT:\n{e.stdout}", file=sys.stderr)
        print(f"STDERR:\n{e.stderr}", file=sys.stderr)
        # Exit the script upon failure
        sys.exit(1)
    except FileNotFoundError:
        # Catch errors if the executable (bgzip, tabix, bcftools) is not found
        print(f"\n--- ERROR: Executable not found. Check if bgzip, tabix, and bcftools are in your PATH. ---", file=sys.stderr)
        sys.exit(1)

def compress_and_merge_vcfs(input_vcf_paths, output_vcf_path):
    """
    Compresses input VCF files, indexes them, and then merges them.
    """
    compressed_files = []

    # --- 1. Compress and Index Each Input VCF File ---
    print("\n### Starting Compression and Indexing ###")
    for original_vcf in input_vcf_paths:
        if not os.path.exists(original_vcf):
            print(f"Warning: Input file not found: {original_vcf}. Skipping.", file=sys.stderr)
            continue

        compressed_vcf = f"{original_vcf}.gz"
        compressed_files.append(compressed_vcf)

        # 1a. Compress using bgzip (keeps the original VCF with -k)
        # Assuming the VCF is sorted, otherwise add a bcftools sort step here.
        bgzip_cmd = ["bgzip", "-k", original_vcf]
        run_command(bgzip_cmd)
        
        # 1b. Index using tabix (or bcftools index)
        # -p vcf specifies the VCF file type
        tabix_cmd = ["tabix", "-p", "vcf", compressed_vcf]
        run_command(tabix_cmd)
        
    if not compressed_files:
        print("\n--- ERROR: No VCF files were successfully compressed. Exiting. ---", file=sys.stderr)
        sys.exit(1)

    # --- 2. Merge all Compressed VCFs using bcftools merge ---
    print("\n### Starting bcftools merge ###")
    
    # Base command: bcftools merge [options] file1.vcf.gz file2.vcf.gz ...
    merge_cmd = [
        "bcftools",
        "merge",
        *compressed_files, # Unpacks the list of compressed file names
        "-Oz",             # Output compressed VCF
        "-o", output_vcf_path
        # You can add more bcftools merge options here, e.g., "-i", "AC:sum"
    ]
    
    run_command(merge_cmd)
    
    print("\nMERGE COMPLETE!")
    print(f"Final merged VCF saved to: {output_vcf_path}")

def main():
    parser = argparse.ArgumentParser(
        description="A script to compress VCF files with bgzip/tabix and then merge them with bcftools merge."
    )
    # Use 'nargs='+' to accept one or more input VCF files
    parser.add_argument(
        "input_vcf_paths", 
        nargs='+', 
        help="One or more paths to the VCF files to be compressed and merged."
    )
    parser.add_argument(
        "-o", "--output", 
        required=True, 
        help="The path and filename for the final merged VCF output (e.g., merged_data.vcf.gz)."
    )

    args = parser.parse_args()
    
    # Ensure the output file name ends with the expected compression suffix
    if not args.output.endswith(".vcf.gz"):
        print("Warning: Output file name should end with .vcf.gz for proper compression.", file=sys.stderr)
        
    compress_and_merge_vcfs(args.input_vcf_paths, args.output)

if __name__ == "__main__":
    main()