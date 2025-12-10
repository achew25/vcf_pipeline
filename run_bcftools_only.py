# run_alignment_and_variant_caller.py

import subprocess
import argparse
import os

def run_bcftools_alignment_and_variant_calling(args):
    """
    Executes the alignment and variant calling workflow using external tools,
    optionally supporting paired-end data.
    """
    
    # File Names and Paths
    fasta_path = os.path.join(args.data_raw_dir, args.fasta_name)
    
    # Required: FASTQ R1 path
    fastq_r1_path = os.path.join(args.data_processed_dir, args.fastq_r1_name)
    
    # Optional: FASTQ R2 path
    fastq_r2_path = None
    if args.fastq_r2_name:
        fastq_r2_path = os.path.join(args.data_processed_dir, args.fastq_r2_name)
        mode = "Paired-End"
    else:
        mode = "Single-End"
    
    # Output file names (based on FASTA name)
    prefix = os.path.splitext(args.fasta_name)[0]
    bam_output_name = f"{prefix}.sorted.bam"
    vcf_output_name = f"{prefix}.vcf"
    bam_path = os.path.join(args.data_processed_dir, bam_output_name)
    vcf_path = os.path.join(args.data_processed_dir, vcf_output_name)

    print(f"Starting {mode} alignment for {args.fastq_r1_name}" + (f" and {args.fastq_r2_name}" if fastq_r2_path else "") + "...")

    # STAGE 1 & 2: Alignment, Filtering, Sorting (Piping)
    
    # 1. minimap2 (Alignment)
    minimap2_cmd = [
        'minimap2', 
        '-a',          # Output in SAM format
        '-x', 'sr',    # Preset for short reads
        fasta_path,    # Reference FASTA
        fastq_r1_path, # Query FASTQ R1
    ]
    # Add R2 if present**
    if fastq_r2_path:
        minimap2_cmd.append(fastq_r2_path) # minimap2 expects R2 path after R1 path
        
    p1_minimap2 = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # 2. samtools view (Filter)
    samtools_view_cmd = ['samtools', 'view', '-h', '-F', '0x900', '-']
    p2_samtools_view = subprocess.Popen(samtools_view_cmd, stdin=p1_minimap2.stdout, 
                                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p1_minimap2.stdout.close() 

    # 3. samtools sort (Sort and Write to file)
    samtools_sort_cmd = ['samtools', 'sort', '-O', 'bam', '-', '-o', bam_path]
    p3_samtools_sort = subprocess.Popen(samtools_sort_cmd, stdin=p2_samtools_view.stdout, 
                                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p2_samtools_view.stdout.close() 

    _, sort_err = p3_samtools_sort.communicate() 
    if p3_samtools_sort.returncode != 0:
        print(f"Error during alignment/sort:\n{sort_err.decode()}")
        return

    print(f"Alignment and sorting complete. BAM saved to {bam_path}")

    # 4. BAM Indexing (Generating .bai file)
    print(f"Indexing BAM file to create {bam_output_name}.bai...")
    samtools_index_cmd = ['samtools', 'index', bam_path]
    
    try:
        subprocess.run(samtools_index_cmd, check=True, capture_output=True, text=True)
        print("BAM indexing (.bai file) complete.")
    except subprocess.CalledProcessError as e:
        print(f"Error during BAM indexing:\n{e.stderr}")
        return

    # 5. Variant Calling (Piping)
    # Using BAM index created
    bcftools_mpileup_cmd = ['bcftools', 'mpileup', '-Ou', '-f', fasta_path, bam_path]
    p4_mpileup = subprocess.Popen(bcftools_mpileup_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    bcftools_call_cmd = ['bcftools', 'call', '-vc', '-Ov']
    
    with open(vcf_path, 'w') as vcf_file:
        p5_call = subprocess.Popen(bcftools_call_cmd, stdin=p4_mpileup.stdout, 
                                   stdout=vcf_file, stderr=subprocess.PIPE)
        p4_mpileup.stdout.close() 
        _, call_err = p5_call.communicate()
        if p5_call.returncode != 0:
            print(f"Error during variant calling:\n{call_err.decode()}")
            return
            
    print(f"Variant calling with bcftools complete. VCF saved to {vcf_path}")


def main():
    parser = argparse.ArgumentParser(description="Run minimap2, samtools, and bcftools for alignment and variant calling (supports single-end and paired-end).")
    parser.add_argument("fasta_name", help="Reference FASTA filename (e.g., EcoliK12-MG1655.fasta)")
    parser.add_argument("fastq_r1_name", help="FASTQ R1 (Forward) filename (e.g., NC_000913.fastq or NC_000913_R1.fastq)")
    # **Key Change:** R2 is now optional
    parser.add_argument("fastq_r2_name", nargs='?', default=None, help="Optional FASTQ R2 (Reverse) filename for paired-end data.")
    parser.add_argument("--data_raw_dir", default="data/raw", help="Directory containing the reference FASTA.")
    parser.add_argument("--data_processed_dir", default="data/processed", help="Directory for the FASTQ, BAM, and VCF files.")
    
    args = parser.parse_args()
    run_bcftools_alignment_and_variant_calling(args)

if __name__ == "__main__":
    main()