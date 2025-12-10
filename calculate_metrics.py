import argparse
import csv
import sys
import os

def parse_ground_truth(gt_filepath):
    """
    Parses the custom ground truth file into a set of unique variant keys.
    Key format: (POS, REF, ALT)
    """
    ground_truth = set()
    print(f"Parsing ground truth file: {gt_filepath}...")
    
    try:
        with open(gt_filepath, 'r') as f:
            # Use the csv reader with a tab delimiter
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                # Skip header/comment lines starting with '#'
                if not row or row[0].startswith('#'):
                    continue
                
                # Expected columns: CHROM, POS, REF, ALT, TYPE, LEN
                if len(row) >= 4:
                    # We use a 3-part key (POS, REF, ALT) for unique identification.
                    # Positions are 1-based.
                    pos = row[1]
                    ref = row[2].upper() # Convert to uppercase for robust comparison
                    alt = row[3].upper()
                    
                    # Store the unique identifier
                    ground_truth.add((pos, ref, alt))
    except FileNotFoundError:
        print(f"Error: Ground truth file not found at {gt_filepath}", file=sys.stderr)
        sys.exit(1)
        
    print(f"Loaded {len(ground_truth)} ground truth variants.")
    return ground_truth


def parse_vcf(vcf_filepath, ground_truth):
    """
    Parses the VCF file, compares variants against the ground truth, 
    and calculates True Positives (TP) and False Positives (FP).
    
    The function also removes True Positives from the 'ground_truth' set, 
    which allows the remaining size of the set to equal the False Negatives (FN).
    """
    true_positives = 0
    false_positives = 0
    
    print(f"Comparing VCF file: {vcf_filepath} against ground truth...")
    
    try:
        with open(vcf_filepath, 'r') as f:
            for line in f:
                # Skip header/meta information lines
                if line.startswith('#'):
                    continue
                
                # VCF standard columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE
                parts = line.strip().split('\t')
                
                if len(parts) >= 5:
                    vcf_pos = parts[1]
                    vcf_ref = parts[3].upper()
                    
                    # VCF ALT field can contain multiple alternates separated by comma
                    vcf_alts = parts[4].upper().split(',')
                    
                    # Check each ALT allele reported by the VCF caller
                    for vcf_alt in vcf_alts:
                        vcf_key = (vcf_pos, vcf_ref, vcf_alt)
                        
                        if vcf_key in ground_truth:
                            # True Positive (TP): Caller found a real variant
                            true_positives += 1
                            
                            # Remove the key to calculate False Negatives later
                            ground_truth.remove(vcf_key) 
                        else:
                            # False Positive (FP): Caller reported a variant that is not in GT
                            false_positives += 1
                            
    except FileNotFoundError:
        print(f"Error: VCF file not found at {vcf_filepath}", file=sys.stderr)
        sys.exit(1)

    return true_positives, false_positives


def calculate_metrics(tp, fp, fn):
    """Calculates Recall, Precision, and F1 score."""
    
    # Total Variants in Ground Truth (all actual positive cases)
    total_ground_truth = tp + fn
    
    # Total Variants Called by VCF (all predicted positive cases)
    total_called = tp + fp
    
    # Calculate Recall (Sensitivity)
    # Recall = TP / (TP + FN)
    if total_ground_truth > 0:
        recall = tp / total_ground_truth
    else:
        recall = 0.0

    # Calculate Precision (PPV)
    # Precision = TP / (TP + FP)
    if total_called > 0:
        precision = tp / total_called
    else:
        precision = 0.0
        
    # Calculate F1 Score (Harmonic mean of Precision and Recall)
    if (precision + recall) > 0:
        f1_score = 2 * (precision * recall) / (precision + recall)
    else:
        f1_score = 0.0

    return recall, precision, f1_score


def main():
    parser = argparse.ArgumentParser(description="Calculate Recall and Precision by comparing a simulated ground truth file against a VCF file.")
    parser.add_argument("gt_file", help="Path to the simulated ground truth file (e.g., ground_truth_variants.txt).")
    parser.add_argument("vcf_file", help="Path to the VCF file output by the variant caller.")
    
    args = parser.parse_args()

    # 1. Parse Ground Truth (Total Actual Variants)
    ground_truth = parse_ground_truth(args.gt_file)
    initial_total_variants = len(ground_truth)
    
    # 2. Parse VCF and Determine TP and FP
    true_positives, false_positives = parse_vcf(args.vcf_file, ground_truth)
    
    # 3. Determine False Negatives (FN)
    # The remaining variants in the 'ground_truth' set are those that were missed (FN)
    false_negatives = len(ground_truth)
    
    # 4. Calculate Final Metrics
    recall, precision, f1_score = calculate_metrics(true_positives, false_positives, false_negatives)
    
    print("\n" + "="*40)
    print("      Variant Caller Performance Metrics")
    print("="*40)
    
    # Output Confusion Matrix components
    print("\nConfusion Matrix Components:")
    print(f"  True Positives (TP):   {true_positives}")
    print(f"  False Positives (FP):  {false_positives}")
    print(f"  False Negatives (FN):  {false_negatives}")
    
    # Output Total Counts
    print("\nTotal Counts:")
    print(f"  Total GT Variants:     {initial_total_variants} (TP + FN)")
    print(f"  Total Called Variants: {true_positives + false_positives} (TP + FP)")
    
    # Output Final Scores
    print("\nFinal Scores:")
    print(f"  Recall (Sensitivity):    {recall:.4f}")
    print(f"  Precision (PPV):         {precision:.4f}")
    print(f"  F1 Score:                {f1_score:.4f}")
    print("="*40)


if __name__ == "__main__":
    main()