def create_fastq(file_path_and_name, mutated_fragments, quality_score_ascii):
    #for key in sequence_details:
    import os
    output_directory = os.path.dirname(file_path_and_name)
    if output_directory:
        os.makedirs(output_directory, exist_ok=True)
        
    max_read_id = len(mutated_fragments)

    with open(file_path_and_name, 'w') as f:
        for i in range(max_read_id):
            read_number = i + 1
            f.write(f"@{read_number}-->{fasta_details[0]}\n")
            #Sequence
            f.write(f"{mutated_fragments[read_number]}\n")
            #+ sign
            f.write("+\n")                
            #Quality Score
            #for base in range(len(mutated_fragments[read_id - 1])):
            quality_score_line = quality_score_ascii * len(mutated_fragments[read_number])
            f.write(f"{quality_score_line}\n")

            if len(mutated_fragments[i]) != len(quality_score_ascii[i]):
                print(f"Skipping read {read_number}: Sequence and quality lengths do not match.")
                continue
            
    
    return f"\nSuccessfully created FASTQ file in {file_path_and_name}"