def create_fastq(file_path_and_name, mutated_fragments, quality_score_ascii):
    #for key in sequence_details:
    with open(file_path_and_name, 'w') as f:
        for i in range(max_read_id):
            read_id = 1
            f.write(f"@{read_id[i]}-->fasta_details[0]\n")
            #Sequence
            f.write(f"{mutated_fragments[i]}\n")
            #+ sign
            f.write("+\n")                
            #Quality Score
            #for base in range(len(mutated_fragments[read_id - 1])):
            quality_score_line = quality_score_ascii * len(mutated_fragments[read_id - 1])
            f.write(f"{quality_score_line}\n")

            if len(mutated_fragments[i]) != len(quality_score_ascii[i]):
                print(f"Skipping read {read_id[i]}: Sequence and quality lengths do not match.")
                continue
    
    return f"\nSuccessfully created FASTQ file in {file_path_and_name}"