def make_mutations(sequence, no_snp, no_indel, max_indel_length, output_gt_path):
    #no_snp is number of SNPs to mutate genome with
    #no_indel is number of indels to mutate genome with
    #sequence is sequence returned from parse_fasta
    print(f"Initialising mutations on sequence {sequence[:30]}..., with {no_snp} SNPs and {no_indel} indels of maximum {max_indel_length} bases long.")
    import random
    random.seed(123) #Maintain seed 123 throughout
    genome_length = len(sequence)
    if genome_length < (no_indel * max_indel_length + no_snp):
        error_message = "Sequence input error. Sequence too short for desired mutations."
        raise ValueError(error_message)
    print("Sequence length is sufficient. Proceeding with mutation generation...")

    indel_dictionary = {}
    # indel_dictionary structure is 
    # indel_id:[insertion/deletion, left_pos_index, indel_length]
    # Note indexes already mutated 
    already_mutated = {}
    snp_list = [] #Store all snp positions here
    deletion_positions = [] #All deletion positions 1-base
    ground_truth_variants = [] #List to store ground truth mutation details for precision and recall later

    # Track indels
    for indel_id in range(1, no_indel + 1):
        indel_dictionary[indel_id] = [] #Initialise dictionary with relevant indel id
        # Produce random indel positions
        # Insertion=1 or deletion=0
        insertion_not_deletion = random.randint(0, 1)

        # Find random index position
        left_pos_index = random.randrange(genome_length) + 1
        while left_pos_index in already_mutated:
            left_pos_index = random.randrange(genome_length) + 1
             
        # Find random index length
        indel_length = random.randint(1, max_indel_length + 1)

        # Document indel info into indel_dictionary
        # If insertion:
        if insertion_not_deletion: #"i" for insertion
            indel_dictionary[indel_id].extend(["i", left_pos_index, indel_length])
            # Send mutated left_pos as anchor to already_mutated
            #for index in range(left_pos_index, left_pos_index + indel_length):
            already_mutated[left_pos_index] = "i"

            # Record Insertion Details
            ground_truth_variants.append({
                'POS': left_pos_index,
                'REF': sequence[left_pos_index - 1],  # The base to be replaced/anchored
                'ALT': 'INSERT_LATER',                # Will be filled after sequence generation
                'TYPE': 'INS',
                'LEN': indel_length
            })

        elif insertion_not_deletion == 0: #"d" for deletion
            indel_dictionary[indel_id].extend(["d", left_pos_index, indel_length])
            # Send left most index those in following deletion to 
            # already_mutated to prevent it from getting mutated and lost
            for position in range(left_pos_index, left_pos_index + indel_length):
                already_mutated[position] = "d"
                deletion_positions.append(position) #Add deletions to list of deletion_positions

            # Record Deletion Details
            # REF is the deleted sequence, ALT is the base at the anchor position
            ref_seq = sequence[left_pos_index - 1 : left_pos_index - 1 + indel_length]
            
            ground_truth_variants.append({
                'POS': left_pos_index,
                'REF': ref_seq,
                'ALT': sequence[left_pos_index - 1], # Base that remains (becomes the anchor)
                'TYPE': 'DEL',
                'LEN': indel_length - 1 #Off by one error as initial base not counted in deletion length
            })


    # Track SNPs
    for snp_id in range(1, no_snp + 1):
        snp_index = random.randrange(genome_length) + 1
        while snp_index in already_mutated:
            snp_index = random.randrange(genome_length) + 1
        # Add SNPs to already_mutated
        already_mutated[snp_index] = "snp"
        # Add SNPs to snp_list
        snp_list.append(snp_index)

        # Record SNP Details
        ground_truth_variants.append({
            'POS': snp_index,
            'REF': sequence[snp_index - 1], # Original base
            'ALT': 'MUTATE_LATER',          # Will be filled after mutation
            'TYPE': 'SNP',
            'LEN': 1
        })


    # Make SNPs
    current_sequence = sequence
    sequence_list = list(current_sequence) #Mutable list form of sequence
    nucleotides = ["A", "G", "T", "C"]
    for snp_pos in snp_list:
        original_base = sequence_list[snp_pos - 1]
        possible_replacements = [base for base in nucleotides if base != original_base]
        new_base = random.choice(possible_replacements)
        sequence_list[snp_pos - 1] = new_base
        
        # Update the ground truth record:
        for variant in ground_truth_variants:
            if variant['POS'] == snp_pos and variant['TYPE'] == 'SNP':
                variant['ALT'] = new_base
                break


    # Make indels
    for indel in indel_dictionary:
        insertion = "i"
        deletion = "d"
        
        #Insertions
        if insertion in indel_dictionary[indel]:
            left_pos = indel_dictionary[indel][1]
            current_indel_length = indel_dictionary[indel][2]
            #indel_base_new = sequence_list[left_pos - 1]
            # Find left_pos base
            left_base = sequence_list[left_pos - 1]
            # Make random sequence of length current_indel_length
            insertion_bases = [left_base]
            insertion_bases.extend(random.choices(nucleotides, k=current_indel_length))
            # Then make list of bases into string with .join()
            inserted_bases = "".join(insertion_bases)
            # Update sequence at index of left_pos
            sequence_list[left_pos - 1] = inserted_bases     

            # Update the ground truth record for insertions:
            for variant in ground_truth_variants:
                if variant['POS'] == left_pos and variant['TYPE'] == 'INS':
                    # The ALT sequence is the new string that replaced the REF base
                    variant['ALT'] = inserted_bases
                    break

        else:
            print("Deletion detected or Indel key does not exist in indel_dictionary.")

        #Deletions
        if deletion in indel_dictionary[indel]:
            left_pos = indel_dictionary[indel][1] - 1 #0 base conversion (deletion includes left_pos base)
            current_indel_length = indel_dictionary[indel][2]
            # Retrieve left_pos base
            for position in deletion_positions:
                # Replace all deletion_positions in sequence_list with ""
                sequence_list[position - 1] = ""


    # Concatenate resulting sequence_list and return
    updated_string = "".join(sequence_list)

    # Write the ground truth file
    write_ground_truth(ground_truth_variants, "ground_truth_variants.txt")
    
    # Give preview of updated string
    print("Check here for string similarity in first 30 bases:")
    print(updated_string[:30])
    print(f"Indels: {indel_dictionary}")
    print(f"List of positions already mutated: \n{already_mutated}")
    print(f"SNP positions: \n{snp_list}")
    print(f"Deletion positions: \n{deletion_positions}")

    write_ground_truth(ground_truth_variants, output_gt_path) #Write {ground_truth}.txt file
    return updated_string



def write_ground_truth(variants, filename):
    
   #  Writes the list of ground truth variants to a tab-separated text file.
    
    # Define the header for the output file
    header = ["CHROM", "POS", "REF", "ALT", "TYPE", "LEN"]
    
    # Set a placeholder chromosome name (as FASTA often lacks one)
    chrom_name = "SIM_GENOME" 
    
    print(f"\nWriting ground truth to {filename}...")
    
    try:
        with open(filename, 'w') as f:
            # Write the header line and description
            f.write("##fileformat=Simulated_Ground_Truth_v1.0\n")
            f.write(f"##INFO=<ID=LEN,Number=1,Type=Integer,Description=\"Length of the indel\">\n")
            f.write(f"##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Type of variant (SNP, INS, DEL)\">\n")
            f.write("#" + "\t".join(header) + "\n")
            
            # Write each variant record
            for var in variants:
                # Format the output line:
                line = [
                    chrom_name, 
                    str(var['POS']), 
                    var['REF'], 
                    var['ALT'], 
                    var['TYPE'], 
                    str(var['LEN'])
                ]
                f.write("\t".join(line) + "\n")
        print(f"Successfully wrote {len(variants)} variants to {filename}.")
    except IOError as e:
        print(f"Error writing ground truth file: {e}")