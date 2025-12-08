def make_mutations(sequence, no_snp, no_indel, max_indel_length):
    #no_snp is number of SNPs to mutate genome with
    #no_indel is number of indels to mutate genome with
    #sequence is sequence returned from parse_fasta

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
        elif insertion_not_deletion == 0: #"d" for deletion
            indel_dictionary[indel_id].extend(["d", left_pos_index, indel_length])
            # Send left most index those in following deletion to 
            # already_mutated to prevent it from getting mutated and lost
            for position in range(left_pos_index, left_pos_index + indel_length):
                already_mutated[position] = "d"
                deletion_positions.append(position) #Add deletions to list of deletion_positions
        

    # Track SNPs
    for snp_id in range(1, no_snp + 1):
        snp_index = random.randrange(genome_length) + 1
        while snp_index in already_mutated:
            snp_index = random.randrange(genome_length) + 1
        # Add SNPs to already_mutated
        already_mutated[snp_index] = "snp"
        # Add SNPs to snp_list
        snp_list.append(snp_index)


    # Make SNPs
    current_sequence = sequence
    sequence_list = list(current_sequence) #Mutable list form of sequence
    nucleotides = ["A", "G", "T", "C"]
    for snp_pos in snp_list:
        original_base = sequence_list[snp_pos - 1]
        possible_replacements = [base for base in nucleotides if base != original_base]
        new_base = random.choice(possible_replacements)
        sequence_list[snp_pos - 1] = new_base
        

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
    # Give preview of updated string
    print(updated_string[:30])
    print(indel_dictionary)
    print(already_mutated)
    print(snp_list)
    print(deletion_positions)
    return updated_string