def get_depth(mutated_sequence, depth_level, avg_fragment_length):
    import random
    random.seed(123)
    genome_length = len(mutated_sequence)
    # avg_fragment_length will be assumed to be 100, but can modify in function argument
    fragments_per_depth_level = genome_length // avg_fragment_length
    mutated_sequence_list = list(mutated_sequence)
    mutated_sequence_fragments = [] #List of all sequence fragments of length avg_fragment_length at all depth_level
    fragment_start_pos = []

    for _ in range(0, depth_level):
        for _ in range(fragments_per_depth_level):
            fragment_start_pos = random.randint(1, genome_length + 1 - avg_fragment_length + 1) #1-base and fencepost
            current_fragment = mutated_sequence_list[fragment_start_pos : fragment_start_pos + avg_fragment_length]
            complete_current_fragment = "".join(current_fragment)
            mutated_sequence_fragments.append(complete_current_fragment)


    # Report on first 4 fragments
    print(f"Getting fragments for {depth_level} depth:")
    print(f"At {depth_level} depth, there were {fragments_per_depth_level} fragments of length {avg_fragment_length}.")
    print("Examples of fragments as follows: ", mutated_sequence_fragments[0:3])
    #print(f"Check corresponding fragment start positions at: {fragment_start_pos[0]}, 
          #{fragment_start_pos[1]}, {fragment_start_pos[2]}, {fragment_start_pos[3]}
          #")
    print(f"Check corresponding fragment start positions at: {fragment_start_pos[0]}, {fragment_start_pos[1]}, {fragment_start_pos[2]}, {fragment_start_pos[3]}")
    return mutated_sequence_fragments