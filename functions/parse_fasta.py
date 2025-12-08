"""""
def parse_fasta(file_path):
    import os
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
    with open(file_path, "r") as file:
        lines = file.read().strip().split("\n")
        fasta_dictionary = {}
        current_id = None

    for line in lines:
        if line.startswith('>'): # Find identifier
            current_id = line[1:].strip()  # remove >
            fasta_dictionary[current_id] = ""
        else:
            fasta_dictionary[current_id] += line.strip()
            # append sequence lines
    print("Check the first 30 characters of sequence:")
    print(fasta_dictionary[current_id][:30])
    print("Check the FASTA id:")
    print(fasta_dictionary.keys())
    print("\n")
    return fasta_dictionary[current_id]
"""
#New version 16:50
def parse_fasta(file_path):
    import os
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
    with open(file_path, "r") as file:
        lines = file.read().strip().split("\n")
        fasta_dictionary = {}
        current_id = None

    for line in lines:
        if line.startswith('>'): # Find identifier
            current_id = line[1:].strip()  # remove >
            fasta_dictionary[current_id] = ""
        else:
            fasta_dictionary[current_id] += line.strip()
            # append sequence lines
            
    print("Check the first 30 characters of sequence:")
    print(fasta_dictionary[current_id][:30])
    print("Check the FASTA ids:")
    print(fasta_dictionary.keys())
    print("\n")
    return fasta_dictionary