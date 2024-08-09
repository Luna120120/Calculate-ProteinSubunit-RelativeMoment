import re, sys, os 
import argparse
from Bio import Align
import pandas as pd
#---------------------------------------------------------------------------------------------------------
# Set up the program

welcome = '" Welcome to use the buried area calculation tool ! "'
print(f'\n{welcome: ^100}')
text_1 = " INFORMATION "
print(f'\n{text_1:-^100}')

print('''
Description:

This tool links the information of the UniPort sequence and the data of the respective PDB sequence''')

print('\n- Use -e to see example usage') 
print('- Use -h to see all arguments')

text_2 = " EXECUTION "
print(f'\n{text_2:-^100}\n')
#---------------------------------------------------------------------------------------------------------
# Defining executing functions
def print_example_usage():
    print('''
Example usage:

To run the program on the commend line, it requires:
   (1) This script
   (2) The .csv file contain the subtracted REL values (Δr) of each amino acid in the PDB sequence
   (3) a file contain the UniPort sequence, you can have a default file in name uniport_seq.txt
To process different UniPort sequences, you can keep refill new UniPort sequence into (3)
   
Input the script name, the filepath of the .csv file, -u and the filepath of the UniPort sequence
   e.g. python cal.py ./5w3f_A-AB.csv -u ./uniport_seq.txt ''')
    
    print('''
Function:

1. Input the data contain subtracted REL value of the respective PDB chain from file 
  (e.g. data for PDB model 5w3f, chain A, from file 5w3f_A-AB.csv),
   
   then create a pandas.DataFrame, translate the three_letter-code in to one-letter-code,
   and generate dictionary_1 contain key : value --- (pdb_index, pdb_aa) : rel_value
   e.g.
    (2, 'R') 0.0
    (3, 'E') 0.0
    ...
    (10, 'G') 0.0
    (11, 'Q') 33.2
    
2. Align the respect PDB sequence (query) to the UniPort sequence (target),
   e.g.
    target            0 MREVISINVGQAGCQIGNACWELYSLEHGIKPDGHLEDGLSKPKGGEEGFSTFFHETGYG
                      0 -|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    query             0 -REVISINVGQAGCQIGNACWELYSLEHGIKPDGHLEDGLSKPKGGEEGFSTFFHETGYG
                        ...
    target          420 EAREDLAALERDYIEVGADSYAEEEEF 447
                    420 ||||||||||||||||||||------- 447
    query           419 EAREDLAALERDYIEVGADS------- 439
    
   then use the alignment relationship to create dictionary_2, 
   which contain key : value --- (real_index, uniport_aa) : (aa_index, pdb_aa)
   e.g.
    (1, 'M') (1, '-')
    (2, 'R') (2, 'R')
    (3, 'E') (3, 'E')
    ...
    (446, 'E') (446, '-')
    (447, 'F') (447, '-')
    
3. Generate the strings of aligned target sequence and query sequence,
   e.g.
    aligned_target
    MREVISINVGQAGCQIGNACWELYSLEHGIKPDGHLEDGLSKPKGGEEGFSTFFHETGYG...EAREDLAALERDYIEVGADSYAEEEEF
    
    aligned_query
    -REVISINVGQAGCQIGNACWELYSLEHGIKPDGHLEDGLSKPKGGEEGFSTFFHETGYG...EAREDLAALERDYIEVGADS-------
    
   then use the (real_index, uniport_aa) from the aligned_target (the key of dictionary_2), 
   which is the uniport_key, to extract the value from dictionary_2, which is the pdb_key,
   then use the pdb_key to extract the value from dictionary_1, to get the corresponding REL value.
   
   Finally:
   (1) write everything into a file in name e.g. P09733_5w3f_cal.cav and export the file, in form
   e.g.
        Index UniProt_AA PDB_AA  REL
    0        1          M      -    -
    1        2          R      R  0.0
    2        3          E      E  0.0
    ..     ...        ...    ...  ...
    445    446          E      -    -
    446    447          F      -    -
    
   (2) calculate the Moment(M) of the Uniport sequence encoded protein
   by equation: 
                    M = sum of [r * (i - m)]
   where:
    i --- the index of each amino acid
    m --- the index of the middle amino acid, len(sequence)/2 if even, len(sequence)/2 + 1 if odd
    r --- the REL value (Δr), 
      get by subtracting the REL of aa in combined chains with the REL of aa in independent chain.
      ''')
    
# Read the table, extract the PDB name and values in the 'AA' column
def extract_pdb_name_and_generate_df(filepath):
    # Read the table from the file, assuming comma separation
    # Skip the first two rows which contain 'REM' and 'All-atoms'
    df = pd.read_csv(filepath, delimiter=',', skiprows=1)

    # Extract the PDB name from the file name
    basename = os.path.basename(filepath)
    pdb_name = basename.split('_')[0]
    rest_name = basename.split('_')[1]
    chain_id = rest_name.split('-')[0]
    
    # Create a new DataFrame that only retains information from column index 1, 3, and 5
    # which are the aa, index, REL
    selected_indices = [1, 3, 5]
    data_df = df.iloc[: ,selected_indices]
    # Swap the aa and the index columns
    cols = list(data_df.columns)
    cols[0], cols[1] = cols[1], cols[0]  # Swap the column names
    data_df = data_df[cols]  # Reindex the DataFrame with the new column order
    # Assign new column names to the DataFrame
    data_df.columns = ['Position', 'AA', 'REL']
    data_df = data_df[['Position', 'AA', 'REL']]
    
    return basename, pdb_name, chain_id, data_df


# Translate three letter aa code into one letter aa code
def convert_three_letter_code_to_one_letter_code(three_letter_list):
    aa_dict = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
        'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        # Additional potential codes
        'SEC': 'U', 'PYL': 'O', 'ASX': 'B', 'GLX': 'Z',
        'XLE': 'J', 'XAA': 'X', 'TER': '*', 'UNK': 'X'}
    one_letter_list = []
    for aa in three_letter_list:
        one_letter_code = aa_dict.get(aa, 'X')  # Default to 'X' if aa is not found
        one_letter_list.append(one_letter_code)
        
    return one_letter_list


# Change the code from three-letters to one-letter, generate a sequence string
def process_df_and_generate_sequence_string(data_df):
    # Convert the second column (amino acids) to a list of three-letter codes
    three_letter_list = data_df.iloc[:, 1].tolist()
    # Convert the three-letter codes to one-letter codes using the function
    one_letter_list = convert_three_letter_code_to_one_letter_code(three_letter_list)
    # Replace the second column with the one-letter codes
    data_df.iloc[:, 1] = one_letter_list
    
    # Generate the PDB sequence 
    PDB_sequence = ''.join(one_letter_list)
    
    return PDB_sequence, data_df


def extract_uniprot_data(uniport_filepath):
    # Initialize an empty string to hold the sequence
    uniport_sequence = ''
    # Initialize a variable to hold the UniProt name
    uniport_name = ''
    
    # Compile the regular expressions for the header and sequence lines
    header_pattern = re.compile(r'^>sp\|(\w+)\|')
    sequence_pattern = re.compile(r'^[A-Z]+$')
    
    # Open the file for reading
    with open(uniport_filepath, 'r') as file:
        # Iterate over each line in the file
        for line in file:
            # Check if the line is a header (starts with '>')
            if line.startswith('>'):
                # Use the regular expression to extract the UniProt name
                header_match = header_pattern.match(line)
                if header_match:
                    uniport_name = header_match.group(1)
                else:
                    print('Invalid header format:', line)
                    return None, None
            else:
                # Check if the line contains only uppercase letters (valid sequence line)
                if sequence_pattern.match(line.strip()):
                    # Remove any whitespace and concatenate to the sequence
                    uniport_sequence += line.strip()
                else:
                    print('Invalid sequence format:', line)
                    return None, None
    
    return uniport_name, uniport_sequence


def perform_global_alignment(uniprot_sequence, PDB_sequence):
    # Initialize the aligner
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'  # Perform global alignment
    aligner.match_score = 1  # Score for identical characters
    aligner.mismatch_score = 0  # Score for non-identical characters
    aligner.open_gap_score = 0  # Score to open a gap
    aligner.extend_gap_score = 0  # Score to extend a gap# Define two protein sequences to be aligned
    # Set sequence
    target = uniprot_sequence
    query = PDB_sequence
    # Perform the alignment
    alignments = aligner.align(target, query)
    # Get the first alignment (the best alignment)
    alignment = alignments[0]
    
    return alignment

# Extract the (index and aa) of aligned target seq and aligned query seq
def extract_aligned_sequences_from_alignment(alignment):
    # Convert the alignment object to a string
    alignment_str = str(alignment)
    # Split the alignment into lines
    lines = alignment_str.strip().split("\n")
    
    # Initialize variables to hold the aligned sequences
    aligned_target = ''
    aligned_query = ''

    # Process the alignment block by block
    for i in range(0, len(lines), 4):   # Each block has 4 lines
        # Extract parts of the target and query sequences
        # by removes any leading or trailing whitespace with strip()
        # and then splits the line into parts based on whitespace with split()
        target_line_parts = lines[i].strip().split() 
        query_line_parts = lines[i+2].strip().split()
        
        # Check if the line starts with 'target' or 'query' and has at least 3 parts
        if target_line_parts[0].startswith('target') and len(target_line_parts) > 2:
            aligned_target += target_line_parts[2]  # The sequence part is the third element
        if query_line_parts[0].startswith('query') and len(query_line_parts) > 2:
            aligned_query += query_line_parts[2]  # The sequence part is the third element
    
    return aligned_target, aligned_query


def create_aa_rel_dict(data_df):
    # Initialize an empty dictionary
    aa_rel_dict = {}
    # Iterate over each row in the DataFrame
    for pdb_index, row in data_df.iterrows():
        # Add the key-value pair to the dictionary
        aa_rel_dict[(row['Position'], row['AA'])] = row['REL']
        
    return aa_rel_dict


def create_mapping_dict_from_alignment(aligned_target, aligned_query):
    mapping_dict = {}
    # Initialize indices for target and query, start counting from index = 1
    target_index, query_index = 1, 1

    # zip function iterates through pairs of amino acids from the aligned target and query sequences
    # Iterate over aligned sequences
    for target_aa, query_aa in zip(aligned_target, aligned_query):
        # Create a key-value pair for the current amino acids and their indices
        key = (target_index, target_aa)
        value = (query_index, query_aa)

        # Add the key-value pair to the mapping dictionary
        mapping_dict[key] = value

        # Always increment both indices, regardless of gaps
        target_index += 1
        query_index += 1

    return mapping_dict


def create_sequence_df(chain_id, uniprot_name, pdb_name, aligned_target, mapping_dict, aa_rel_dict):
    # Initialize a list to store the data for each row
    data_list = []
    
    # Iterate over the aligned target sequence
    for index, uniport_aa in enumerate(aligned_target):
        real_index = index+1
        # print(f'Index: {real_index}, UniProt_AA: {uniport_aa}') # This is the target key
        
        # Use the current target index and amino acid to get the corresponding query index and amino acid from mapping_dict
        query_key = mapping_dict.get((real_index, uniport_aa))
        # print(f'Query key: {query_key}')
        
        # If the query_key is found in position_aa_rel_dict, get the corresponding REL value
        # Otherwise, set the REL value to '-'
        rel_value = aa_rel_dict.get(query_key, '-') if query_key else '-'
        # print(f'REL value: {rel_value}')
        # If the query_key exists, get the corresponding PDB amino acid, otherwise use '-'
        pdb_aa = query_key[1] if query_key else '-'
        # print(f'pdb_aa: {pdb_aa}')
        
        # Create a dictionary for the current row with the amino acid and its index (starting from 1)
        row = {'Index': real_index, 'UniProt_AA': uniport_aa, 'PDB_AA': pdb_aa, 'REL': rel_value}
        # Add the row dictionary to the data list
        data_list.append(row)

    # Create a DataFrame from the list of row dictionaries
    sequence_df = pd.DataFrame(data_list)
    # Reorder the columns to have 'PDB_AA' before 'REL'
    sequence_df = sequence_df[['Index', 'UniProt_AA', 'PDB_AA', 'REL']]
    filename = f'{uniprot_name}_{pdb_name}_{chain_id}_cal.csv'
    with open(filename,'w') as file:
        sequence_df.to_csv(file, sep = ',', index=False, header=True, mode = 'a')
    
    return sequence_df


# Calculate the moment(M) of proteins (The truncated protein in the elute)
def calculate_moment(sequence_df):
    # Determine the middle position 'm' of the sequence
    n = len(sequence_df)
    m = (n // 2) if n % 2 == 0 else (n // 2) + 1
    
    # Initialize the sum
    total_sum = 0
    
    # Iterate through each row of the DataFrame
    for index, row in sequence_df.iterrows():
        # Check if the REL value is numeric
        if pd.notnull(row['REL']) and row['REL'] != '-':
            r = float(row['REL'])  # Convert REL to float
            i = row['Index']       # Position index
            # Calculate the product of r and (i - m) and add it to the sum
            total_sum += r * (i - m)
            total_sum = round(total_sum, 2)
            
    return total_sum

#---------------------------------------------------------------------------------------------------------
# Set the parser
parser = argparse.ArgumentParser(description='Process csv file. Align PDB protein sequence to the UniPort protein sequence. Calculate the Moment(protein)')
parser.add_argument('path', type=str,
                    help='Path to a .csv file or a folder of .csv files')
parser.add_argument('-u', '--uniport_sequence_path', type=str,
                    help='the UniPort sequence, which is the target sequence of alignment')
parser.add_argument('-e', '--example', action='store_true',
                    help='show example usage')
args = parser.parse_args()

# Define the main function
def main():
    if args.example:
        print_example_usage()
        processed = True
        return
    
    # Populate paths_to_process
    paths_to_process = []
    if os.path.isfile(args.path) and args.path.lower().endswith('.csv'):    # Check if the provided path is a .csv file
        paths_to_process.append(args.path)
    # elif os.path.isdir(args.path):                                          # Check if the provided path is indeed a directory
    #     for filename in os.listdir(args.path):                              # Iterate over files in the directory
    #         file_path = os.path.join(args.path, filename)                   # Construct the full file path
    #         if filename.lower().endswith('.csv'):
    #             paths_to_process.append(file_path)
    if not paths_to_process:
        print(f'No valid .csv files found in path: {args.path}. The program is quit, please try again.')
        return
    
    # Process each pathway
    for filepath in paths_to_process:
        # 1. Get PDB name (query name) and df contain relevant info
        basename, pdb_name, chain_id, data_df = extract_pdb_name_and_generate_df(filepath)
        print(f'Loaded relevant information of {pdb_name}, chain {chain_id}, from {basename}')
        print(data_df)
        # 2. Change code
        pdb_sequence, data_df = process_df_and_generate_sequence_string(data_df)
        print('\nChanged letter code')
        print(data_df)
        print(f'\nGot the respective sequence from the PDB model {pdb_name} as the target sequence')
        # 3. Get UniPort name (target name) and string
        uniprot_name, uniprot_sequence = extract_uniprot_data(args.uniport_sequence_path)
        print(f'Got the UniPort sequence {uniprot_name} as the target sequence')
        # 4. Perform a global alignment
        alignment = perform_global_alignment(uniprot_sequence, pdb_sequence)
        print(f'\nAlignment: \n{alignment}')
        # 5. Extract aligned sequences from the alignment
        aligned_target, aligned_query = extract_aligned_sequences_from_alignment(alignment)
        print(f'Got aligned target sequence')
        print(f'Got aligned query sequence')
        # 6. create aa_rel_dict
        aa_rel_dict = create_aa_rel_dict(data_df)
        # 7. create mapping_dict
        mapping_dict = create_mapping_dict_from_alignment(aligned_target, aligned_query)
        # 8. Write sequence_df 
        sequence_df = create_sequence_df(chain_id, uniprot_name, pdb_name, aligned_target, mapping_dict, aa_rel_dict)
        print(f'\nExported result in file {uniprot_name}_{pdb_name}_{chain_id}_cal.csv ')
        print(sequence_df)
        # 9. Calculate the moment(M) of proteins (The truncated protein in the elute)
        M = calculate_moment(sequence_df)
        print(f'The Moment(M) of PDB model {pdb_name}, chain {chain_id} encoded protein, with respect to the UniPort sequence {uniprot_name}: {M}')



if __name__ == '__main__':
    main()
else:
	print("running as module\n")