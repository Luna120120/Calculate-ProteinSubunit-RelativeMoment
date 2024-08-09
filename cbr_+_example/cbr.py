import re, sys, os 
import argparse
from itertools import combinations
import subprocess
import pandas as pd

#---------------------------------------------------------------------------------------------------------
# Set up the program

welcome = '" Welcome to use the buried area calculation tool ! "'
print(f'\n{welcome: ^100}')
text_1 = " INFORMATION "
print(f'\n{text_1:-^100}')

print('\nThis tool has 2 functions:')
print('''
1. Generate .pdb files for individual chain and combined chains, for:
   (1) all
   (2) a part of
   (3) 2 designated chains
For one or a folder of PDB models.''')
print('''
2. (1) Run an external program "naccess" that calculate the absolute (ABS) and relative (REL) 
       accessibility of each residues with in the amino acid sequence, in both (a) individual chain 
       and (2) complex status, and save the results in (individuals and combined) .rsa files.
   (2) Subtract the data in the [combined.rsa file] by the respective data in [individual.rsa files] 
       of its consist two chains to calculate the buried residues in complex, save in .csv files.
   (3) --- After subtraction, positive values in data indicates buried residues in the complex
       If there is any buried residues in each subunit chain, write the result in .csv files.''')
print('''
Note: 
1. The PDB model name can not contain ".", "-", or "_", otherwise would cause error.
   The conventional 4 letter + number name is recommended, e.g. 5w3f
2. You have to make sure the "naccess" program export .rsa files to the right directory that this 
   script can read in, otherwise would cause error. By default, this script read in .rsa files from 
   the same directory it is located in. You can change the output pathway in the "naccess" scrip.
      ''') 
print('\n- Use -e to see example usage') 
print('- Use -h to see all arguments')

text_2 = " EXECUTION "
print(f'\n{text_2:-^100}\n') 

#---------------------------------------------------------------------------------------------------------
# Defining executing functions

# Part 0: Read file and extract basic information: (1)pdb_id, (2)chain numbers, (3)chain IDs, (4)data for each chain

# Printing example usage
def print_example_usage():
    print('\nExample usage:\n')
    print('''
This tool has 3 modes:

1. Mode [-a --all]: Process all chains in each model
   Input example: "python ./cbr.py 2uv8.pdb -a"

2. Mode [-p --part]: Process a part of chains in each model. 
   This mode run similar to the -a --all mode, but only process data for partial chains.
   It requires to specify the chains IDs after -p flag (case sensitive).
   e.g. PDB model 6oa9 has 10 chains ['A','B','C','D','E','F','G','H','I','J']
        Input example: "python ./cbr.py 2uv8.pdb -p A B E F"
        This mode will only perform calculation on the unique two-chian-complex 
        combinations between these four chains.

3. Mode [-d --designated]: Process two designated chains in each model.
   This mode only run on designated two chains of each PDB model.
   e.g. PDB model 6oa9 has 10 chains ['A','B','C','D','E','F','G','H','I','J']
        Input example: "python ./cbr.py 2uv8.pdb -d"
        The IDs of chains in this .pdb file will be print to screen, then an input request 
        will popup. Enter e.g. A B (case sensitive) in input popup, the program will only 
        perform calculation on the unique combination of these 2 chains.


This tool perform calculations in two stages:

Stage 1: 
   Generates .pdb files for each individual chain 
   and all possible unique two-chain-complex combinations.
   e.g. for PDB model 4ffb has 3 chains ['A','B','C'], generate:
        "4ffb_A.pdb", "4ffb_B.pdb", "4ffb_C.pdb",
        "4ffb_AB.pdb", "4ffb_AC.pdb", "4ffb_BC.pdb"
        
Stage 2: --- will only perform if select "y" in the popup: "Type "y" to run or "n" to quit:"

   (1) Calculate the absolute and relative accessibility of residues by the .pdb files, get:
       "4ffb_A.rsa", "4ffb_B.rsa", "4ffb_C.rsa",
       "4ffb_AB.rsa", "4ffb_AC.rsa", "4ffb_BC.rsa"
       
   (2) Calculate the subtracted data of each unique two-chain-complex by the data of 
       its consisted two individual chains, get:
       e.g. "4ffb_A-AB.csv", "4ffb_B-AB.csv",
            "4ffb_A-AC.csv", "4ffb_C-AC.csv",
            "4ffb_B-BC.csv", "4ffb_C-BC.csv"
            
   (3) If there is any, write the data of buried residues, get:
       e.g. "4ffb_A-AB_buried_residues.csv", "4ffb_B-AB_buried_residues.csv",
            "4ffb_A-AC_buried_residues.csv", "4ffb_C-AC_buried_residues.csv",
            "4ffb_B-BC_buried_residues.csv", "4ffb_C-BC_buried_residues.csv"
''') 

# Populate paths_to_process_list
def populate_paths_to_list(path):
    paths_to_process_list = []

    if os.path.isfile(path) and path.lower().endswith('.pdb'):    # Check if the provided path is a .pdb file
        paths_to_process_list.append(path)
        return paths_to_process_list
    elif os.path.isdir(path):                              # Check if the provided path is indeed a directory
        for filename in os.listdir(path):                  # Iterate over files in the directory
            filepath = os.path.join(path, filename)       # Construct the full file path
            if filename.lower().endswith('.pdb'):
                paths_to_process_list.append(filepath)
        return paths_to_process_list
    
    if not paths_to_process_list:
        print(f'No valid PDB files found in path: {path}. The program is quit, please try again.')
        return

# Extract the model id from the .pdb file
def extract_pdb_id(filepath):
    id, _ = os.path.splitext(os.path.basename(filepath)) # Get the basename form the file path, split the name a the extension, just get the name
    return id

# Extract the chains id and information lines of chains from the .pdb file
def read_file_and_extract_chain_id(pdb_file):
    with open(pdb_file, 'r') as file:
        lines = file.readlines()
        
    chain_set = set()
    for line in lines:
        if line.startswith('ATOM'):
            chain_id = line[21]
            chain_set.add(chain_id)
    # print(f'The type of chains: {type(chains)}')
    # print(f'The type of sorted chains: {type(sorted (chains))}')
    return chain_set, lines

# Part 1: Generate all individual and combined files for each PDB model

# Write info of each chain into individual files
def write_individual_files(pdb_id, chain_set, lines):
    filenames = []
    for chain in chain_set:
        filename = f'{pdb_id}_{chain}.pdb'
        with open(filename, 'w') as file:
            for line in lines:
                if line.startswith('ATOM') and line[21] == chain:
                    file.write(line)
        print(f'Written individual-chain-file: {filename}')
        filenames.append(filename)
    return filenames

# For a part of, or all chains, write info of every pair of chains into combined files
def write_part_or_all_combined_files(pdb_id, chain_set, lines):
    pdb_filename_list = []
    # Get all unique combinations of two chains
    chain_combinations = combinations(sorted(chain_set), 2)  
    for chain_pair in chain_combinations:
        # Skip if the chains are the same or if the reverse pair has already been processed
        if chain_pair[0] == chain_pair[1] or (chain_pair[1], chain_pair[0]) in pdb_filename_list:
            continue
        pdb_filename = f'{pdb_id}_{"".join(chain_pair)}.pdb'
        with open(pdb_filename, 'w') as file:
            for line in lines:
                if line.startswith('ATOM') and line[21] in chain_pair:
                    file.write(line)
        print(f'Written combined-chain-file: {pdb_filename}')
        pdb_filename_list.append(pdb_filename)
    return pdb_filename_list

# For two designated chains, write info of this pair into a combined file
def write_designated_combined_files(pdb_id, chain_set, lines):
    pdb_filename_list = []
    if len(chain_set) == 2 and chain_set[0] != chain_set[1]:
        pdb_filename = f'{pdb_id}_{chain_set[0]}{chain_set[1]}.pdb'
        with open(pdb_filename, 'w') as file:
            for line in lines:
                if line.startswith('ATOM') and (line[21] == chain_set[0] or line[21] == chain_set[1]):
                    file.write(line)
        print(f"Written combined-chains-file: {pdb_filename}")
        pdb_filename_list.append(pdb_filename)
    else:
        print("Error: Please specify two different chains.")
    return pdb_filename_list

# Part 2

# 2.1 Perform naccess on individual-chain-files and combined-chains-files
def perform_naccess_on_files(pdb_filename):
    list_rsa_filename = []
    basename = pdb_filename.split('.')[0]
    rsa_filename = f'{basename}.rsa'
    try:
        subprocess.run(["./naccess", pdb_filename], check=True)
        list_rsa_filename.append(rsa_filename)
        print(f'"naccess" successfully processed {pdb_filename}')
    except subprocess.CalledProcessError as e:
        print(f'Error running "naccess" on {pdb_filename}: {e}')
    except FileNotFoundError as e:
        print(f'The "naccess" program was not found: {e}')
    return list_rsa_filename

# 2.2 Extract data from .rsa file and return as a pandas DataFrame
def extract_rsa_data_into_df(rsa_filename):
    # Define a regular expression pattern to match the lines starting with 'RES'
    # and capture the different parts of the line.
    line_pattern = re.compile(r'^RES\s+(\w{3})\s+([A-Za-z])\s*(\d+)\s+(-?[\d\.]+)\s+(-?[\d\.]+)\s+(-?[\d\.]+)\s+(-?[\d\.]+)\s+(-?[\d\.]+)\s+(-?[\d\.]+)\s+(-?[\d\.]+)\s+(-?[\d\.]+)\s+(-?[\d\.]+)\s+(-?[\d\.]+)')
    
    with open(rsa_filename, 'r') as file:
        data = []
        for line in file:
            match = line_pattern.match(line)
            if match:
                # If a match is found, extract all the captured groups
                row = match.groups()
                # Prepend 'RES' to the row tuple
                new_row = ('RES',) + row
                data.append(new_row)
    
    # Create a DataFrame with the correct column names
    df = pd.DataFrame(data, columns=None)
    
    # Extract the base name of the file without the extension
    df_name = rsa_filename.split('.')[0]
    # print(f'The DataFrame {df_name} is created, shape: {df.shape}')
    # print(f'{df}')
    return df, df_name

# 2.3 Match names
def find_matching_individual_dfs_to_combined_df(df_combined_name):
    """
    Given the name of a combined DataFrame, extract and return the names of the 
    individual DataFrames that correspond to the chains in the combined DataFrame.
    """
    # Construct the names of the corresponding individual DataFrames
    pdb_id, chain_set = df_combined_name.split('_')
    chain1_id, chain2_id = chain_set[0], chain_set[1]
    
    df_individual_1_name = f"{pdb_id}_{chain1_id}"
    df_individual_2_name = f"{pdb_id}_{chain2_id}"
    
    return df_individual_1_name, df_individual_2_name

# Sort decimal places
def format_decimal_places(df):
    # Identify ABS and REL columns by their position (assuming ABS and REL alternate after the first 2 columns)
    abs_columns = df.columns[4::2]  # Starting from the 5rd column, every second column is 'ABS'
    rel_columns = df.columns[5::2]  # Starting from the 6th column, every second column is 'REL'
    # Format ABS columns with 2 decimal places and REL columns with 1 decimal place
    for col in abs_columns:
        df.loc[:, col] = pd.to_numeric(df[col], errors='coerce')
        df.loc[:, col] = df[col].apply(lambda x: f'{x:.2f}' if pd.notnull(x) else x)
    for col in rel_columns:
        df.loc[:, col] = pd.to_numeric(df[col], errors='coerce')
        df.loc[:, col] = df[col].apply(lambda x: f'{x:.1f}' if pd.notnull(x) else x)
    return df

# 2.4 Perform the subtraction of ABS and REL values by the individual DataFrame to the combined DataFrame(each chain part)
#     write the results to file
def subtract_rsa_values_and_write_to_file(df_chain_1, df_chain_2, df_combined, df_combined_name):
    # Play with names
    basename = df_combined_name.split('.')[0]
    pdb_id, chain_set = basename.split('_')
    chain1_id, chain2_id = chain_set[0], chain_set[1]
    
    # Copy DataFrames to avoid altering original data and convert all values to float
    df1 = df_chain_1.copy()
    df2 = df_chain_2.copy()
    df_combined_1 = df_combined[df_combined.iloc[:, 2] == chain1_id].copy().reset_index(drop=True)  # Filter for chain1
    df_combined_2 = df_combined[df_combined.iloc[:, 2] == chain2_id].copy().reset_index(drop=True)  # Filter for chain2
    
    # print(f'The df_combined_1 is prepared, shape: {df_combined_1.shape} \n{df_combined_1}')
    # print(f'The df_combined_2 is prepared, shape: {df_combined_2.shape} \n{df_combined_2}')
    # # Define file names for the output
    # df_combined_1_filename = f'{pdb_id}_{chain1_id}_from_combined.csv'
    # df_combined_2_filename = f'{pdb_id}_{chain2_id}_from_combined.csv'
    # # Write df_combined_1 to a CSV file
    # df_combined_1.to_csv(df_combined_1_filename, index=False)
    # print(f'df_combined_1 written to file: {df_combined_1_filename}')
    # # Write df_combined_2 to a CSV file
    # df_combined_2.to_csv(df_combined_2_filename, index=False)
    # print(f'df_combined_2 written to file: {df_combined_2_filename}')
    
    # Ensure subtraction is only applied to value columns, i.e. excluding identifier columns (the first 4 columns)
    data1 = df1.columns[4:]
    data2 = df2.columns[4:]
    data_combined_1 = df_combined_1.columns[4:]
    data_combined_2 = df_combined_2.columns[4:]
    # Convert all data to numeric for subtraction
    for df in [df1, df2, df_combined_1, df_combined_2]:
        for col in df.columns[4:]:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Perform subtraction
    df_subtracted_1 = df1[data1].subtract(df_combined_1[data_combined_1], fill_value=0)
    df_subtracted_2 = df2[data2].subtract(df_combined_2[data_combined_2], fill_value=0)
    # print(f'\nSubtraction completed by {pdb_id}_{chain1_id} - {basename}({chain1_id} part)')
    # print(f'{df_subtracted_1}')
    # print(f'Subtraction completed by {pdb_id}_{chain2_id} - {basename}({chain2_id} part)')
    # print(f'{df_subtracted_2}')
    
    # Extract identifier columns from df1 and df2
    identifiers_1 = df1.iloc[:, :4].reset_index(drop=True)
    identifiers_2 = df2.iloc[:, :4].reset_index(drop=True)
    # Concatenate identifier columns with subtracted DataFrames
    df_subtracted_1 = pd.concat([identifiers_1, df_subtracted_1.reset_index(drop=True)], axis=1)
    df_subtracted_2 = pd.concat([identifiers_2, df_subtracted_2.reset_index(drop=True)], axis=1)
    # Consist the data decimal places (ABS:2, REL:1)
    df_subtracted_1 = format_decimal_places(df_subtracted_1)
    df_subtracted_2 = format_decimal_places(df_subtracted_2)
    # Reorder columns into numerical sequence by the NUM column
    df_subtracted_1 = df_subtracted_1.sort_index()
    df_subtracted_2 = df_subtracted_2.sort_index()
    
    # Write to .csv file
    df_subtracted_1_filename = f'{pdb_id}_{chain1_id}-{chain_set}.csv'
    df_subtracted_2_filename = f'{pdb_id}_{chain2_id}-{chain_set}.csv'
    # sentence = 'REM,File of subtracted absolute (ABS) and relative (REL) (%) accessibilities for\n'
    header_1 = 'REM,RES_NUM, , ,All-atoms, ,Total-Side, ,Main-Chain, ,Non-polar, ,All-polar\n'
    header_2 = 'REM,AA,CHAIN,NUM,ABS,REL,ABS,REL,ABS,REL,ABS,REL,ABS,REL\n'
    with open(df_subtracted_1_filename, 'w') as file:
        # file.write(sentence)
        file.write(header_1)
        file.write(header_2)
        df_subtracted_1.to_csv(file, sep=',', index=False, header=False, mode='a')  # Append mode in case we want to add more lines
        print(f'\nWritten subtraction file {df_subtracted_1_filename}, aa: {len(df_subtracted_1)}')
    with open(df_subtracted_2_filename, 'w') as file:
        # file.write(sentence)
        file.write(header_1)
        file.write(header_2)
        df_subtracted_2.to_csv(file, sep=',', index=False, header=False, mode='a')  # Append mode in case we want to add more lines
        print(f'Written subtraction file {df_subtracted_2_filename}, aa: {len(df_subtracted_2)}')
    return df_subtracted_1, df_subtracted_2, df_subtracted_1_filename, df_subtracted_2_filename

# 2.5 Check for buried residues and write to files
def check_for_buried_residues_and_write_to_file(df_subtracted_1, df_subtracted_2, df_subtracted_1_filename, df_subtracted_2_filename):
    # Play with names
    basename1 = df_subtracted_1_filename.split('.')[0]
    basename2 = df_subtracted_2_filename.split('.')[0]
    pdb_id, subtracted_name1 = basename1.split('_')
    pdb_id, subtracted_name2 = basename2.split('_')
    chain1_id, combined_id = subtracted_name1.split('-')
    chain2_id, combined_id = subtracted_name2.split('-')
    # Convert the relevant columns to numeric for comparison in the next step
    df_subtracted_1.iloc[:, 4] = pd.to_numeric(df_subtracted_1.iloc[:, 4], errors='coerce')
    df_subtracted_1.iloc[:, 5] = pd.to_numeric(df_subtracted_1.iloc[:, 5], errors='coerce')
    df_subtracted_2.iloc[:, 4] = pd.to_numeric(df_subtracted_2.iloc[:, 4], errors='coerce')
    df_subtracted_2.iloc[:, 5] = pd.to_numeric(df_subtracted_2.iloc[:, 5], errors='coerce')
    # Check: counted as a buried residue if the All-atoms-ABS or All-atoms-REL is greater than 0
    buried_residues_in_chain1 = df_subtracted_1.loc[(df_subtracted_1.iloc[:, 4] > 0) | (df_subtracted_1.iloc[:, 5] > 0)]
    buried_residues_in_chain2 = df_subtracted_2.loc[(df_subtracted_2.iloc[:, 4] > 0) | (df_subtracted_2.iloc[:, 5] > 0)]
    
    # Write to .csv files
    filename1 = f'{pdb_id}_{chain1_id}-{combined_id}_buried_residues.csv'
    filename2 = f'{pdb_id}_{chain2_id}-{combined_id}_buried_residues.csv'
    # sentence = 'REM,File of subtracted absolute (ABS) and relative (REL) (%) accessibilities for\n'
    header_1 = 'REM,RES_NUM, , ,All-atoms, ,Total-Side, ,Main-Chain, ,Non-polar, ,All-polar\n'
    header_2 = 'REM,AA,CHAIN,NUM,ABS,REL,ABS,REL,ABS,REL,ABS,REL,ABS,REL\n'
    if not buried_residues_in_chain1.empty:
        with open(filename1, 'w') as file:
            # file.write(sentence)
            file.write(header_1)
            file.write(header_2)
            buried_residues_in_chain1.to_csv(file, sep=',', index=False, header=False, mode='a')
        print(f'Written buried residues to file, PDB model {pdb_id} complex {combined_id}, by chain {chain1_id}, aa: {len(buried_residues_in_chain1)}')
        # print(f'{buried_residues_in_chain1}')
    else:
        print(f'No residues buried in PDB model {pdb_id} complex {combined_id}, by chain {chain1_id}')
    if not buried_residues_in_chain2.empty:
        with open(filename2, 'w') as file:
            # file.write(sentence)
            file.write(header_1)
            file.write(header_2)
            buried_residues_in_chain2.to_csv(file, sep=',', index=False, header=False, mode='a')
        print(f'Written buried residues to file, PDB model {pdb_id} complex {combined_id}, by chain {chain2_id}, aa: {len(buried_residues_in_chain2)}')
        # print(f'{buried_residues_in_chain2}')
    else:
        print(f'No residues buried in PDB model {pdb_id} complex {combined_id}, by chain {chain2_id}')
        
#---------------------------------------------------------------------------------------------------------
# Set the parser
parser = argparse.ArgumentParser(description='Process one or a folder of PDB files.')
parser.add_argument('path', type=str,
                    help='Path to a PDB file or a folder containing PDB files')
parser.add_argument('-a', '--all', action='store_true',
                    help='Process all chains in each model')
parser.add_argument('-p', '--part', nargs='+', 
                    help='Process a part of chains in each model, require at least two chain IDs, e.g. -p A B C')
parser.add_argument('-d', '--designated', action='store_true',
                    help='Process two designated chains in each model')
parser.add_argument('-e', '--example', action='store_true',
                    help='Show example usage')
args = parser.parse_args()

# Define the main function
def main():
    # Part 0: Read file and extract basic information: (1)pdb_id, (2)chain numbers, (3)chain IDs, (4)data for each chai
    path = args.path
    paths_to_process_list = populate_paths_to_list(path)
    for filepath in paths_to_process_list:
        pdb_id = extract_pdb_id(filepath)
        chain_set, lines = read_file_and_extract_chain_id(filepath)
        print(f'File of PDB model {pdb_id} is loaded for processing, contain {len(chain_set)} chains: {sorted(chain_set)}')
        # print(f'File type: pdb_id: {type(pdb_id)} pdb_chains: {type(chain_set)} pdb_lines: {type(lines)}')
    
    # Part 1: Generate all individual and combined files for each PDB model
    if args.example:
            print_example_usage()
            processed = True
            return
    else:
        print('\nPart 1:')
        
    individual_pdb_file_for_all_model_dic = {}
    combined_pdb_file_for_all_model_dic = {}
    
    individual_pdb_file_list = []
    combined_pdb_file_list = []
    
    processed = False
    
    print(f'\n----- Action on PDB model: {pdb_id} -----\n')
    
    if args.all:
        if len(chain_set) >= 2:
            # Generate .pdb individual-chain-files and combined-chains-files and save them into lists
            individual_pdb_file_list.extend(write_individual_files(pdb_id, chain_set, lines))
            combined_pdb_file_list.extend(write_part_or_all_combined_files(pdb_id, chain_set, lines))
            # print(f'individual_file_all_list: {individual_file_all_list}')
            # print(f'combined_file_all_list: {combined_file_all_list}')
            individual_pdb_file_for_all_model_dic[pdb_id] = individual_pdb_file_list
            combined_pdb_file_for_all_model_dic[pdb_id] = combined_pdb_file_list
            # print(f'dic_individual_pdb_file_for_all_model: {dic_individual_pdb_file_for_all_model}')
            # print(f'dic_combined_pdb_file_for_all_model: {dic_combined_pdb_file_for_all_model}')
            processed = True
        else:
            print(f'\The model {pdb_id} has less than two chains IDs. Skipping to the next model if possible.')
    
    elif args.part:
        input_chains = args.part
        if len(input_chains) >= 2:
            normalized_chains = [chain for chain in input_chains]
            if all(chain in chain_set for chain in normalized_chains):
                chain_data = [line for line in lines if line.startswith('ATOM') and line[21] in normalized_chains]
                if chain_data:  # Ensure there are lines for designated chains before proceeding
                    # Generate .pdb individual-chain-files and combined-chains-files and save them into lists
                    individual_pdb_file_list.extend(write_individual_files(pdb_id, normalized_chains, chain_data))
                    combined_pdb_file_list.extend(write_part_or_all_combined_files(pdb_id, normalized_chains, chain_data))
                    # print(f'individual_file_designated_list: {individual_file_designated_list}')
                    # print(f'combined_file_designated_list: {combined_file_designated_list}')
                    individual_pdb_file_for_all_model_dic[pdb_id] = individual_pdb_file_list
                    combined_pdb_file_for_all_model_dic[pdb_id] = combined_pdb_file_list
                    # print(f'individual_pdb_file_for_all_model_dic: {individual_pdb_file_for_all_model_dic}')
                    # print(f'combined_pdb_file_for_all_model_dic: {combined_pdb_file_for_all_model_dic}')
                    processed = True
                else:
                    print(f'\nThe .pdb file of model {pdb_id} has no valid data to process. Skipping to the next model if possible.')
            else:
                print('This model does not have all input chains. Skipping to the next model if possible.')
        else:
            print("Not enough valid chains specified. The program is quit, please check try again.")
            return  # Exit the program
    
    elif args.designated:
        print('Which two chains do you want to have individual and combined .pdb files generated?')
        user_input = input('Please input the chain IDs, e.g. A B: ')
        input_chains = user_input.split()
        # Validate the number of input chain IDs
        if len(input_chains) == 2:
            # Normalize chain IDs (if your IDs are case-sensitive, e.g. typed a b, adjust accordingly)
            normalized_chains = [chain for chain in input_chains]
            # Check if both designated chains are present in the model
            if all(chain in chain_set for chain in normalized_chains):
                chain_data = [line for line in lines if line.startswith('ATOM') and line[21] in normalized_chains]
                if chain_data:  # Ensure there are lines for designated chains before proceeding
                    # Generate .pdb individual-chain-files and combined-chains-files and save them into lists
                    individual_pdb_file_list.extend(write_individual_files(pdb_id, normalized_chains, chain_data))
                    combined_pdb_file_list.extend(write_designated_combined_files(pdb_id, normalized_chains, chain_data))
                    # print(f'individual_file_designated_list: {individual_file_designated_list}')
                    # print(f'combined_file_designated_list: {combined_file_designated_list}')
                    individual_pdb_file_for_all_model_dic[pdb_id] = individual_pdb_file_list
                    combined_pdb_file_for_all_model_dic[pdb_id] = combined_pdb_file_list
                    # print(f'individual_pdb_file_for_all_model_dic: {individual_pdb_file_for_all_model_dic}')
                    # print(f'combined_pdb_file_for_all_model_dic: {combined_pdb_file_for_all_model_dic}')
                    processed = True
                else:
                    print(f'\nThe .pdb file of model {pdb_id} has no valid data to process. Skipping to the next model if possible.')
            else:
                print('\nThis model does not have all input chains. The program is quit, please check and try again.')
                return  # Exit the program
        else:
            if len(input_chains) < 2:
                print('\nYou have to input two chain IDs. The program is quit, please check try again.')
                return  # Exit the program
            else:
                print('\nYou can only input two chain IDs. The program is quit, please check try again.')
                return  # Exit the program

    if not processed:
        print(f'\nError: Missing or inappropriate action argument. Please check your options.')
    
    # Part 2: (1) Run an external program "naccess" to calculate the absolute and relative accessability of each residue in individual chain and complex (combined chains)
    #         (2) Use the values in individual-chain-files to subtract the corresponding values in the combined-chains-files to calculate the buried area between two subunits if there is any
    if individual_pdb_file_list and combined_pdb_file_list:
        print('\nPart 2:\n')
        user_input = input('Do you want to check for the buried residues in two-subunit-complex? \nType "y" to run or "n" to quit: ').strip().lower()
    
        if user_input == 'y':
            
            for filepath in paths_to_process_list:
                pdb_id = extract_pdb_id(filepath)

                print(f'\n----- Action on PDB model: {pdb_id} -----\n')
            
                # 2.1 Generate .rsa files and save them into lists
                individual_pdb_file_list = individual_pdb_file_for_all_model_dic[pdb_id]
                combined_pdb_file_list = combined_pdb_file_for_all_model_dic[pdb_id]
                individual_rsa_file_list = []
                combined_rsa_file_list = []
                # For individual chain
                if individual_pdb_file_list:  # Check if the list is not empty
                    print(f'The individual .pdb files for [naccess] to process: {individual_pdb_file_list}')
                    for individual_pdb_file in individual_pdb_file_list:
                        individual_rsa_file_list.extend(perform_naccess_on_files(individual_pdb_file))
                    print(f'The individual .rsa files generated: {individual_rsa_file_list}\n')
                # For combined chains
                if combined_pdb_file_list:  # Check if the list is not empty
                    print(f'The combined .pdb files for [naccess] to process: {combined_pdb_file_list}')
                    for combined_pdb_file in combined_pdb_file_list:
                        combined_rsa_file_list.extend(perform_naccess_on_files(combined_pdb_file))
                    print(f'The combined .rsa files generated: {combined_rsa_file_list}')

                # 2.2 Extract data from .rsa files and return as pandas DataFrames
                df_individual_dic = {}
                df_combined_dic = {}    
                # For individual chain
                for individual_rsa_file in individual_rsa_file_list:
                    df_individual, df_individual_name = extract_rsa_data_into_df(individual_rsa_file)
                    df_individual_dic[df_individual_name] = df_individual
                # print(f'The dic_df_individual has: \n{dic_df_individual.keys()}\n')
                # For combined chains
                for combined_rsa_file in combined_rsa_file_list:
                    df_combined, df_combined_name = extract_rsa_data_into_df(combined_rsa_file)
                    df_combined_dic[df_combined_name] = df_combined
                # print(f'The dic_df_combined has: \n{dic_df_combined.keys()}\n')
                
                # 2.3 Find the matching individual DataFrames for each combined DataFrame, perform the subtraction
                for df_combined_name, df_combined in df_combined_dic.items():
                    df1_name, df2_name = find_matching_individual_dfs_to_combined_df(df_combined_name)
                    # Fetch the corresponding individual DataFrames from dic_df_individual
                    df1 = df_individual_dic.get(df1_name)
                    df2 = df_individual_dic.get(df2_name)
                    
                    # 2.4 P.4 Perform the subtraction and write to files
                    if df1 is not None and df2 is not None: # Ensure both individual DataFrames were found before proceeding
                        df_subtracted_1, df_subtracted_2, df_subtracted_1_filename, df_subtracted_2_filename = subtract_rsa_values_and_write_to_file(df1, df2, df_combined, df_combined_name)
                        
                        # 2.5 Check for buried residues and write to files
                        check_for_buried_residues_and_write_to_file(df_subtracted_1, df_subtracted_2, df_subtracted_1_filename, df_subtracted_2_filename)
                    else:
                        print('\nError: one or both individual DataFrames were not found')
        elif user_input == 'n':
            print("\nThe program is quit")
        else:
            print('\nInvalid input. The program is quit, please start again.')



if __name__ == '__main__':
    main()
else:
	print("running as module\n")
