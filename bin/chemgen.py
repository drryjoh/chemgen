#!/usr/bin/env python3

import argparse
import sys
import cantera as ct
import numpy as np
import shutil
from pathlib import Path
from modules.process_files import *
from modules.configuration import *
from modules.headers import *
from modules.tests import *

def find_chemical_mechanism(file_name):
    # Get the path to the current script (chemgen.py)
    script_path = Path(__file__).resolve()
    # Get the parent directory of the script (assumed to be /path_to_code/bin)
    base_directory = script_path.parent.parent
    # Define the path to the chemical_mechanisms directory
    chemical_mechanisms_dir = base_directory / 'chemical_mechanisms'

    # Check if the file has the .yaml extension; if not, append it
    if not file_name.endswith('.yaml'):
        file_name += '.yaml'

    # Convert the input file_name into a Path object
    chemical_mechanism = Path(file_name)

    # Case 1: Direct path provided by the user
    if chemical_mechanism.is_file():
        return chemical_mechanism

    # Case 2: The file is located in the chemical_mechanisms directory
    chemical_mechanism_in_dir = chemical_mechanisms_dir / chemical_mechanism.name
    if chemical_mechanism_in_dir.is_file():
        return chemical_mechanism_in_dir

    # If neither case works, print an error and return None
    print(f"Error: The file '{file_name}' does not exist.")
    return None

def generate_code(chemical_mechanism, destination_folder):
    """
    Function to process the file.
    """
    
    # Remove the 'yaml' extension from the chemical mechanism file name
    chemical_mechanism_name = chemical_mechanism.stem

    destination_path = destination_folder / chemical_mechanism.name
    shutil.copy(chemical_mechanism, destination_path)

    print(f"File '{chemical_mechanism}' has been copied to '{destination_path}'.")

# Define functions or classes here
def main():
    """
    Main function where the core logic of the script runs.
    """
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="A brief description of the script")
    parser.add_argument('chemical_mechanism', type=str, help="The name of the file")
    parser.add_argument('destination', type=str, help="The destination folder where the files will be generated")
    args = parser.parse_args()
    
    # Convert arguments to Path objects
    chemical_mechanism = find_chemical_mechanism(args.chemical_mechanism)
    destination_folder = Path(args.destination)
    
    # Check if the destination folder exists, if not, create it
    if not destination_folder.exists():
        print(f"Destination folder '{destination_folder}' does not exist. Creating it...")
        destination_folder.mkdir(parents=True, exist_ok=True)

    # Core logic of the script
    print(f"Processing file: {args.chemical_mechanism}")
    
    # Call process_file with the correct paths
    generate_code(chemical_mechanism, destination_folder)
    gas = ct.Solution(chemical_mechanism)
    configuration = get_configuration(configuration_filename='configuration.yaml')

    headers = process_cantera_file(gas, configuration)
    if True: #replace with run time argument
        test_file = 'chemgen_test.cpp'
        create_test(gas, headers, test_file, configuration, decorators = 'decorators')
        compile_header_test(test_file)

# Entry point
if __name__ == "__main__":
    main()

