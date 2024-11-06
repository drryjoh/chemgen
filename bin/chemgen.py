#!/usr/bin/python3

import argparse
import sys
import cantera as ct
import numpy as np
import shutil
from pathlib import Path
from modules.process_files import *
from modules.configuration import *
from modules.headers import *
from modules.compile_and_run import *

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

# Define functions or classes here
def main():
    """
    Main function where the core logic of the script runs.
    """
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="A brief description of the script")
    parser.add_argument('chemical_mechanism', type=str, help="The name of the file")
    parser.add_argument('destination', type=str, help="The destination folder where the files will be generated")
    parser.add_argument("--custom-source", type=str, help="Path to custom source writer file")
    parser.add_argument("--custom-test", type=str, help="Path to custom test case writer file")
    parser.add_argument("--compile", action="store_true", help="Compile the source writer code")

    args = parser.parse_args()
    
    # Convert arguments to Path objects
    chemical_mechanism = find_chemical_mechanism(args.chemical_mechanism)
    destination_folder = Path(args.destination)/'src'
    
    # Check if the destination folder exists, if not, create it
    if not destination_folder.exists():
        print(f"Destination folder '{destination_folder}' does not exist. Creating it...")
        destination_folder.mkdir(parents=True, exist_ok=True)

    # Core logic of the script
    print(f"Processing file: {args.chemical_mechanism}")
    

    gas = ct.Solution(chemical_mechanism)
    [configuration, configuration_file] = get_configuration(configuration_filename='configuration.yaml')
    headers = process_cantera_file(gas, configuration, destination_folder,args)
    print("***headers***")

    if "types_inl.h" in headers:
        headers.remove("types_inl.h")
        headers.insert(0,"types_inl.h")
    

    test_file = ''
    
    if args.custom_test:
        try:
            test_file = 'chemgen.cpp'
            # Load the custom SourceWriter
            create_test = load_custom_test(args.custom_test)
            create_test(gas, args.chemical_mechanism, headers, test_file, configuration, destination_folder)

        except (FileNotFoundError, AttributeError) as e:
            print(f"Error loading custom test writer: {e}")
            sys.exit(1)

    else: #replace with run time argument
        test_file = 'chemgen.cpp'
        from modules.default_test import create_test
        create_test(gas, args.chemical_mechanism, headers, test_file, configuration, destination_folder)
    
    if args.compile:
        compile(test_file, configuration_file, destination_folder)

# Entry point
if __name__ == "__main__":
    main()

