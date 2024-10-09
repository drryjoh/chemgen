#!/usr/bin/env python3

import argparse
import sys
import cantera as ct
import numpy as np
import shutil
from pathlib import Path

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
    chemical_mechanism = Path(args.chemical_mechanism)
    destination_folder = Path(args.destination)
    
    # Check if the destination folder exists, if not, create it
    if not destination_folder.exists():
        print(f"Destination folder '{destination_folder}' does not exist. Creating it...")
        destination_folder.mkdir(parents=True, exist_ok=True)

    # Check if the chemical mechanism file exists before processing
    if not chemical_mechanism.is_file():
        print(f"Error: The file '{chemical_mechanism}' does not exist.")
        return

    # Core logic of the script
    print(f"Processing file: {chemical_mechanism}")
    
    # Call process_file with the correct paths
    process_file(chemical_mechanism, destination_folder)

def process_file(chemical_mechanism, destination_folder):
    """
    Function to process the file.
    """
    # Example file processing logic
    try:
        with chemical_mechanism.open('r') as file:
            content = file.read()
            print(f"File content:\n{content}")
    except FileNotFoundError:
        print(f"File {chemical_mechanism} not found.")
        sys.exit(1)
    
    # Remove the 'yaml' extension from the chemical mechanism file name
    chemical_mechanism_name = chemical_mechanism.stem

    destination_path = destination_folder / chemical_mechanism.name
    shutil.copy(chemical_mechanism, destination_path)

    print(f"File '{chemical_mechanism}' has been copied to '{destination_path}'.")

# Entry point
if __name__ == "__main__":
    main()

