import yaml
import textwrap
from pathlib import Path
import subprocess
import sys
import os
import glob

from .configuration import *

def get_text_to_format(code_directory, file_to_format):
    # Get the current directory
    current_dir = Path(__file__).resolve().parent
    file_path = current_dir.parent.parent / code_directory / file_to_format

    with open(file_path, 'r') as file:
        content = file.read()
    content = content.strip()
    content = content.replace('\t', '    ')
    return content

def write_formatted_code(code_directory, file_to_format, configuration, destination_folder, target_file = None, append = False):
    content = get_text_to_format(code_directory, file_to_format)
    try:
        new_content = content.format(**vars(configuration))
    except:
        # Split the text into lines for easier tracking
        lines = content.splitlines()
        for k, line in enumerate(lines):
            try:
                line.format(**vars(configuration))
            except:
                print(f"while formatting: {file_to_format}")
                print(f"unescaped braces found in line {k+1}")
                print(line)
        sys.exit(f"unescaped braces found in file {file_to_format}")

    
    if target_file == None:
        target_file = file_to_format.replace('.in','')

    if append:
        with open(destination_folder/target_file, 'a') as file:
            file.write(new_content)
    else:
        with open(destination_folder/target_file, 'w') as file:
            file.write(new_content)      
