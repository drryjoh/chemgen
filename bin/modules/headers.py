
import yaml
import textwrap
from pathlib import Path
import subprocess
import sys
import os
import glob
from .format_code import *

def make_headers(code_directory, file_names, headers, configuration = None, decorators = 'decorators'):
    for file in file_names:
        headers.append(file.replace('.in',''))
        write_formatted_code(code_directory, file, configuration = configuration, decorators = decorators)

def create_headers(decorators = 'decorators', configuration = None):
    code_directory = Path('src') / 'math_headers'
    file_names = ['exp_gen.h.in', 'multiply_divide.h.in', 'pow_gen.h.in']
    headers = []
    make_headers(code_directory, file_names, headers, decorators = decorators, configuration = configuration)

    code_directory = Path('src') / 'thermophysics'
    file_names = ['constants.h.in']
    make_headers(code_directory, file_names, headers, decorators = decorators, configuration = configuration)

    code_directory = Path('src') / 'reaction_headers'
    file_names = ['arrhenius.h.in']
    make_headers(code_directory, file_names, headers, decorators = decorators, configuration = configuration)
    return headers

def clear_headers(directory):
    # Find all .cpp and .h files in the directory (and optionally subdirectories)
    h_files = glob.glob(os.path.join(directory, '*.h'))

    # Combine the lists of files
    files_to_delete = h_files

    # Delete each file
    for file_path in files_to_delete:
        try:
            os.remove(file_path)
            print(f"Deleted: {file_path}")
        except OSError as e:
            print(f"Error deleting {file_path}: {e}")

