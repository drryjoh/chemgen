
import yaml
import textwrap
from pathlib import Path
import subprocess
import sys
import os
import glob
from .format_code import *

def write_headers(file, headers):
    for header in headers:
        file.write(f"#include \"{header}\"\n")

def make_headers(code_directory, file_names, headers, configuration, destination_folder):
    for file in file_names:
        headers.append(file.replace('.in',''))
        write_formatted_code(code_directory, file, configuration, destination_folder)

def create_headers(configuration, chemistry_solver, destination_folder):
    code_directory = Path('src') / 'math_headers'
    file_names = ['multiply_divide.h.in', 'pow_gen.h.in','exp_gen.h.in', 'array_handling.h.in']
    headers = []
    make_headers(code_directory, file_names, headers, configuration, destination_folder)

    code_directory = Path('src') / 'thermophysics'
    file_names = ['constants.h.in',"thermally_perfect.h.in"]
    make_headers(code_directory, file_names, headers, configuration, destination_folder)

    code_directory = Path('src') / 'reaction_headers'
    file_names = ['arrhenius.h.in', 'third_body.h.in','falloff_troe.h.in', 'falloff_lindemann.h.in', 'falloff_sri.h.in', 'pressure_dependent_arrhenius.h.in']
    make_headers(code_directory, file_names, headers, configuration, destination_folder)

    if chemistry_solver:
        code_directory = Path('src') / 'solvers'
        if chemistry_solver == "rk4":
            file_names = ['rk4.h.in']
            make_headers(code_directory, file_names, headers, configuration, destination_folder)
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

