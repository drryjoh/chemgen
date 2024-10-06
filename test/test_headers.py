#!/usr/bin/env python3
import yaml
import textwrap
from pathlib import Path

class Config:
    pass

def get_configuration(configuration_filename='configuration_header.yaml'):
    with open(configuration_filename, 'r') as file:
        config = yaml.safe_load(file)
    # Create an instance of Config and set attributes dynamically
    config_obj = Config()
    for key, value in config['decorators'].items():
        setattr(config_obj, key, value)

    return config_obj    

def get_text_to_format(code_directory, file_to_format):
    # Get the current directory
    current_dir = Path(__file__).resolve().parent
    file_path = current_dir.parent / code_directory / file_to_format

    with open(file_path, 'r') as file:
        content = file.read()
    content = content.strip()
    content = content.replace('\t', '    ')
    return content

def write_formatted_code(code_directory, file_to_format, target_file = None, append = True):
    content = get_text_to_format(code_directory, file_to_format)
    config = get_configuration("configuration_header.yaml")
    new_content = content.format(**vars(config))
    
    if target_file == None:
        target_file = file_to_format.replace('.in','')

    if append:
        with open(target_file, 'a') as file:
            file.write(new_content)
    else:
        with open(target_file, 'w') as file:
            file.write(new_content)      

def make_headers(code_directory, file_names, headers):
    for file in file_names:
        headers.append(file)
        write_formatted_code(code_directory, file)
    return headers

def create_headers():
    code_directory = Path('src') / 'reaction_headers'
    file_names = ['arrhenius.h.in']
    headers = []
    make_headers(code_directory, file_names, headers)

    code_directory = Path('src') / 'math_headers'
    file_names = ['exp_gen.h.in', 'multiply_divide.h.in', 'pow_gen.h.in']
    make_headers(code_directory, file_names, headers)

def main():
    create_headers()

if __name__ == "__main__":
    main()
