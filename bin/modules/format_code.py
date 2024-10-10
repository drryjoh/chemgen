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

def write_formatted_code(code_directory, file_to_format, configuration = None, decorators = 'decorators', target_file = None, append = True):
    content = get_text_to_format(code_directory, file_to_format)
    
    if configuration == None:
        configuration = get_configuration("configuration.yaml", decorators=decorators)
    
    new_content = content.format(**vars(configuration))
    
    if target_file == None:
        target_file = file_to_format.replace('.in','')

    if append:
        with open(target_file, 'a') as file:
            file.write(new_content)
    else:
        with open(target_file, 'w') as file:
            file.write(new_content)      
