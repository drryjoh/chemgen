#!/usr/bin/env python3
import yaml
import textwrap
from pathlib import Path
import subprocess
import sys
import os
import glob


class Config:
    pass

def get_configuration(configuration_filename='configuration_header.yaml'):
    with open(configuration_filename, 'r') as file:
        configuration = yaml.safe_load(file)
    # Create an instance of Config and set attributes dynamically
    config_obj = Config()
    for key, value in config[decorators].items():
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
    configuration = get_configuration("configuration_header.yaml")
    new_content = content.format(**vars(configuration))
    
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
        headers.append(file.replace('.in',''))
        write_formatted_code(code_directory, file, decorators)

def create_headers(decorators = 'decorators'):
    code_directory = Path('src') / 'math_headers'
    file_names = ['exp_gen.h.in', 'multiply_divide.h.in', 'pow_gen.h.in']
    headers = []
    make_headers(code_directory, file_names, headers, decorators)

    code_directory = Path('src') / 'thermophysics'
    file_names = ['constants.h.in']
    make_headers(code_directory, file_names, headers, decorators)

    code_directory = Path('src') / 'reaction_headers'
    file_names = ['arrhenius.h.in']
    make_headers(code_directory, file_names, headers, decorators)
    return headers

def run_command(command):
    """Run a shell command and check for errors."""
    try:
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(result.stdout.decode())  # Print standard output
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running: {command}")
        print(e.stderr.decode())  # Print error output
        sys.exit(1)  # Exit with error

def compile_cpp_code(build_dir, source_files):
    """Compile C++ code using g++ or other compilers."""
    os.makedirs(build_dir, exist_ok=True)
    
    # Command to compile C++ code
    compile_command = f"g++ -std=c++14 -o {build_dir}/output_program {' '.join(source_files)}"
    print(f"Compiling C++ files: {source_files}")
    run_command(compile_command)

def run_tests(build_dir):
    """Run tests on the compiled binary."""
    test_command = f"./{build_dir}/output_program"
    print("Running tests...")
    run_command(test_command)

def compile_header_test():
    # Define directories and C++ source files
    build_directory = "build"
    cpp_source_files = ["main.cpp"]

    # Step 1: Compile the C++ code
    compile_cpp_code(build_directory, cpp_source_files)

    # Step 2: Run the tests
    run_tests(build_directory)

def create_header_cpp(headers):
    configuration = get_configuration("configuration_header.yaml")
    with open('./main.cpp', 'w') as file: 
        file.write("#include <cmath>\n")
        for header in headers:
            file.write(f"#include \"{header}\"\n")
        content = """
#include <iostream>  // For printing the result to the console

int main() {{
    // Call the arrhenius function with the specified parameters
    {scalar} result = arrhenius({scalar_cast}(100),  {scalar_cast}(1.5), {scalar_cast}(1.3e6), {scalar_cast}(1800));
    {scalar} dresult_dtemperature = darrhenius_dtemperature({scalar_cast}(100), {scalar_cast}(1.5), {scalar_cast}(1.3e6),  {scalar_cast}(1800));

    // Output the result
    std::cout << "Result of arrhenius(100, 1.3e6, 1.5, 1800): " << result << std::endl;
    std::cout << "Result of darrhenius_dtemperature(100, 1.3e6, 1.5, 1800): " << dresult_dtemperature << std::endl;

    return 0;
}}
            """
        file.write(content.format(**vars(configuration)))

def clear_files(directory):
    # Find all .cpp and .h files in the directory (and optionally subdirectories)
    cpp_files = glob.glob(os.path.join(directory, '*.cpp'))
    h_files = glob.glob(os.path.join(directory, '*.h'))

    # Combine the lists of files
    files_to_delete = cpp_files + h_files

    # Delete each file
    for file_path in files_to_delete:
        try:
            os.remove(file_path)
            print(f"Deleted: {file_path}")
        except OSError as e:
            print(f"Error deleting {file_path}: {e}")

def main():
    decorators_to_test = ['decorators','decorators_float']
    for decorators in decorators_to_test:
        clear_files('./')
        headers = create_headers()
        create_header_cpp(headers)
        compile_header_test()

if __name__ == "__main__":
    main()
