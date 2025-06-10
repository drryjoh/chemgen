import importlib.util

from .configuration import *
from .headers import *
import cantera as ct
from .cmake import *

def load_custom_test(filepath):
    # Load a module from the given file path
    spec = importlib.util.spec_from_file_location("create_test", filepath)
    custom_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(custom_module)
    # Return the SourceWriter class from the custom module
    return custom_module.create_test

def run_command(command):
    """Run a shell command and check for errors."""
    try:
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(result.stdout.decode())  # Print standard output
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running: {command}")
        print(e.stderr.decode())  # Print error output
        sys.exit(1)  # Exit with error

def compile_cpp_code(build_dir, source_files, config):
    """Compile C++ code using g++ or other compilers."""
    os.makedirs(build_dir, exist_ok=True)
    os.makedirs(build_dir/"bin", exist_ok=True)
    build_option = ''
    if config['build'].get('chemgen_input'):
        build_option = config['build'].get('chemgen_input')
        print(f"building with opentions: {build_option}")
    # Command to compile C++ code
    compile_command = f"clang++ -std=c++17 {build_option} -O3 -o {build_dir}/bin/chemgen {' '.join(source_files)}"
    print(compile_command)
    print(f"Compiling C++ files: {source_files}")
    run_command(compile_command)

def run_tests(build_dir):
    """Run tests on the compiled binary."""
    test_command = f"./{build_dir}/bin/chemgen"
    print("Running tests...")
    run_command(test_command)

def compile(test_file, configuration_file, destination_folder, third_parties, compile = True):
    # Define directories and C++ source files
    build_directory = destination_folder.parent
    cpp_source_files = ['src'+'/'+test_file]
    # generate cmake
    
    generate_cmake_file(configuration_file, build_directory, third_parties)
    # Compile the C++ code
    if compile:
        compile_cpp_code(build_directory, cpp_source_files, configuration_file)
    # Run the tests
    run_tests(build_directory)
