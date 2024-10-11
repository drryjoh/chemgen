from .configuration import *
from .headers import *

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

def compile_header_test(test_file):
    # Define directories and C++ source files
    build_directory = "build"
    cpp_source_files = [test_file]

    # Step 1: Compile the C++ code
    compile_cpp_code(build_directory, cpp_source_files)

    # Step 2: Run the tests
    run_tests(build_directory)

def get_test_conditions():
    config_path = Path('test_configuration.yaml')
    
    if config_path.exists():
        with config_path.open('r') as file:
            configuration = yaml.safe_load(file)
    else:
        current_dir = Path(__file__).resolve().parent
        configuration_filename = current_dir.parent.parent/ 'test/test_configuration.yaml'
        print("**Test file test_configuration.yaml not found using /test/test_configuration.yaml **")
        with open(configuration_filename, 'r') as file:
            configuration = yaml.safe_load(file)

    # Extract values from the parsed YAML data
    temperature = configuration['test_conditions']['temperature']
    pressure = configuration['test_conditions']['pressure']
    species_list = configuration['test_conditions']['species']

    # Print out the parsed values
    print("Test configuration:")
    print(f"Temperature: {temperature}")
    print(f"Pressure: {pressure}")
    print("Species:")

    for species in species_list:
        print(f"  - Name: {species['name']}, MoleFraction: {species['MoleFraction']}")
    species_string  = ' '.join([f"{species['name']}:{species['MoleFraction']} " for species in species_list])
    return [temperature, pressure, species_string]

def create_test(gas, headers, test_file, configuration):
    with open(test_file, 'w') as file:
        file.write("#include <cmath>\n")
        file.write("#include <array>\n")
        file.write("#include <iostream>  // For printing the result to the console\n")
        write_headers(file, headers)
        [temperature, pressure, species_string] = get_test_conditions()
        gas.TPX = temperature, pressure, species_string
        concentrations = gas.concentrations
        concentration_test = '{species} species  = {{{array}}};'.format(array = ','.join(["{scalar_cast}({c})".format(c=c, **vars(configuration)) for c in concentrations]),**vars(configuration)) 
        
        content = """
// Overload << operator for std::array
template <typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& arr) {{
    os << "[ ";
    for (const auto& value : arr) 
    {{
        os << value << " ";
    }}
    os << "]";
    return os;
}}

int main() {{
    // Call the arrhenius function with the specified parameters
    {concentration_test}
    {scalar} temperature =  {temperature};
    {species} result = source(species, temperature);

    // Output the result
    std::cout << "Source test result:  " << result << std::endl;
    std::cout << "Cantera test result: " <<"{cantera_net_production_rates}"<<std::endl;

    return 0;
}}
            """
        file.write(content.format(**vars(configuration), concentration_test = concentration_test, temperature = temperature, cantera_net_production_rates = ' '.join([f"{npr}" for npr in gas.net_production_rates])))