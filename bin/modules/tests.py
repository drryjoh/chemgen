from .configuration import *
from .headers import *
import cantera as ct

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
    print(compile_command)
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

def get_test_conditions(chemical_mechanism):
    config_path = Path('test_configuration.yaml')
    if "simple_tb" in chemical_mechanism:
        print("**simple TB**")
        test_file = 'test_tb_configuration.yaml'
    else:
        test_file = 'test_configuration.yaml'
    config_path = Path(test_file)
    if config_path.exists():
        with config_path.open('r') as file:
            configuration = yaml.safe_load(file)
    else:
        current_dir = Path(__file__).resolve().parent
        configuration_filename = current_dir.parent.parent/ 'test/{0}'.format(test_file)
        print("**Test file test_configuration.yaml not found using /test/test_configuration.yaml **")
        with open(configuration_filename, 'r') as file:
            configuration = yaml.safe_load(file)
    print(config_path)
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

def create_test(gas, chemical_mechanism, headers, test_file, configuration):
    with open(test_file, 'w') as file:
        file.write("#include <cmath>\n")
        file.write("#include <array>\n")
        file.write("#include <iostream>  // For printing the result to the console\n")
        write_headers(file, headers)
        [temperature, pressure, species_string] = get_test_conditions(chemical_mechanism)
        gas.TPX = temperature, pressure, species_string
        concentrations = gas.concentrations
        concentration_test = '{species} species  = {{{array}}};'.format(array = ','.join(["{scalar_cast}({c})".format(c=c, **vars(configuration)) for c in concentrations]),**vars(configuration)) 
        enthalpies = gas.standard_enthalpies_RT * gas.T * ct.gas_constant/gas.molecular_weights
        entropies = gas.standard_entropies_R * ct.gas_constant/gas.molecular_weights
        energies = gas.standard_int_energies_RT * gas.T * ct.gas_constant/gas.molecular_weights
        gibbs = gas.standard_gibbs_RT * gas.T * ct.gas_constant
        equilibrium_constants = gas.equilibrium_constants
        gas.TPX = temperature, pressure, species_string

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
    {scalar} pressure_return = pressure(species, temperature);

    // Output the result
    std::cout << "Source test result:  " << result << std::endl;
    std::cout << "Cantera test result: " <<"{cantera_net_production_rates}"<<std::endl;

    std::cout << "Cantera species cps: " <<"{cantera_species_cp}"<<std::endl;
    std::cout << "Chemgen species cps: " << species_specific_heat_constant_pressure_mass_specific(temperature) <<std::endl;

    std::cout << "Cantera species enthalpies: " <<"{cantera_species_enthalpy}"<<std::endl;
    std::cout << "Chemgen species enthalpies: " << species_enthalpy_mass_specific(temperature) <<std::endl;

    std::cout << "Cantera species internal energies: " <<"{cantera_species_energies}"<<std::endl;
    std::cout << "Chemgen species internal energies: " << species_internal_energy_mass_specific(temperature) <<std::endl;

    std::cout << "Cantera species internal entropies: " <<"{cantera_species_entropy}"<<std::endl;
    std::cout << "Chemgen species internal entropies: " << species_entropy_mass_specific(temperature) <<std::endl;

    std::cout << "Cantera species gibbs energy: " <<"{cantera_species_gibbs}"<<std::endl;
    std::cout << "Chemgen species gibbs energy: " << species_gibbs_energy_mole_specific(temperature) <<std::endl;
    std::cout << "Chemgen species gibbs energy direct: " << species_enthalpy_mass_specific(temperature) - scale_gen(temperature, species_entropy_mass_specific(temperature)) <<std::endl;

    std::cout << "Cantera equilibrium constants: " <<"{equilibrium_constants}"<<std::endl;
    std::cout << "Chemgen equilibrium constants: " << equilibrium_constants(temperature) <<std::endl;
    

    std::cout << "Pressure: " <<pressure_return <<std::endl;
    std::cout << "Temperature Monomial at 300           : " <<temperature_monomial({scalar_cast}(300)) <<std::endl;
    std::cout << "Temperature Energy Monomial at 300           : " <<temperature_energy_monomial({scalar_cast}(300)) <<std::endl;
    std::cout << "Temperature Entropy Monomial at 300           : " <<temperature_entropy_monomial({scalar_cast}(300)) <<std::endl;
    std::cout << "Temperature Gibbs Monomial at 300           : " <<temperature_gibbs_monomial({scalar_cast}(300)) <<std::endl;
    

    std::cout << "Temperature Monomial Derivative at 300: " <<dtemperature_monomial_dtemperature({scalar_cast}(300)) <<std::endl;

    return 0;
}}
            """
        file.write(content.format(**vars(configuration), 
        concentration_test = concentration_test, 
        temperature = temperature, 
        cantera_net_production_rates = ' '.join([f"{npr}" for npr in gas.net_production_rates]),
        cantera_species_cp = ' '.join([f"{scp}" for scp in gas.standard_cp_R * ct.gas_constant / gas.molecular_weights]),
        cantera_species_enthalpy = ' '.join([f"{enth}" for enth in enthalpies]),
        cantera_species_entropy = ' '.join([f"{ent}" for ent in entropies]),
        cantera_species_energies = ' '.join([f"{energy}" for energy in energies]),
        cantera_species_gibbs = ' '.join([f"{gibb}" for gibb in gibbs]),
        equilibrium_constants = ' '.join([f"{eqcon}" for eqcon in equilibrium_constants])))