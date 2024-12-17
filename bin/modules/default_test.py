import importlib.util

from .configuration import *
from .headers import *
import cantera as ct
from .cmake import *

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
    
    # Extract values from the parsed YAML data
    test_conditions = configuration.get('test_conditions', {})

    # Handle the optional 'random' key
    random = test_conditions.get('random', None)  # Default to None if 'random' is not present
    if random:
        return [0, 0, 0, random]

    # Extract other required values
    temperature = test_conditions['temperature']
    pressure = test_conditions['pressure']
    species_list = test_conditions['species']

    # Print out the parsed values
    print("Test configuration:")
    print(f"Temperature: {temperature}")
    print(f"Pressure: {pressure}")
    print("Species:")

    for species in species_list:
        print(f"  - Name: {species['name']}, MoleFraction: {species['MoleFraction']}")
    species_string  = ' '.join([f"{species['name']}:{species['MoleFraction']} " for species in species_list])
    return [temperature, pressure, species_string, False]

def create_test(gas, chemical_mechanism, headers, test_file_name, configuration, destination_folder, n_points = 0):
    test_file = destination_folder/test_file_name
    with open(test_file, 'w') as file:
        file.write("#include <cmath>\n")
        file.write("#include <algorithm>\n")
        file.write("#include <array>\n")
        file.write("#include <iostream>  // For printing the result to the console\n")
        write_headers(file, headers)

        [temperature, pressure, species_string, random] = get_test_conditions(chemical_mechanism)
        if random:
            import random
            # Generate n random numbers
            random_X = [random.random() for _ in range(gas.n_species)]
            total = sum(random_X)

            # Normalize the numbers so they sum to 1
            random_X = [x / total for x in random_X]
            gas.TPX = 300 + random.random() * 5000, 10132.5 + 1013250 * random.random(), random_X
            temperature = gas.T
            pressure = gas.P
            print("Random state generated")
            print(f"Temperature: {gas.T}")
            print(f"Pressure: {gas.P}")
            print(f"Mole Fractions: {gas.X}")
        else:
            gas.TPX = temperature, pressure, species_string
        concentrations = gas.concentrations

        concentration_test = '{species} species  = {{{array}}};'.format(array = ','.join(["{scalar_cast}({c})".format(c=c, **vars(configuration)) for c in concentrations]),**vars(configuration)) 
        enthalpies = gas.standard_enthalpies_RT * gas.T * ct.gas_constant/gas.molecular_weights
        entropies = gas.standard_entropies_R * ct.gas_constant/gas.molecular_weights
        energies = gas.standard_int_energies_RT * gas.T * ct.gas_constant/gas.molecular_weights
        gibbs = gas.standard_gibbs_RT * gas.T * ct.gas_constant
        gibbs_reactions = gas.delta_standard_gibbs/gas.T/ct.gas_constant
        equilibrium_constants = gas.equilibrium_constants

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

{index} main() {{
    std::cout << "*** ChemGen ***" <<std::endl;
    {concentration_test}
    {scalar} temperature_ =  {temperature};
    {species} result = source(species, temperature_);


    {scalar} pressure_return = pressure(species, temperature_);
    {scalar} int_energy = internal_energy_volume_specific(species, temperature_);
    std::cout << "temperature: " << temperature_<<std::endl;
    for({index} i=0; i<10; ++i)
    {{
        std::cout << "temperature_ for "<< i <<" iterations: " << (temperature_ - temperature(int_energy, species, i)) / (temperature_)<<std::endl;
    }}
    // Output the result
    std::cout << "Source test result:  " << result << std::endl;
    std::cout << "Cantera test result: " <<"{cantera_net_production_rates}"<<std::endl;

    
    std::cout << "ChemGen internal energy: "<< int_energy <<std::endl;
    std::cout << "Cantera internal energy: " <<"{cantera_int_energy}"<<std::endl;

    
    std::cout << "Chemgen species cps: " << species_specific_heat_constant_pressure_mass_specific(temperature_) <<std::endl;
    std::cout << "Cantera species cps: " <<"{cantera_species_cp}"<<std::endl;

    
    std::cout << "Chemgen species enthalpies: " << species_enthalpy_mass_specific(temperature_) <<std::endl;
    std::cout << "Cantera species enthalpies: " <<"{cantera_species_enthalpy}"<<std::endl;
    
    std::cout << "Chemgen species internal energies: " << species_internal_energy_mass_specific(temperature_) <<std::endl;
    std::cout << "Cantera species internal energies: " <<"{cantera_species_energies}"<<std::endl;

    std::cout << "Chemgen species internal entropies: " << species_entropy_mass_specific(temperature_) <<std::endl;
    std::cout << "Cantera species internal entropies: " <<"{cantera_species_entropy}"<<std::endl;

    
    std::cout << "Chemgen species gibbs energy: " << species_gibbs_energy_mole_specific(temperature_) <<std::endl;
    std::cout << "Cantera species gibbs energy: " <<"{cantera_species_gibbs}"<<std::endl;
    
    std::cout << "Chemgen species gibbs reactions : " << gibbs_reaction(log_gen(temperature_)) <<std::endl;
    std::cout << "Cantera species gibbs reactions : " <<"{cantera_gibbs_reactions}"<<std::endl;

    std::cout << "Pressure: " <<pressure_return <<std::endl;
    std::cout << "Temperature Monomial at 300           : " <<temperature_monomial({scalar_cast}(300)) <<std::endl;
    std::cout << "Temperature Energy Monomial at 300           : " <<temperature_energy_monomial({scalar_cast}(300)) <<std::endl;
    std::cout << "Temperature Entropy Monomial at 300           : " <<temperature_entropy_monomial({scalar_cast}(300)) <<std::endl;
    std::cout << "Temperature Gibbs Monomial at 300           : " <<temperature_gibbs_monomial({scalar_cast}(300)) <<std::endl;
    
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
        cantera_gibbs_reactions = ' '.join([f"{gibb}" for gibb in gibbs_reactions]),
        equilibrium_constants = ' '.join([f"{eqcon}" for eqcon in equilibrium_constants]),
        cantera_int_energy = gas.int_energy_mass * gas.density))
