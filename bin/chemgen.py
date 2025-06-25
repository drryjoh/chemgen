#!python3

import argparse
import sys
import cantera as ct
import numpy as np
import shutil
from pathlib import Path
from modules.process_files import *
from modules.configuration import *
from modules.headers import *
from modules.compile_and_run import *
from modules.create_pybind import *
 
def find_chemical_mechanism(file_name):
    # Get the path to the current script (chemgen.py)
    script_path = Path(__file__).resolve()
    # Get the parent directory of the script (assumed to be /path_to_code/bin)
    base_directory = script_path.parent.parent
    # Define the path to the chemical_mechanisms directory
    chemical_mechanisms_dir = base_directory / 'chemical_mechanisms'

    # Check if the file has the .yaml extension; if not, append it
    if not file_name.endswith('.yaml'):
        file_name += '.yaml'

    # Convert the input file_name into a Path object
    chemical_mechanism = Path(file_name)

    # Case 1: Direct path provided by the user
    if chemical_mechanism.is_file():
        return chemical_mechanism

    # Case 2: The file is located in the chemical_mechanisms directory
    chemical_mechanism_in_dir = chemical_mechanisms_dir / chemical_mechanism.name
    if chemical_mechanism_in_dir.is_file():
        return chemical_mechanism_in_dir

    # If neither case works, print an error and return None
    print(f"Error: The file '{file_name}' does not exist.")
    return None

# Define functions or classes here
def main():
    cantera_version = ct.__version__
    major_version = cantera_version[0]
    if major_version != '3':
        print("We support cantera versions >3 please install cantera >3\n pip3 install cantera 3.0.0")
        exit()
    """
    Main function where the core logic of the script runs.
    """
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="A brief description of the script")
    parser.add_argument('chemical_mechanism', type=str, help="The name of the file")
    parser.add_argument('destination', type=str, help="The destination folder where the files will be generated")
    parser.add_argument("--custom-source", type=str, help="Path to custom source writer file")
    parser.add_argument("--custom-test", type=str, help="Path to custom test case writer file")
    parser.add_argument("--compile", action="store_true", help="Compile the source writer code")
    parser.add_argument("--cmake", action="store_true", help="Compile the source writer code")
    parser.add_argument("--n-points-test", type=int, default=1000,  help="Number of points for testing (default: 1000)")
    parser.add_argument("--verbose", action="store_true", help="Verbose code generation")
    parser.add_argument("--remove_reactions", action="store_true", help="Generate the ability to remove single reaction from jacobian")
    parser.add_argument("--fit-gibbs-reaction", action="store_true", default=True, help="Fit the gibbs free energy per reaction") # TODO: this can never be False
    parser.add_argument("--ignore-temp-dependence", action="store_true", help="Ignore temperature dependence in Jacobian if internal energy is state variable")
    parser.add_argument("--force", action="store_true", help="Force code generation despite warnings")
    parser.add_argument("--pybind", action="store_true", help="Create pybind linker")
    parser.add_argument("--skip", action="store_true", help="Skip code generation")
    parser.add_argument("--skip-tests", action="store_true", help="Skip running tests")
    parser.add_argument("--get-sparsity", action="store_true", help="Calculate Jacobian sparsity")
    parser.add_argument("--plot-sparsity", action="store_true", help="Plot Jacobian sparsity pattern")
    parser.add_argument("--temperature-equation", action="store_true", help="Solve temperature equation instead of assuming constant internal energy")
    parser.add_argument("--ignore-other-species", action="store_true", help="Ignore third-body and pressure dependence when calculating derivatives with respect to species")

    args = parser.parse_args()
    
    # Convert arguments to Path objects
    chemical_mechanism = find_chemical_mechanism(args.chemical_mechanism)
    destination_folder = Path(args.destination)/'src'
    n_points_test = args.n_points_test
    fit_gibbs_reaction  = True
    if args.fit_gibbs_reaction == False:
        fit_gibbs_reaction  = False
        print("Gibbs free energies will be fitted per species and then summation will be performed according to stoicheimetry.\n Warning, this has shown to cause some errors when compared to cantera.")
    
    temperature_jacobian = True
    if args.ignore_temp_dependence:
        if args.temperature_equation:
            exit("Cannot ignore temperature dependence if solving temperature equation")
        temperature_jacobian = False
        print("Source Jacobian will be created without temperature dependence")

    force  = False
    if args.force == True:
        force  = True
        print("ChemGen will continue despite warnings")
    
    # Check if the destination folder exists, if not, create it
    if not destination_folder.exists():
        print(f"Destination folder '{destination_folder}' does not exist. Creating it...")
        destination_folder.mkdir(parents=True, exist_ok=True)

    # Core logic of the script
    print(f"Processing file: {args.chemical_mechanism}")
    

    gas = ct.Solution(chemical_mechanism)
    [configuration, configuration_file] = get_configuration(configuration_filename='configuration.yaml')

    check_configuration(configuration, args)

    if args.temperature_equation:
        setattr(configuration, "internal_energy_or_temperature",  "CHEMGEN_TEMPERATURE_EQUATION")
    else:
        setattr(configuration, "internal_energy_or_temperature",  "CHEMGEN_INTERNAL_ENERGY_EQUATION")

    use_third_parties = False
    third_party_path = Path(__file__).resolve().parent.parent/'third_party'
    libraries = []
    if configuration_file['build'].get('chemgen_smp','').lower() == 'tbb':
        use_third_parties = True
        libraries.append('tbb')
    if configuration_file['build'].get('chemgen_mpi',''):
        use_third_parties = True
        libraries.append('mpi')

    generate_chemistry_solver = False
    chemistry_solver = configuration_file.get('solver', {}).get('chemistry_solver', None)
    chemistry_solver_preconditioner = configuration_file.get('solver', {}).get('preconditioner', None)
    direct_solver = ''
    if chemistry_solver:
        generate_chemistry_solver = True
        if chemistry_solver.lower() == "none":
            print("none was specified for solver: chemistry_solver in configuration file, no solver will be generated")
            generate_chemistry_solver = False
        elif chemistry_solver.lower() == "rk4":
            print("RK4 chemistry solver chosen")
        elif chemistry_solver.lower() == "backwards_euler":
            linear_solver = configuration_file.get('solver', {}).get('linear_solver', None)
            print("Backwards Euler chemistry solver chosen")
            if linear_solver!=None and linear_solver.lower() == "gmres":
                print("GMRES linear solver chosen")
            elif linear_solver!=None and linear_solver.lower() == "direct":
                print("Direct solver chosen for Jacobian inversion")
                direct_solver = "#define CHEMGEN_DIRECT_SOLVER"
            else:
                print("linear solver not recognized, defaulting to GMRES")
        elif chemistry_solver.lower() == "all":
            linear_solver = configuration_file.get('solver', {}).get('linear_solver', None)
            print("All solver options will be compiled in")
            if linear_solver!=None and linear_solver.lower() == "gmres":
                print("GMRES linear solver chosen")
            elif linear_solver!=None and linear_solver.lower() == "direct":
                print("Direct solver chosen for Jacobian inversion")
                direct_solver = "#define CHEMGEN_DIRECT_SOLVER"
            else:
                print("linear solver not recognized, defaulting to GMRES")
        else:
            print("Chemistry solver unsupported. Please choose from [rk4, backwards_euler].")
            exit()
    else:
        generate_chemistry_solver = False
        print("Not generating with a chemgen chemistry solver.")
    
    preconditioner = ""
    if chemistry_solver_preconditioner and chemistry_solver:
        
        if chemistry_solver_preconditioner.lower() == "none":
            print("none was specified for solver: preconditioner in configuration file, no preconditioner will be used")
            preconditioner = ""
        elif chemistry_solver_preconditioner.lower() == "gauss_seidel":
            preconditioner = "#define CHEMGEN_PRECONDITIONER_GAUSS_SEIDEL"
        elif chemistry_solver_preconditioner.lower() == "jacobi":
            preconditioner = "#define CHEMGEN_PRECONDITIONER_JACOBI"
        else:
            print("Chemistry solver preconditioner unsupported. Please choose from [none, gauss_seidel, jacobi].")
            exit()
    
    setattr(configuration, "preconditioner",  preconditioner)
    setattr(configuration, "direct_solver",  direct_solver)

    third_parties = [use_third_parties, third_party_path, libraries]

    test_file = 'chemgen.cpp'

    if args.skip:
        print("Skipping code generation")

        if args.get_sparsity:
            print("Warning: cannot get sparsity pattern if skipping code generation")
    else:
        if args.get_sparsity:
            if temperature_jacobian or args.remove_reactions:
                raise NotImplementedError

            sparsity_pattern = np.zeros([gas.n_species, gas.n_species], dtype=int)
        else:
            sparsity_pattern = None
            if args.plot_sparsity:
                print("Warning: Cannot plot sparsity pattern if no get-sparsity argument")

        headers = process_cantera_file(gas, configuration, destination_folder,args, chemistry_solver, verbose = args.verbose, fit_gibbs_reaction = fit_gibbs_reaction, temperature_jacobian = temperature_jacobian, remove_reactions = args.remove_reactions, sparsity_pattern=sparsity_pattern)

        if args.get_sparsity:
            n_nonzeros = np.count_nonzero(sparsity_pattern)
            n_entries = sparsity_pattern.size
            n_zeros = n_entries - n_nonzeros
            sparsity_percent = 100.*(1. - n_nonzeros/n_entries)
            print("Jacobian: # nonzero entries = {}, # zero entries = {}, total # entries = {}".format(n_nonzeros, n_zeros, n_entries))
            print("Sparsity percentage = %g%%" % sparsity_percent)
            if args.plot_sparsity:
                import matplotlib.pyplot as plt
                plt.figure()
                plt.spy(sparsity_pattern)
                plt.show()

        if "types_inl.h" in headers:
            headers.remove("types_inl.h")
            headers.insert(0,"types_inl.h")

        if "chemical_state_functions.h" in headers:
            headers.remove("chemical_state_functions.h")
            headers.append("chemical_state_functions.h")

        if "rk4.h" in headers:
            headers.remove("rk4.h")
            headers.append("rk4.h")

        if "linear_solvers.h" in headers:
            headers.remove("linear_solvers.h")
            headers.append("linear_solvers.h")

        if "backwards_euler.h" in headers:
            headers.remove("backwards_euler.h")
            headers.append("backwards_euler.h")

        if "sdirk.h" in headers:
            headers.remove("sdirk.h")
            headers.append("sdirk.h")

        if "rosenbroc.h" in headers:
            headers.remove("rosenbroc.h")
            headers.append("rosenbroc.h")

        if "yass.h" in headers:
            headers.remove("yass.h")
            headers.append("yass.h")

        if args.custom_test:
            try:
                # Load the custom SourceWriter
                create_test = load_custom_test(args.custom_test)
                create_test(gas, args.chemical_mechanism, headers, test_file, configuration, destination_folder, n_points = n_points_test)

            except (FileNotFoundError, AttributeError) as e:
                print(f"Error loading custom test writer: {e}")
                sys.exit(1)
        else: #replace with run time argument
            from modules.default_test import create_test
            create_test(gas, args.chemical_mechanism, headers, test_file, configuration, destination_folder)

    if args.pybind:
        create_pybind(gas, headers, args, configuration, destination_folder, remove_reactions = args.remove_reactions)
    else:
        compile_and_run(test_file, configuration_file, destination_folder, third_parties, args.cmake, args.compile, args.skip_tests)

if __name__ == "__main__":
    main()

