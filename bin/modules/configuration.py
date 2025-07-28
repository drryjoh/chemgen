import yaml
from pathlib import Path

class Config:
    pass

def add_derived_attributes(config_obj):
    setattr(config_obj, "device_sum_inl",  "{device_sum}".format(**vars(config_obj)).format(**vars(config_obj)))
    setattr(config_obj, "device_fill_inl", "{device_fill}".format(**vars(config_obj)).format(**vars(config_obj)))
    setattr(config_obj, "device_transform_inl", "{device_transform}".format(**vars(config_obj)).format(**vars(config_obj)))
    setattr(config_obj, "device_element_sum_inl", "{device_element_sum}".format(**vars(config_obj)).format(**vars(config_obj)))
    setattr(config_obj, "device_element_multiply_inl", "{device_element_multiply}".format(**vars(config_obj)).format(**vars(config_obj)))

    setattr(config_obj, "device_element_sum_offset_inl", config_obj.device_element_sum_inl.replace("a_s.begin()", "a_s.begin()+1").replace("c_s.begin()", "c_s.begin()+1"))

def get_configuration(configuration_filename = 'configuration.yaml', decorators = 'decorators'):
    config_path = Path(configuration_filename)
    
    if config_path.exists():
        with config_path.open('r') as file:
            configuration = yaml.safe_load(file)
    else:
        default_configuration = get_default_configuration()
        with default_configuration.open('r') as file:
            configuration = yaml.safe_load(file)
    # Create an instance of Config and set attributes dynamically
    config_obj = Config()
    for key, value in configuration[decorators].items():
        setattr(config_obj, key, value)
    
    add_derived_attributes(config_obj)
    return [config_obj, configuration]

#update with other checks later on
def check_configuration(configuration, args):
    if args.ignore_temp_dependence:
        # jacobian_type  = f"{configuration.jacobian_typedef}"
        # if "Species, n_species" in jacobian_type and not force:
        #     exit(f"{jacobian_type} is probably incorrect for jacobian_typedef in configuraiton file\n Consider one with size <n_species +1, n_species + 1> such as std::array<ChemicalState, n_species + 1>\n to continue use --force")
        # elif "Species, n_species" in jacobian_type and force:
        #     print(f"{jacobian_type} is probably incorrect ")
        setattr(configuration, "temperature_jacobian", "off")
    else:
        setattr(configuration, "temperature_jacobian", "on")

    jacobian_type  = f"{configuration.jacobian_typedef}"
    jacobian_type_no_whitespace = jacobian_type.replace(" ", "")

    if jacobian_type_no_whitespace == "std::array<Species,n_species>":
        # Directly modify jacobian_type if straightforward
        configuration.jacobian_typedef = "std::array<ChemicalState, n_species + 1>"
        print("Changing jacobian_type from {} to {}".format(jacobian_type, configuration.jacobian_typedef))
    elif "Species" in jacobian_type_no_whitespace or ("n_species" in jacobian_type_no_whitespace and "n_species+" not in jacobian_type_no_whitespace):
        exit(f"{jacobian_type} is probably incorrect for jacobian_typedef in configuration file\n Consider one with size <n_species +1, n_species + 1> such as std::array<ChemicalState, n_species + 1>\n")

    jacobian_end  = f"{configuration.jacobian_end}"
    jacobian_end_no_whitespace = jacobian_end.replace(" ", "")
    if jacobian_end_no_whitespace != "n_variables" and jacobian_end_no_whitespace != "n_species+1":
        exit("jacobian_end should be 'n_variables' or 'n_species + 1'\n")

def get_default_configuration():
    current_dir = Path(__file__).resolve().parent
    configuration_filename = current_dir.parent/ 'configuration.yaml'
    print("** No configuration file detected, using decorators defaults in /bin/configuration.yaml **")

    return configuration_filename

def update_configuration_eigen(configuration):
    if not configuration.eigen:
        species_eigen = configuration.species
        species_function_eigen = configuration.species_function
        species_parameter_eigen = configuration.species_parameter

        chemical_state_eigen = configuration.chemical_state
        chemical_state_function_eigen = configuration.chemical_state_function
        chemical_state_parameter_eigen = configuration.chemical_state_parameter

        jacobian_eigen = configuration.jacobian
        jacobian_function_eigen = configuration.jacobian_function
        jacobian_parameter_eigen = configuration.jacobian_parameter
        jacobian_eigen_dynamic = configuration.jacobian
        jacobian_eigen_dense = configuration.jacobian

        left = "["
        mid = "]["
        right = "]"

        if configuration.fixed_jacobian:
            raise NotImplementedError

        linear_solver_direct_eigen = ""
        linear_solver_iterative_eigen = ""
        linear_solver_eigen = ""

        declare_linear_solver_direct_eigen = ""
        declare_linear_solver_iterative_eigen = ""
        declare_linear_solver_eigen = ""

        declare_compute_flag = ""
        initialize_compute_flag = ""
        update_compute_flag = ""
        if_compute_flag = ""
        solver_and_compute_flag_arguments = ""
    else:
        scalar = configuration.scalar
        index = configuration.index

        species_eigen = f"Matrix<{scalar}, n_species, 1>"
        species_function_eigen = f"Matrix<{scalar}, n_species, 1>"
        species_parameter_eigen = f"const Matrix<{scalar}, n_species, 1>&"

        chemical_state_eigen = f"Matrix<{scalar}, n_variables, 1>"
        chemical_state_function_eigen = f"Matrix<{scalar}, n_variables, 1>"
        chemical_state_parameter_eigen = f"const Matrix<{scalar}, n_variables, 1>&"

        jacobian_eigen = f"Matrix<{scalar}, n_variables, n_variables>"
        jacobian_function_eigen = f"Matrix<{scalar}, n_variables, n_variables>"
        jacobian_parameter_eigen = f"const Matrix<{scalar}, n_variables, n_variables>&"
        jacobian_eigen_dynamic = f"Matrix<{scalar}, Dynamic, Dynamic>"
        jacobian_eigen_dense = jacobian_eigen

        left = "("
        mid = ","
        right = ")"

        linear_solver_direct_eigen = f"PartialPivLU<{jacobian_eigen}>"

        if configuration.eigen_sparse:
            jacobian_eigen = f"SparseMatrix<{scalar}>"
            jacobian_function_eigen = f"SparseMatrix<{scalar}>"
            jacobian_parameter_eigen = f"const SparseMatrix<{scalar}>&"
            jacobian_eigen_dynamic = f"SparseMatrix<{scalar}>"

            linear_solver_direct_eigen = f"SparseLU<SparseMatrix<{scalar}>, NaturalOrdering<{index}>>"

        declare_linear_solver_direct_eigen = f"{linear_solver_direct_eigen} solver;"

        linear_solver_iterative_eigen = f"GMRES<{jacobian_eigen_dynamic}, Preconditioner>"
        declare_linear_solver_iterative_eigen = f"{linear_solver_iterative_eigen} solver;"

        if not configuration.direct_solver:
            linear_solver_eigen = linear_solver_iterative_eigen
            declare_linear_solver_eigen = declare_linear_solver_iterative_eigen
        else:
            linear_solver_eigen = linear_solver_direct_eigen
            declare_linear_solver_eigen = declare_linear_solver_direct_eigen

        if configuration.fixed_jacobian:
            declare_compute_flag = "bool compute_jacobian;"
            initialize_compute_flag = "compute_jacobian = true;"
            update_compute_flag = """
            if (iter > 0) compute_jacobian = false;"""

            if_compute_flag = "if (compute_jacobian)"

            solver_and_compute_flag_arguments = ", solver, compute_jacobian"
        else:
            declare_compute_flag = ""
            initialize_compute_flag = ""
            update_compute_flag = ""

            if_compute_flag = ""

            solver_and_compute_flag_arguments = ", solver, true"


    setattr(configuration, "species_eigen", species_eigen)
    setattr(configuration, "species_function_eigen", species_function_eigen)
    setattr(configuration, "species_parameter_eigen", species_parameter_eigen)

    setattr(configuration, "chemical_state_eigen", chemical_state_eigen)
    setattr(configuration, "chemical_state_function_eigen", chemical_state_function_eigen)
    setattr(configuration, "chemical_state_parameter_eigen", chemical_state_parameter_eigen)

    setattr(configuration, "jacobian_eigen", jacobian_eigen)
    setattr(configuration, "jacobian_function_eigen", jacobian_function_eigen)
    setattr(configuration, "jacobian_parameter_eigen", jacobian_parameter_eigen)
    setattr(configuration, "jacobian_eigen_dynamic", jacobian_eigen_dynamic)
    setattr(configuration, "jacobian_eigen_dense", jacobian_eigen_dense)

    setattr(configuration, "left", left)
    setattr(configuration, "mid", mid)
    setattr(configuration, "right", right)

    setattr(configuration, "linear_solver_direct_eigen", linear_solver_direct_eigen)
    setattr(configuration, "linear_solver_iterative_eigen", linear_solver_iterative_eigen)
    setattr(configuration, "linear_solver_eigen", linear_solver_eigen)

    setattr(configuration, "declare_linear_solver_direct_eigen", declare_linear_solver_direct_eigen)
    setattr(configuration, "declare_linear_solver_iterative_eigen", declare_linear_solver_iterative_eigen)
    setattr(configuration, "declare_linear_solver_eigen", declare_linear_solver_eigen)

    setattr(configuration, "declare_compute_flag", declare_compute_flag)
    setattr(configuration, "initialize_compute_flag", initialize_compute_flag)
    setattr(configuration, "update_compute_flag", update_compute_flag)
    setattr(configuration, "if_compute_flag", if_compute_flag)
    setattr(configuration, "solver_and_compute_flag_arguments", solver_and_compute_flag_arguments)

