import yaml
from pathlib import Path

class Config:
    pass

def add_derived_attributes(config_obj):
    """ adds derived attributes from configuration decorators """
    setattr(config_obj, "device_sum_inl",  "{device_sum}".format(**vars(config_obj)).format(**vars(config_obj)))
    setattr(config_obj, "device_fill_inl", "{device_fill}".format(**vars(config_obj)).format(**vars(config_obj)))
    setattr(config_obj, "device_transform_inl", "{device_transform}".format(**vars(config_obj)).format(**vars(config_obj)))
    setattr(config_obj, "device_element_sum_inl", "{device_element_sum}".format(**vars(config_obj)).format(**vars(config_obj)))
    setattr(config_obj, "device_element_multiply_inl", "{device_element_multiply}".format(**vars(config_obj)).format(**vars(config_obj)))
    setattr(config_obj, "device_element_sum_offset_inl", config_obj.device_element_sum_inl.replace("a_s.begin()", "a_s.begin()+1").replace("c_s.begin()", "c_s.begin()+1"))

def get_configuration(configuration_filename = 'configuration.yaml', decorators = 'decorators'):
    """ opens and returns configuration file """
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
    """ tests configuraiton file compatability with some hard-set """
    if args.ignore_temp_dependence:
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
