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
def check_configuration(configuration, temperature_jacobian, force = False): 
    if temperature_jacobian:
        jacobian_type  = f"{configuration.jacobian_typedef}"
        if "Species, n_species" in jacobian_type and not force:
            exit(f"{jacobian_type} is probably incorrect for jacobian_typedef in configuraiton file\n Consider one with size <n_species +1, n_species + 1> such as std::array<ChemicalState, n_species + 1>\n to continue use --force")
        elif "Species, n_species" in jacobian_type and force:
            print(f"{jacobian_type} is probably incorrect ")
        setattr(configuration, "temperature_jacobian", "on")
    else:
        setattr(configuration, "temperature_jacobian", "off")



def get_default_configuration():
    current_dir = Path(__file__).resolve().parent
    configuration_filename = current_dir.parent/ 'configuration.yaml'
    print("** No configuration file detected, using decorators defaults in /bin/configuration.yaml **")

    return configuration_filename
