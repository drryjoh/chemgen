import yaml
from pathlib import Path

class Config:
    pass

def get_configuration(configuration_filename='configuration_header.yaml', decorators = 'decorators'):
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
    
    print("** decorators used **")
    for key, value in vars(config_obj).items():
        print(f"{key}: {value}")
    return config_obj 

def get_default_configuration():
    current_dir = Path(__file__).resolve().parent
    configuration_filename = current_dir.parent/ 'configuration.yaml'
    print("** No configuration file detected, using decorators defaults in /bin/configuration.yaml **")
    return configuration_filename
