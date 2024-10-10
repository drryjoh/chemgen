import yaml

class Config:
    pass

def get_configuration(configuration_filename='configuration_header.yaml', decorators = 'decorators'):
    with open(configuration_filename, 'r') as file:
        config = yaml.safe_load(file)
    # Create an instance of Config and set attributes dynamically
    config_obj = Config()
    for key, value in config[decorators].items():
        setattr(config_obj, key, value)
    
    print("** decorators used **")
    for key, value in vars(config_obj).items():
        print(f"{key}: {value}")
    return config_obj    