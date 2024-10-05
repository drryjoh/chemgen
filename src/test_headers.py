import yaml
import textwrap

class Config:
    pass

def get_configuration():
    with open('configuration.yaml', 'r') as file:
        config = yaml.safe_load(file)
    # Create an instance of Config and set attributes dynamically
    config_obj = Config()
    for key, value in config['decorators'].items():
        print(key)
        setattr(config_obj, key, value)

    return config_obj    

f = open('reaction_headers/arrhenius.h.in','r')
config = get_configuration()

with f as file:
    content = file.read()

print(content)
print(vars(config))
new_content = content.format(**vars(config))
#new_content  = textwrap.dedent(eval(f'f"""{content}"""')).strip()

print(new_content)
