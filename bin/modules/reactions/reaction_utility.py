import numpy as np
def get_efficiencies(reaction):
    import cantera as ct
    cantera_version = ct.__version__.split('.')
    major_version = cantera_version[0]
    minor_version = cantera_version[1]
    efficiencies = None
    if major_version != '3':
        print("We support cantera versions >3 please install cantera >3\n pip3 install cantera 3.0.0")
    if float(minor_version) > 0:
        efficiencies = reaction.third_body.efficiencies
    else:
        efficiencies = reaction.efficiencies
    return efficiencies
def get_mixture_concentration(efficiencies, species_names, configuration):
    mixture_concentration_array = []
    # TODO: Need to account for default efficiency
    if all(np.abs(eff-1.0) < 0.001 for eff in efficiencies.values()) or len(efficiencies) == 0:
        return "mixture_concentration"
    else:
        if not efficiencies:
            efficiencies = {specie: 1.0 for specie in species_names}
        for species_index, specie in enumerate(species_names):
            if specie in efficiencies:
                if efficiencies[specie] != 1:
                    mixture_concentration_array.append(f"({configuration.scalar_cast}({efficiencies[specie]})-{configuration.scalar_cast}(1))*{configuration.species_element.format(i = species_index)}")

        return "mixture_concentration + {0}".format(' + '.join(mixture_concentration_array))

def get_default_efficiency(reaction):
    try:
        return reaction.third_body.default_efficiency
    except AttributeError:
        return reaction.default_efficiency

def get_mixture_concentration_derivatives(reaction, efficiencies, species_names, configuration):
    mixture_concentration_array = ['{scalar_cast}({default_efficiency})'.format(**vars(configuration), default_efficiency=get_default_efficiency(reaction))] * len(species_names)
    if all(np.abs(eff-1.0) < 0.001 for eff in efficiencies.values()) or len(efficiencies) == 0:
        return "{0}".format(','.join(mixture_concentration_array))
    else:
        if not efficiencies:
            efficiencies = {specie: 1.0 for specie in species_names} # TODO: is this necessary?
        for species_index, specie in enumerate(species_names): # TODO: species is both singular and plural
            if specie in efficiencies:
                mixture_concentration_array[species_index] = f"{configuration.scalar_cast}({efficiencies[specie]})"
        return "{0}".format(', '.join(mixture_concentration_array))

def get_mixture_concentration_derivatives_array(reaction, efficiencies, species_names, configuration):
    array = np.full(len(species_names), get_default_efficiency(reaction))

    if len(efficiencies) == 0 or not efficiencies:
        return array
    else:
        for species_index, specie in enumerate(species_names):
            if specie in efficiencies:
                array[species_index] = efficiencies[specie]
        return array
