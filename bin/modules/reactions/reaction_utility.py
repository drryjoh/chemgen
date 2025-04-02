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