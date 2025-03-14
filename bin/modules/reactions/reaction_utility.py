import numpy as np
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

def get_mixture_concentration_derivatives(efficiencies, species_names, configuration):
    mixture_concentration_array = ['{scalar_cast}(1)'.format(**vars(configuration))] * len(species_names)
    if all(np.abs(eff-1.0) < 0.001 for eff in efficiencies.values()) or len(efficiencies) == 0:
        return "{0}".format(','.join(mixture_concentration_array))
    else:
        if not efficiencies:
            efficiencies = {specie: 1.0 for specie in species_names}
        for species_index, specie in enumerate(species_names):
            if specie in efficiencies:
                if efficiencies[specie] != 1:
                    mixture_concentration_array[species_index] = f"{configuration.scalar_cast}({efficiencies[specie]})"
        return "{0}".format(', '.join(mixture_concentration_array))