from .third_body import *
from .arrhenius import *
from .falloff import *
from .pressure_dependent_arrhenius import *
def write_reaction_rates(file, reaction_rates):
    for reaction in reaction_rates:
        file.write(f"    {reaction}\n")
