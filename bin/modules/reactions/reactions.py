from .third_body import *
from .arrhenius import *
from .falloff import *
from .pressure_dependent_arrhenius import *
def write_reaction_rates(file, reaction_rates):
    for reaction in reaction_rates:
        file.write(f"    {reaction}\n")
def write_reaction_rates_derivatives(file, reaction_rates_derivatives):
    for reaction in reaction_rates_derivatives:
        file.write(f"    {reaction}\n")
