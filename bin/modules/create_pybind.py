import os
import re
from pathlib import Path
def pybind_extract_functions_from_header(file_path):
    FUNC_PATTERN = re.compile(r'^\s*(?:inline\s+)?([\w:<>]+)\s+(\w+)\s*\(([^)]*)\)\s*(?:;|\{)?',  re.MULTILINE)
    with open(file_path, "r") as f:
        content = f.read()
        return FUNC_PATTERN.findall(content)

def pybind_format_binding(return_type, func_name, args):
    return f'    m.def("{func_name}", &{func_name}, "{func_name} function");'

def create_pybind(gas, headers, configuration, destination_folder, remove_reactions = False):
    includes = headers
    bindings = []
    bindings_file = Path(destination_folder) / Path("chemgen_pybind.cpp")     # generated binding file\
    setup_file = Path(destination_folder) / Path("setup_chemgen.py")     # generated binding file

    for path in Path(destination_folder).rglob("*.h"):
        for return_type, func_name, args in pybind_extract_functions_from_header(path):
            bindings.append(pybind_format_binding(return_type, func_name, args))

    with open(bindings_file, "w") as f:
        f.write('#include <pybind11/pybind11.h>\n')
        f.write('#include <pybind11/stl.h>\n\n')
        f.write("#include <cmath>\n")
        f.write("#include <algorithm>\n")
        f.write("#include <array>\n")
        f.write("#include <iostream>  // For printing the result to the console\n")
        for header in headers:
            f.write(f"#include \"{header}\"\n")
        remove_reactions_text = ""
        remove_reactions_call_text = "" 
        if remove_reactions:
            remove_reactions_text = """
std::vector<std::vector<{scalar}>> source_jacobian_remove_reaction_py(const std::vector<{scalar}>& species, {scalar} temperature, std::vector<std::vector<{scalar}>>& J, {index} i) 
{{
    Species sp;
    std::copy(species.begin(), species.end(), sp.begin());

    {jacobian} Js;
    for (int i = 0; i < n_species; ++i)
        for (int j = 0; j < n_species; ++j)
            Js[i][j] = J[i][j];  // Manual element-wise copy

    {jacobian} jac = source_jacobian_R(sp, temperature, Js, i);

    for (int i = 0; i < n_species; ++i) std::copy(J[i].begin(), J[i].end(), Js[i].begin());

    std::vector<std::vector<{scalar}>> jac_out(n_species, std::vector<{scalar}>(n_species));
    for (int i = 0; i < n_species; ++i)
        for (int j = 0; j < n_species; ++j)
            jac_out[i][j] = jac[i][j];

    return jac_out;
}}
            """.format(**vars(configuration))
            remove_reactions_call_text = """    m.def("source_jacobian_remove_reaction", &source_jacobian_remove_reaction_py, "source_jacobian with reaction removed function");"""


        f.write("""
namespace py = pybind11;

std::vector<{scalar}> source_py(const std::vector<{scalar}>& species, {scalar} temperature) 
{{
    Species sp;
    std::copy(species.begin(), species.end(), sp.begin());
    auto result = source(sp, temperature);
    return std::vector<{scalar}>(result.begin(), result.end());
}}

std::vector<std::vector<{scalar}>> source_jacobian_py(const std::vector<{scalar}>& species, {scalar} temperature) 
{{
    Species sp;
    std::copy(species.begin(), species.end(), sp.begin());

    SpeciesJacobian jac = source_jacobian(sp, temperature);

    std::vector<std::vector<{scalar}>> jac_out(n_species, std::vector<{scalar}>(n_species));
    for (int i = 0; i < n_species; ++i)
        for (int j = 0; j < n_species; ++j)
            jac_out[i][j] = jac[i][j];

    return jac_out;
}}

std::vector<{scalar}> sdirk4_py(const std::vector<{scalar}>& species, {scalar} temperature, {scalar} dt, {scalar} norm, {index} max_iter) 
{{
    Species sp;
    std::copy(species.begin(), species.end(), sp.begin());
    {scalar} int_energy = internal_energy_volume_specific(sp, temperature);
    {chemical_state} y = set_chemical_state(int_energy, sp);
    
    auto result = sdirk4(y, dt, norm, max_iter);

    return std::vector<{scalar}>(result.begin(), result.end());
}}

{scalar} temperature_from_internal_energy_py(const std::vector<{scalar}>& species, {scalar} internal_energy) 
{{
    Species sp;
    std::copy(species.begin(), species.end(), sp.begin());
    return temperature(internal_energy, sp);
}}


{remove_reactions}
PYBIND11_MODULE(chemgen, m)
{{ 
    m.def("source", &source_py, "source function");
    m.def("source_jacobian", &source_jacobian_py, "source_jacobian function");
    m.def("sdirk4", &sdirk4_py, "sdirk 4");
    m.def("temperature_from_internal_energy", &temperature_from_internal_energy_py, "temperature 4");
    {remove_reactions_call}
}}

        """.format(**vars(configuration), remove_reactions = remove_reactions_text, remove_reactions_call = remove_reactions_call_text))
    
    with open(setup_file, "w") as f:
        f.write("""
from setuptools import setup, Extension
import sys
import pybind11

ext_modules = [
    Extension(
        'chemgen',
        sources=['src/chemgen_pybind.cpp'],
        include_dirs=[pybind11.get_include()],
        language='c++',
        extra_compile_args=['-std=c++17']
    )
]

setup(
    name='chemgen',
    version='1.0',
    author='Ryan F. Johnson',
    description='Pybind11 bindings for C++ ChemGen',
    ext_modules=ext_modules,
)
        """)
    

    print(f"[+] Generated {bindings_file} with {len(bindings)} bindings.")

