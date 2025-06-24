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

def create_pybind(gas, headers, args, configuration, destination_folder, remove_reactions = False):
    includes = headers
    bindings = []
    bindings_file = Path(destination_folder) / Path("chemgen_pybind.cpp")     # generated binding file\
    setup_file = Path(destination_folder) / Path("setup_chemgen.py")     # generated binding file

    for path in Path(destination_folder).rglob("*.h"):
        for return_type, func_name, args_ in pybind_extract_functions_from_header(path):
            bindings.append(pybind_format_binding(return_type, func_name, args_))

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
    for (int i = 0; i < n_variables; ++i)
        for (int j = 0; j < n_variables; ++j)
            Js[i][j] = J[i][j];  // Manual element-wise copy

    {jacobian} jac = source_jacobian_R(sp, temperature, Js, i);

    for (int i = 0; i < n_variables; ++i) std::copy(J[i].begin(), J[i].end(), Js[i].begin());

    std::vector<std::vector<{scalar}>> jac_out(n_variables, std::vector<{scalar}>(n_variables));
    for (int i = 0; i < n_variables; ++i)
        for (int j = 0; j < n_variables; ++j)
            jac_out[i][j] = jac[i][j];

    return jac_out;
}}
            """.format(**vars(configuration))
            remove_reactions_call_text = """    m.def("source_jacobian_remove_reaction", &source_jacobian_remove_reaction_py, "source_jacobian with reaction removed function");"""


        f.write("""
namespace py = pybind11;

bool ignore_temp_dependence_py()
{{
    return {ignore_temp_dependence};
}}

bool temperature_equation_py()
{{
    return {temperature_equation};
}}

#define {internal_energy_or_temperature}

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

    std::vector<std::vector<{scalar}>> jac_out(n_variables, std::vector<{scalar}>(n_variables));
    for (int i = 0; i < n_variables; ++i)
        for (int j = 0; j < n_variables; ++j)
            jac_out[i][j] = jac[i][j];

    return jac_out;
}}


std::vector<{scalar}> sdirk4_py(const std::vector<{scalar}>& species, {scalar} temperature, {scalar} dt, {scalar} norm, {index} max_iter, {scalar} linear_abs_tol, {scalar} linear_rel_tol) 
{{
    Species sp;
    std::copy(species.begin(), species.end(), sp.begin());
#if defined(CHEMGEN_INTERNAL_ENERGY_EQUATION)
    {scalar} int_energy = internal_energy_volume_specific(sp, temperature);
    {chemical_state} y = set_chemical_state(int_energy, sp);
#else
    {chemical_state} y = set_chemical_state(temperature, sp);
#endif
    
    auto result = sdirk4(y, dt, norm, max_iter, linear_abs_tol, linear_rel_tol);

    return std::vector<{scalar}>(result.begin(), result.end());
}}

std::vector<{scalar}> rosenbroc_py(const std::vector<{scalar}>& species, {scalar} temperature, {scalar} dt, {scalar} linear_abs_tol, {scalar} linear_rel_tol) 

{{
    Species sp;
    std::copy(species.begin(), species.end(), sp.begin());
#if defined(CHEMGEN_INTERNAL_ENERGY_EQUATION)
    {scalar} int_energy = internal_energy_volume_specific(sp, temperature);
    {chemical_state} y = set_chemical_state(int_energy, sp);
#else
    {chemical_state} y = set_chemical_state(temperature, sp);
#endif
    
    auto result = rosenbroc(y, dt, linear_abs_tol, linear_rel_tol);

    return std::vector<{scalar}>(result.begin(), result.end());
}}


std::vector<{scalar}> yass_py(const std::vector<{scalar}>& species, {scalar} temperature, {scalar} dt, {scalar} max_norm, {scalar} min_dt, {index} max_iter, {scalar} linear_abs_tol, {scalar} linear_rel_tol) 
{{
    Species sp;
    std::copy(species.begin(), species.end(), sp.begin());
#if defined(CHEMGEN_INTERNAL_ENERGY_EQUATION)
    {scalar} int_energy = internal_energy_volume_specific(sp, temperature);
    {chemical_state} y = set_chemical_state(int_energy, sp);

#else
    {chemical_state} y = set_chemical_state(temperature, sp);
#endif

    auto result = yass(y, dt, max_norm, min_dt, max_iter, linear_abs_tol, linear_rel_tol);

    return std::vector<{scalar}>(result.begin(), result.end());
}}

#if defined(CHEMGEN_INTERNAL_ENERGY_EQUATION)
{scalar} temperature_from_internal_energy_py(const std::vector<{scalar}>& species, {scalar} internal_energy)
{{
    Species sp;
    std::copy(species.begin(), species.end(), sp.begin());
    return temperature(internal_energy, sp);
}}
#endif

{scalar} temperature_py(const std::vector<{scalar}>& species, {scalar} energy)
{{
    // energy is internal energy or temperature
    Species sp;
    std::copy(species.begin(), species.end(), sp.begin());
    return temperature(energy, sp);
}}

{scalar} internal_energy_volume_specific_py(const std::vector<{scalar}>& species, {scalar} temperature)
{{
    Species sp;
    std::copy(species.begin(), species.end(), sp.begin());
    return internal_energy_volume_specific(sp, temperature);
}}

std::vector<{scalar}> dtemperature_dspecies_py(const std::vector<{scalar}>& species, {scalar} temperature)
{{
    Species sp;
    std::copy(species.begin(), species.end(), sp.begin());
    auto result = dtemperature_dspecies(sp, temperature);
    return std::vector<{scalar}>(result.begin(), result.end());
}}

{remove_reactions}
PYBIND11_MODULE(chemgen, m)
{{ 
    m.def("ignore_temp_dependence", &ignore_temp_dependence_py, "ignore_temp_dependence function");
    m.def("temperature_equation", &temperature_equation_py, "temperature_equation function");
    m.def("source", &source_py, "source function");
    m.def("source_jacobian", &source_jacobian_py, "source_jacobian function");
    m.def("sdirk4", &sdirk4_py, "SDIRK 4");
    m.def("rosenbroc", &rosenbroc_py, "Rosenbroc 2");
    m.def("yass", &yass_py, "YASS");
#if defined(CHEMGEN_INTERNAL_ENERGY_EQUATION)
    m.def("temperature_from_internal_energy", &temperature_from_internal_energy_py, "temperature 4");
#endif
    m.def("temperature", &temperature_py, "temperature");
    m.def("internal_energy_volume_specific", &internal_energy_volume_specific_py, "internal_energy_volume_specific");
    m.def("dtemperature_dspecies", &dtemperature_dspecies_py, "dtemperature_dspecies");
    {remove_reactions_call}
}}

#undef {internal_energy_or_temperature}

        """.format(**vars(configuration), remove_reactions = remove_reactions_text, remove_reactions_call = remove_reactions_call_text,
                   ignore_temp_dependence=int(args.ignore_temp_dependence), temperature_equation=int(args.temperature_equation)))
    
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

