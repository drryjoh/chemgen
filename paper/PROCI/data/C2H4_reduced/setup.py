#!python3
from setuptools import setup, Extension
import sys
import pybind11

formulations  = ["chemgen_conservative","chemgen_temperature","chemgen_dtemperature_ignore", "chemgen_conservative_ignore_species", "chemgen_temperature_ignore_species"]
for formulation in formulations:
    # Open, replace, and overwrite the file
    file_path = f"{formulation}/src/chemgen_pybind.cpp"   # change to your filename

    with open(file_path, "r") as f:
        content = f.read()

    # Replace the line
    old = "PYBIND11_MODULE(chemgen, m)"
    new = f'PYBIND11_MODULE({formulation}, m)'

    content = content.replace(old, new)

    with open(file_path, "w") as f:
        f.write(content)

    ext_modules = [
        Extension(
            formulation,
            sources=[f'{formulation}/src/chemgen_pybind.cpp'],
            include_dirs=[pybind11.get_include()],
            language='c++',
            extra_compile_args=['-std=c++17']
        )
    ]

    setup(
        name=formulation,
        version='1.0',
        author='Ryan F. Johnson',
        description='Pybind11 bindings for C++ ChemGen',
        ext_modules=ext_modules,
    )
            