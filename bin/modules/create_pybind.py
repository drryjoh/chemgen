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

def create_pybind(gas, headers, configuration, destination_folder):
    includes = headers
    bindings = []
    bindings_file = Path(destination_folder) / Path("chemgen_bindings.cpp")     # generated binding file\
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

        f.write('namespace py = pybind11;\n\n')
        f.write(f'PYBIND11_MODULE(chemgen, m)\n')
        f.write("{")
        f.write('\n'.join(bindings))
        f.write('\n}\n')
    
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
    author='Ryan Johnson',
    description='Pybind11 bindings for C++ ChemGen',
    ext_modules=ext_modules,
)
        """)
    

    print(f"[+] Generated {bindings_file} with {len(bindings)} bindings.")

