from setuptools import setup, Extension
import sys
import pybind11

ext_modules = [
    Extension(
        'chemwrapper',
        sources=['chem.cpp', 'bindings.cpp'],
        include_dirs=[pybind11.get_include()],
        language='c++',
        extra_compile_args=['-std=c++17']
    )
]

setup(
    name='chemwrapper',
    version='0.1',
    author='Your Name',
    description='Pybind11 bindings for C++ chemistry code',
    ext_modules=ext_modules,
)

