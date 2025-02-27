# Chemgen


## Table of Contents

- [Introduction Tutorials](#Introduction-Tutorials)
- [Advanced Examples](#dependencies)
- [CFD Examples](#usage)
- [License](#license)

## Introduction Tutorials 

These tutorials were created to give a general understanding of the code mechanics in ChemGen. 

1) [Mechanism Creation](./mechanism_creation/README.md)

This tutorial generates source code using a the most basic ChemGen features. The mechanism is an arbitrary single reaction and several species so that the generated source code is readable.

2) [ChemGen Error Assessment](./chemgen_error_assessment/README.md) 

This tutorial creates a custom chemgen.cpp function that tests random chemical states and generates L2-norms for the specificied mechanism. We then plot their errors in various ways.

3) [Mechansim Creation with Decorators](./decorators/README.md)

This tutorial gives an overview of how ChemGen can be used to target difference C++ code using the decorators. Decorators can enable GPU device execuction via kokkos or simply change how variables are passed into functions. 

## Dependencies

The project requires the following Python packages to run:

- **[Cantera](https://cantera.org)**: Cantera is an open-source suite of tools for problems involving chemical kinetics, thermodynamics, and transport processes.
  
  Install it using:

  ```bash
  python3 -m pip install cantera
  ```

- **[PyYAML](https://pyyaml.org/)**: A Python library for parsing and emitting YAML.
  
  Install it using:

  ```bash
  python3 -m pip install pyyaml
  ```

Ensure you are using Python 3 and `pip` is invoked through `python3 -m pip` for consistent environment management.

## Usage

Provide instructions on how to use the project here. Include example commands or code snippets to guide users on how to run the project.

## License


