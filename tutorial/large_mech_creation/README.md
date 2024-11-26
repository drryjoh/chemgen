# Chemgen


## Table of Contents

- [Introduction Tutorials](#Introduction Tutorials)
- [Advanced Examples](#dependencies)
- [CFD Examples](#usage)
- [License](#license)

## Introduction Tutorials 

These tutorials were created to give a general understanding of the code mechanics in ChemGen. 

1) [Simple Mechanism Creation](./simple_mech_creation/README.md)

This tutorial generates source code using a simple mechanism. The mechanism is an arbitrary single reaction and several species so that the generated source code is readable.

2) [Large Mechanism Creation](./large_mech_creation/README.md)

This tutorial generates source code using a large mechanism (FFCM2). The mechanism is large enough that the readability of all the reaction terms accrued together can be quite cumbersome, but, neverless is still usable code

3) [Mechansim Creation with Decorators](./decorators/README.md)

This tutorial gives an overview of how ChemGen can be used to target difference C++ code. For instance a simple function

```cpp
double my_function(const double& a) const {return a*a;}
```

can be represnted in chemgen as

```cpp
{device_function}
{scalar_function} my_function({scalar_parameter} a) {const_option} {{return a*a;}}
```

where the configuration.yaml file is used to specify the decorators
```yaml
decorators:
  scalar_function: double
  scalar_parameter: const double&
  const_option: const
```

this function and all other scalar functions can be changed for simgle precision  via
```yaml
decorators:
  scalar_function: float
  scalar_parameter: const float&
  const_option: const
```

to yield

```cpp
float my_function(const float& a) const {return a*a;}
```

throughout the code.
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


