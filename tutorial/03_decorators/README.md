# Chemgen


## Table of Contents



- [Chemgen](#chemgen)
  - [Table of Contents](#table-of-contents)
  - [Preparation](#preparation)
  - [The concept of decorators](#the-concept-of-decorators)
  - [Precision](#precision)
  - [Other uses](#other-uses)
    - [Device functions](#device-functions)

## Preparation

For simplicity, we'll shorten references to ChemGen's paths. To access ChemGen in this directory, run the following command:

```bash
export PATH="$(cd ../../bin && pwd):$PATH"
```

Now, ChemGen can be executed from any directory by simply calling `chemgen.py`. Ensure that all [prerequisites are installed](../../README.md).

ChemGen provides a `--custom-test` option that allows you to override the default `write_test` function to create a custom `chemgen.cpp`. This tutorial includes a `custom_test.py` file for that purpose.

## The concept of decorators

Consider the target function that we would like generate as function in `myfunc.h` file

```cpp
double my_function(const double& a) const {return a*a;}
```

In chemgen,  a `myfunc.h.in` file would be created that looks like

```cpp
{device_function}
{scalar_function} my_function({scalar_parameter} a) {const_option} {{return a*a;}}
```

where the configuration.yaml (located in [bin](../bin/configuration.yaml) or locally) file is used to specify the decorators
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

## Precision
Using the above example we have provided a [configuration_float](configuration_float.yaml) and [configuration_double](configuration_double.yaml) file, along with a [test_configuration.yaml](test_configuration.yaml) state for this directory's `one_reaction` mechanism. This tutorial has no `custom_test.py`, so drop `--custom-test` and use the default test writer with `--run-tests` (chemgen's actual flag -- there is no `--run`):

```bash
cp configuration_float.yaml configuration.yaml
chemgen.py one_reaction . --compile --run-tests

cp configuration_double.yaml configuration.yaml
chemgen.py one_reaction . --compile --run-tests
```

Both runs print a `Source test result:` line. `std::cout`'s default precision (6 significant figures) hides the float/double difference at a glance -- to actually see it, increase the printed precision (e.g. `std::cout << std::setprecision(15)`) or inspect `src/reactions.h`, where the scalar type substitution is visible directly in the generated signatures (see below). `tests/test_tutorial_03_decorators.py` exercises both configs and asserts they each still agree with Cantera within their precision's tolerance.

## Other uses

The decorators:
```yaml
  scalar_parameter: "const double&"
```

Are meant to give ability to pass by reference and call functions in a certain way. By removing the const-reference,
```yaml
  scalar_parameter: "double"
```
the generated functions change from
```cpp
double call_forward_reaction_0(const double& temperature, const double& log_temperature)  { return arrhenius(double(103990000.00000003), double(0.0), double(64057040.0), temperature, log_temperature);}
```
to

```cpp
double call_forward_reaction_0(double temperature, double log_temperature)  { return arrhenius(double(103990000.00000003), double(0.0), double(64057040.0), temperature, log_temperature);}
```

This can be done using the [configuration_point](configuration_double.yaml) file

In addition if these functions are to go into a struct and you desire them to be const you can change
```yaml
const_option: ""
```
to

```yaml
const_option: "const"
```

```cpp
double call_forward_reaction_0(const double& temperature, const double& log_temperature) const { return arrhenius(double(103990000.00000003), double(0.0), double(64057040.0), temperature, log_temperature);}
```

### Device functions

A later tutorial will demonstrate this in better detail, but

```yaml
 device_option: ""
 species_typedef: "std::array<double, n_species>"
```
can be changed to 
```yaml
 device_option: "KOKKOS_INLINE_FUNCTION"
 species_typedef: "Kokkos::View<double[n_species]>"
```

Which yields
```
using Species = Kokkos::View<double[n_species]>;
KOKKOS_INLINE_FUNCTION
double call_forward_reaction_0(double temperature, double log_temperature)  { return arrhenius(double(103990000.00000003), double(0.0), double(64057040.0), temperature, log_temperature);}
```

Which definitely requires more nuances, but can be used to enable Kokkos
