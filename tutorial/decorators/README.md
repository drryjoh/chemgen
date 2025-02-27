# Chemgen


## Table of Contents

- [The concept of decorators](#The-concept-of-deocrators)
- [Single point vs Double precision](#Precision)
- [Other Uses](#Other-Uses)

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
Using the above example we have provided a [configuration_float](../bin/configuration_float.yaml) [configuration_double](../bin/configuration_double.yaml). By performing the following

```
cp configuration_float.yaml configuration.yaml
chemgen.py one_reaction . --custom-test custom_test.py --compile
cp configuration_double.yaml configuration.yaml
chemgen.py one_reaction . --custom-test custom_test.py --compile
```

The resulting output demonstrating the difference in precision is found
```
Source test result for scalar float:   [ -0.0027410432230681 -0.0013705216115341 0.0027410432230681 0.0000000000000000 ]
Source test result for scalar double:  [ -0.0027410426276326 -0.0013705213138163 0.0027410426276326 0.0000000000000000 ]
```

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
```
Which definitely requires more nuances, but can be used to enable Kokkos


