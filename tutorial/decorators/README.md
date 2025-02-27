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


