# Chemgen


## Table of Contents

- [The concept of deocrators](#The-concept-of-deocrators)
- [Single point vs Double precision](#Precision)
- [Other Uses](#Other-Uses)

## The concept of deocrators

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

Provide instructions on how to use the project here. Include example commands or code snippets to guide users on how to run the project.

## Other uses


