# ChemGen Tutorial

## Table of Contents
- [ChemGen Tutorial](#chemgen-tutorial)
  - [Table of Contents](#table-of-contents)
  - [Description](#description)
  - [Preparation](#preparation)
  - [Example 1: Homogeneous Reactor](#example-1-homogeneous-reactor)
  - [Example 2: Order of convergence](#example-2-order-of-convergence)

## Description

This tutorial demonstrates how to use the various chemistry implicit solvers produced by ChemGen. An implicit solver solves the system:
$$
\frac{y^{n+1} - y^{n}}{\Delta t} = S\left(y^{n+1}\right)
$$

Newton's method can be used to find $y^{n+1}$ via iteration by solving

$$
\mathcal{G}\left(y_{k}^{n+1}\right)\left(y_{k+1}^{n+1}-y_{k}^{n+1}\right)
= -f\left(y_{k}^{n+1}\right),
\quad
f\left(y_{k}^{n+1}\right)
= \frac{y^{n+1} - y^{n}}{\Delta t}
- S\left(y^{n+1}\right)
$$

for $y_{k+1}^{n+1}$ , the state at the $(k+1)$ th Newton iteration.  $\mathcal{G}\left(y_{k}^{n+1}\right)$ is the Jacobian for the time integration. The Jacobian inversion is utilized incoordination with various integration schemes provided by ChemGen are listed below.


| Strategy        | Order | Stages | # Linear Solves                  |
|-----------------|:------:|:-------:|----------------------------------|
| SDIRK-2         | 2nd   | 2       | up to $n_{\mathrm{max}}$ per stage |
| SDIRK-4         | 4th   | 5       | up to $n_{\mathrm{max}}$ per stage |
| Rosenbrock      | 2nd   | 2       | one per stage                      | 
| YASS            | 1st   | 1       | one                                | 
| Backward Euler  | 1st   | 1       | up to $n_{\mathrm{max}}$           | 

More details available in our paper under review and available [here](https://arxiv.org/abs/2510.10005).

---

## Preparation

For simplicity, we'll shorten references to ChemGen's paths. To access ChemGen in this directory, run the following command:

```bash
export PATH="$(cd ../../bin && pwd):$PATH"
```

Now, ChemGen can be executed from any directory by simply calling `chemgen.py`. Ensure that all [prerequisites are installed](../../README.md).

ChemGen provides a `--custom-test` option that allows you to override the default `write_test` function to create a custom `chemgen.cpp`. This tutorial includes a `custom_test.py` file for that purpose.

Beyond example 1, pybind is utilized to demonstrate certain properties which we recommend using the [pybind tutorial](/tutorial/05_pybind/README.md).

## Example 1: Homogeneous Reactor

To execute this tutorial, use the following command:
```bash
chemgen.py FFCM2_model.yaml . --custom-test custom_test.py --compile --run
```

A [configuration](configuration.yaml) is provided in this tutorial with an extra field

```yaml
solver:
  chemistry_solver: all
  linear_solver: gmres 
  preconditioner: jacobi
```

This generates all the chemistry solvers utilized in this tutorial. The custom test generated has blocks of code for each implicit solver approach:

```cpp
    for(int i = 0; i < n_run; i++)
    {
        y = backwards_euler(y, dt);
        t = t + dt;
        be_file << t << " " << temperature(y);
        for (const auto& val : get_species(y)) be_file << " " << val;
        be_file << "\n";
    }
```



## Example 2: Order of convergence 

