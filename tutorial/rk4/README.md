# ChemGen Tutorial

## Table of Contents
- [Description](#description)
- [Preparation](#preparation)
- [Generating Data for Error Assessmnet](#generating-data-for-error-assessment)
- [ChemGen Execution Format](#chemgen-execution-format)
- [Generating Custom Tests](#generating-custom-tests)
- [Post-Processing Data](#post-processing-data)

## Description

This tutorial demonstrates how to use ChemGen to generate source code for evaluating species production source terms and simulate a homogeneous reactor
```math
\frac{dy_c}{dt} = S(y_c) 
```

With an initial state $`y_c(t=0)=y_{c,0}`$

---

## Preparation

For simplicity, we'll shorten references to ChemGen's paths. To access ChemGen in this directory, run the following command:

```bash
export PATH="$(cd ../../bin && pwd):$PATH"
```

Now, ChemGen can be executed from any directory by simply calling `chemgen.py`. Ensure that all [prerequisites are installed](../../README.md).

ChemGen provides a `--custom-test` option that allows you to override the default `write_test` function to create a custom `chemgen.cpp`. This tutorial includes a `custom_test.py` file for that purpose.


## Generating Data for Error Assessment

To execute this tutorial, use the following command:

```bash
chemgen.py ffcm2_h2.yaml . --custom-test custom_test.py --compile >> chem_out.txt
```

A [configuration](configuration.yaml) is provided in this tutorial with an extra field

```yaml
solver:
  chemistry_solver: rk4
```

which compiles in the Runge Kutta fourth oder solver found in [rk4.in.h](../../src/solvers/rk4.h.in). The solver is used in the generated chemgen.cpp file from [write_test](write_test.py):

```cpp
   for(int i = 0; i < 40000; i++)
    {
        y = rk4(y, dt);
        t = t + dt;
        std::cout <<t<<" "<<temperature(y) <<" "<< get_species(y) << std::endl;
    }
```

Once compiled and run the `chem_out.txt` can be used to compare to cantera's homogeneous reactor with the script [post_ct](post_ct.py)


```
./post_ct.py
```
.
![Error as a plotted using a box and whisker](bw.png)
![Error distribution](hist.png)
