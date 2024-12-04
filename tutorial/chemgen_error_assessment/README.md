# Chemgen


## Table of Contents

- [Generating Custom Test](#Generating Custom Test)
- [Generating Error Data](#Generating Data)
- [Post Processing Data](#Post Processing)

## Description

Here we will use chemgen to generate the source code for evaluating the species production source term. We then evaluate the accuracy of this error by comparing to cantera pre-calculated source terms. 

## Preparation

In all tutorials we will shorten the use of ChemGen's paths to the repo. To access ChemGen in this folder run the following command:

```bash
export PATH="$(cd ../../bin && pwd):$PATH"
```

Now, ChemGen can be run from any directory by simply specifying `chemgen.py`. Make sure all [prequisites are installed](../../README.md).

ChemGen has an option "--custom-test" where a python function, `write_test` can be overwritten to create a custom `chemgen.cpp`. We've included `custom_test.py` in this directory.


## Generating Data

To run this tutorial the following command should be used
```terminal
chemgen.py FFCM2_model .  --custom-test custom_test.py --n-points-test 1000 --compile
```

The typical fomrmat for ChemGen execution is

```terminal
chemgen.py [path/to/kinetics/file] [path/to/generated/code]
```

the first item in the command references the chemgen python source which executes the code generation, the second command is the name of a chemical kinetic model (in this case FFCM2_model), the . directs the target directory where to generate source. Additional options are used in this tutorial denoted by the `--` command. The `--custom-test` utilizes a custom test generation covered in detailed in the next session. The `--n-points-test` passes in a random number of points to test, in this case 1000, and `--compile` compiles the test cases after generating the source code. The default compilation is

```
clang++ -std=c++17 -O2 -o ./bin/chemgen src/chemgen.cpp
```

However, cmake will be used in future tutorials and can also be generationed using the `--cmake` command instead.


## Custom test

Included in the directory is a python file used to overwrite the default custom_test generated. The custom test generates `n_points` of number of chemical states and randomly creates the chemical makeup (concentrations) and temperature. The source term is then evaluated for each point, which creates `n_species` amount of data per point. Each species specific source term is compared to cantera and an L2-norm is generated per point

```math
l_2 = \sqrt{\sum_{i=1}^{n_s}\frac{1}{n_s}((S(y_i)_{cg}-S(y_i)_{ct})/S(y_i)_{ct})^2}
```

Where `$ct$` is the cantera solution and `$cg$` is the chemgen solution. 
