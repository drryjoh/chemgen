# ChemGen Tutorial

## Table of Contents
- [Description](#description)
- [Preparation](#preparation)
- [Building Kokkos](#Building-kokkos)

## Description

This tutorial demonstrates how to use ChemGen to generate source code using the Kokkos library for fine grain parellism. Specifically the instructions here are for a mac book with an i9 core (8 cores) using POSIX for threading. Over time we will update this tuturial to cover other build configurations.


## Preparation

For simplicity, we'll shorten references to ChemGen's paths. To access ChemGen in this directory, run the following command:

```bash
export PATH="$(cd ../../bin && pwd):$PATH"
```

Now, ChemGen can be executed from any directory by simply calling `chemgen.py`. Ensure that all [prerequisites are installed](../../README.md).

ChemGen provides a `--custom-test` option that allows you to override the default `write_test` function to create a custom `chemgen.cpp`. This tutorial includes a `custom_test.py` file for that purpose.


## Building kokkos
In the root folder for chemgen pull all submodules

```
cd /path/to/chemgen/
git submodule init
git submodule update
```
This will pull the Kokkos library as third party

To compile Kokkos for this tutorial we provided a build and install directory, enter the build directory and build via

```terminal
cd build_kokkos
cmake ../../../third_party/kokkos -DCMAKE_INSTALL_PREFIX=../install_kokkos -DKokkos_ENABLE_THREADS=ON -DKokkos_ENABLE_SERIAL=ON
make -j 8 install
```

Once Kokkos is built we can then generate chemistry code that uses it. 

## Kokkos demonstration case

Here we evaluate the source term for 10,000 different chemical states using the hydrogen mechanism extracted from FFCM2. The command

```
chemgen.py ffcm2_h2.yaml . --custom-test custom_test.py --n-points 10000 
```

creates the custom test case in chemgen.cpp that utilizes Kokkos. Key source code is the creation of views for the chemical states and the parallel for which uses the compiled fine grain parallelism (in this example POSIX). The views are created as their own types (where n_points is changed in the generation process)

```
const int n_points = 10000;
using PointScalar = Kokkos::View<double*>;
using PointSpecies = Kokkos::View<double**>;
```

Inside of main these views are initialized

```        
        PointSpecies concentration_tests_device("concentration_tests", n_points, n_species);
        PointSpecies cantera_sources_device("cantera_sources", n_points, n_species);
        PointScalar temperatures_device("temperatures", n_points);
```
and a parallel for loop is used to evaluate the source term
```

    Kokkos::parallel_for("ProcessPoints", n_points, KOKKOS_LAMBDA(int i) 
    {
        // Extract a subview (row) of species concentrations
        auto conc_subview = Kokkos::subview(concentration_tests_device, i, Kokkos::ALL());

        // Copy subview data to a local array
        Species conc_array;  
        for (int j = 0; j < n_species; j++)
        {
            conc_array[j] = conc_subview(j);
        }

        // Call source() with the extracted data
        Species result = source(conc_array, temperatures_device(i));

        // Add other computations if needed...
    });
```

To compile with posix and kokkos we used the command:

```
g++ -o parallel_kokkos test_kokkos.cpp -I./install_kokkos/include -L./install_kokkos/lib -lkokkoscore -std=c++17 -pthread
```

Once compiled the number of available cores can be found on a mac via 

```
sysctl -n hw.logicalcpu
```

Then the effective speed up of fine grain parallelism can be seen by performing the following:
```terminal
% export KOKKOS_NUM_THREADS=1
% ./parallel_kokkos          
*** ChemGen ***
Total execution time: 0.0368081 seconds
% export KOKKOS_NUM_THREADS=2
% ./parallel_kokkos          
*** ChemGen ***
Total execution time: 0.0191024 seconds
% export KOKKOS_NUM_THREADS=4
% ./parallel_kokkos          
*** ChemGen ***
Total execution time: 0.0109041 seconds
% export KOKKOS_NUM_THREADS=8
% ./parallel_kokkos          
*** ChemGen ***
Total execution time: 0.00625924 seconds
```
