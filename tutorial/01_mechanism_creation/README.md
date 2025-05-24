# Chemgen


## Table of Contents

- [Add ChemGen to Path](#Add ChemGen to path)
- [Run ChemGen](#Run ChemGen)

## Add ChemGen to Path 

In all tutorials we will shorten the use of ChemGen's paths to the repo. To access ChemGen in this folder run the following command:

```bash
export PATH="$(cd ../../../bin && pwd):$PATH"
```

Now, ChemGen can be run from any directory by simply specifying `chemgen.py`

## Run ChemGen

Chemgen will automatically look for cantera a .yaml files located in either /chemical_mechansims/ or local folder. To run Chemgen you need to provide parameters, the chemical mechanism and the source directory where code will be generated.

```bash
chemgen.py <source_mechanism> <directory>
```

we have provided  a simple test case that will evaluate various ChemGen subrouties after compiling a chemgen.cpp file. It is important to note that not all generated code requires this test, we just included the interface to test a variety of surboutines and advanced implementations.

For the first test run

```bash
chemgen.py simple_mech .
```

Both `simple_mech.yaml` and `simple_mech` can be used. The yaml extension is automatically searched for. For this test case, `simple_mech.yaml` is included in `/chemical_mechanisms/`. 

There will be a newly generated directory, `src/` which contains
```
array_handling.h
exp_gen.h
generated_constants.h
reactions.h
third_body.h
arrhenius.h
falloff_lindemann.h
multiply_divide.h
source.h
types_inl.h
falloff_sri.h
pow_gen.h
thermally_perfect.h
constants.h
falloff_troe.h
pressure_dependent_arrhenius.h
thermotransport_fits.h
```
Exploring the generated code one can see all functionality to calculate the source term for `simple_mech` using c++ has been generated. For instance the reaction in the yaml file

```yaml
- equation: 2 H2 + O2 <=> 2 H2O  # Reaction 1
  rate-constant: {A: 1.0399e+14, b: 0.0, Ea: 1.531e+04}
  duplicate: true
```

has been generated as
```cpp
double call_forward_reaction_0(const double& temperature, const double& log_temperature)  { return arrhenius(double(103990000.00000003), double(0.0), double(64057040.0), temperature, log_temperature);}
```
in `reactions.h` which is then used in `source.h`.

running
```bash
chemgen.py simple_mech . --compile
```

compiles and runs a new file, `chemgen.cpp`, which includes all the header files and performs a test with output
```bash
*** ChemGen ***
temperature: 1800
temperature_ for 0 iterations: 0.722222
temperature_ for 1 iterations: -0.101355
temperature_ for 2 iterations: -0.000886011
temperature_ for 3 iterations: -7.4822e-08
temperature_ for 4 iterations: -3.78956e-16
temperature_ for 5 iterations: 0
temperature_ for 6 iterations: 0
temperature_ for 7 iterations: 0
temperature_ for 8 iterations: 0
temperature_ for 9 iterations: 0
Source test result:  [ 0 -0.00727977 0 -0.00363989 0 0.00727977 0 0 0 0 0 0 0 0 ]
Cantera test result: 0.0 -0.007294716394335021 0.0 -0.0036473581971675113 0.0 0.007294716394335021 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
ChemGen internal energy: 229952
Cantera internal energy: 229948.9697167456
Chemgen species cps: [ 20621.2 16582.3 1302.19 1169.21 2001.95 2780.05 1652.17 2129.14 5193.16 1271.49 1739.24 1283.2 1358.35 2074.67 ]
Cantera species cps: 20621.18704899117 16616.78802469927 1301.6863166535868 1167.0799719707766 2004.0508234338902 2781.953463285981 1650.4585771578277 2127.8710008356074 5193.1609851249505 1270.558041870319 1738.5816213185194 1281.8044880331277 1357.1758540814728 2075.9084831417763
Chemgen species enthalpies: [ 2.47237e+08 2.29062e+07 1.75468e+07 1.61511e+06 4.9521e+06 -9.93407e+06 2.5121e+06 -1.28094e+06 7.79935e+06 1.74668e+06 6.22683e+07 -2.17974e+06 -7.13874e+06 2.81289e+07 ]
Cantera species enthalpies: 247236979.76114592 22924287.42497859 17546354.348052036 1614086.106703522 4952835.895921169 -9933592.095499126 2510821.7456610794 -1281946.3720538646 7799348.825509908 1746606.8628090904 62268039.86086621 -2180056.7442304264 -7139104.986552956 28129761.151496165
Chemgen species internal energies: [ 2.3239e+08 1.54825e+07 1.66113e+07 1.1474e+06 4.07211e+06 -1.07648e+07 2.05866e+06 -1.72093e+06 4.06027e+06 1.21244e+06 6.10222e+07 -2.71405e+06 -7.47881e+06 2.72489e+07 ]
Cantera species internal energies: 232389725.08587226 15500660.08734177 16610918.838790463 1146368.3520727376 4072843.38032901 -10764346.062342083 2057388.0453436875 -1721942.6298499438 4060272.9162199423 1212372.811524881 61022012.659744255 -2714367.08741771 -7479172.534388077 27249768.635904007
Chemgen species internal entropies: [ 150883 91638.6 12445.2 8278.82 14038.8 14411.3 9347.74 9911.46 40854.8 8865.15 16277.3 9102.87 6880.28 13803.5 ]
Cantera species internal entropies: 150882.52021782557 91714.6332852643 12441.629169912665 8274.873072719249 14037.717709743141 14408.506261912842 9344.433258665218 9909.851918999964 40854.761228438154 8862.383000673015 16276.6431797554 9099.375704537315 6882.076039934865 13803.098376233454
Chemgen species gibbs energy: [ -2.45464e+07 -2.86362e+08 -7.76736e+07 -4.25146e+08 -3.45556e+08 -6.46295e+08 -4.72431e+08 -6.50385e+08 -2.63128e+08 -3.98105e+08 3.95992e+08 -5.20011e+08 -8.59171e+08 5.58154e+07 ]
Cantera species gibbs energy: -24546369.08398764 -286598697.81681025 -77572401.9464944 -424955372.2032677 -345498157.07935154 -646178294.1554644 -472288072.9066179 -650336789.6101968 -263127938.99676064 -397957990.63080305 396003656.55118823 -519835713.6772565 -859356783.3478973 55854118.55121097
Pressure: 101325
Temperature Monomial at 300           : [ 1 300 90000 2.7e+07 8.1e+09 2.43e+12 7.29e+14 2.187e+17 ]
Temperature Energy Monomial at 300           : [ 1 300 90000 2.7e+07 8.1e+09 2.43e+12 7.29e+14 2.187e+17 6.561e+19 ]
Temperature Entropy Monomial at 300           : [ 1 300 90000 2.7e+07 8.1e+09 2.43e+12 7.29e+14 2.187e+17 5.70378 ]
Temperature Gibbs Monomial at 300           : [ 1 300 90000 2.7e+07 8.1e+09 2.43e+12 7.29e+14 2.187e+17 6.561e+19 1711.13 ]
```

This tests the source code and other thermodyanic properties, consult the `chemgen.cpp` for the source code implementation.

running
```bash
chemgen.py simple_mech . --cmake
```

will make a CMakeLists.txt file for user-imposed compilation strategies. 

```bash
mkdir build
cd build
ccmake ../
make
./bin/chemgen
```

Previously, running
```bash
chemgen.py simple_mech . --compile
```

Generates source code and compile with a hard coded command. Eventually these hard coded commands will be replaced with generating cmake on the fly and compiling. If your machine has specific compilation setting the `cmake` approach is probably best.


## Run Larger Mechanism

The previously generated test, `chemgen.cpp`, utilized `test_configuration.yaml` in this directory which includes the following configuration:

```yaml
test_conditions:
  temperature: 1800
  pressure: 101325.0
  species:
    - name: O2
      MoleFraction: 0.2
    - name: N2
      MoleFraction: 0.6
    - name: H2
      MoleFraction: 0.2
```
to create a test state that generated the output above. We can generate a larger mechanism and perform test for with more species by involving other species in the configuration fie:
```yaml
test_conditions:
  temperature: 1800
  pressure: 101325.0
  species:
    - name: O2
      MoleFraction: 0.2
    - name: N2
      MoleFraction: 0.6
    - name: H2
      MoleFraction: 0.1
    - name: C2H4
      MoleFraction: 0.1
```

or we could generate a random state
```yaml
test_conditions:
  random: On
  temperature: 1800
  pressure: 101325.0
  species:
    - name: O2
      MoleFraction: 0.2
    - name: N2
      MoleFraction: 0.6
    - name: H2
      MoleFraction: 0.1
    - name: C2H4
      MoleFraction: 0.1
```
Note: `random: On` supersedes any other settings in the test configuration file. We can perform this larger test with FFCM2_model.yaml. The FFCM2_model can be download in the chemical mechanisms folder
```bash
cd /chemgen/chemical_mechansims/
./download.sh
```

Then you can run chemgen in any run directory:

```bash
chemgen.py FFCM2_model.yaml . --compile
```

For the explicit `C2H4` addition, the source term is easy to parse, but the random data is large and cumbersome. In further tutorials we perform random scattering of the source term to compare its `$L_2$`-norm and compute ChemGen's accuracy compared to cantera.
