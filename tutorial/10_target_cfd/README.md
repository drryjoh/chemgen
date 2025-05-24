# ChemGen Tutorial

## Table of Contents
- [Description](#description)
- [Preparation](#preparation)
- [Code generation for OpenFOAM](#Code-generation-for-OpenFOAM)

## Description

This tutorial demonstrates how to use ChemGen to generate source code for OpenFOAM 2412.

---

## Preparation

For simplicity, we'll shorten references to ChemGen's paths. To access ChemGen in this directory, run the following command:

```bash
export PATH="$(cd ../../bin && pwd):$PATH"
```

Now, ChemGen can be executed from any directory by simply calling `chemgen.py`. Ensure that all [prerequisites are installed](../../README.md).

## Code generation for OpenFOAM

To execute this tutorial, use the following command:

```bash
chemgen.py chemkin/chem.yaml .
```

A [configuration](configuration.yaml) is provided in this tutorial with the scalar fields decalared not as double or float but as the native `scalar` used by OpenFOAM

```yaml
decorators:
  scalar: "scalar"
  scalar_function: "scalar"
  scalar_cast: "scalar"
  scalar_parameter: "const scalar&"
```

`std::array` is still used over `scalar_list` but that can also be modified.

By generating the source code, it is now compattable to be included superficially in the OpenFOAM solver ChemFOAM where `chemgenInterface.C` and `chemgenInterface.H` are used to interface OpenFOAM with chemgen directly. This abstracts away the need for OpenFOAM's `psiReactionThermo`, `BasicChemistryModel`, `reactingMixture`, `thermoPhysicsTypes`, `basicSpecieMixture`, and other mixture-based classes. After generating, the code can be compiled using OpenFOAM's `wmake` if you have a successfull OpenFOAM 2412 build. Note: unlike the existing ChemFOAM, the mechanism is not read at run time, so if you want to use a different mechanism, you will need to regenerate and re-compile.



