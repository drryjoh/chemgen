#!python3
import numpy as np
eigenvalues = np.load("eigenvalues.npy")
n_points = len(eigenvalues)
header = """
/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  6
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0.0002";
    object      chemical_time_scales;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 1 0 0 0 0];

internalField   nonuniform List<scalar> 
"""+f"{n_points}\n(\n".format(n_points)
bottom = """
);
boundaryField
{
    inlet
    {
        type            zeroGradient;
    }
    outlet
    {   
        type            zeroGradient;
    }
    bottom
    {
        type            cyclicAMI;
        value           unifom 1.0;
    }
    top
    {
        type            cyclicAMI;
        value           unifom 1.0;
    }

    }   
    frontAndBack
    {
        type            empty;
    }
}
"""
f = open("0.0002/chemical_times_scales","w")
f.write(header)
for value in eigenvalues:
    f.write(f"{value:7e}\n")
f.write(bottom)
f.close()
