Info << "For Cantera:" << endl;
Info << T[0]<<", " << p[0]<<", [";
forAll(Y, specieI)
{
    Info << Y[specieI][0]<< "," ;
}
dtChem = chemistry.solve(runTime.deltaTValue());
forAll(Y, i)
{
    Info << "Y" << i <<": "<<Y[i][0] <<endl;
}
scalar Qdot = chemistry.Qdot()()[0]/rho[0];
integratedHeat += Qdot*runTime.deltaTValue();
