    runTime.write();
    Info << "dtChem" << dtChem <<endl;
    forAll(Y, specieI)
    {
        auto S1 = customChemistryModel.sourceField();
        auto S2 = chemistry.RR(specieI);
        Info << Y[specieI].name() << ": "<<S1[specieI][0] <<", "<< S2[0] <<", "<< S2[0]/(S1[specieI][0]+1e-6)<<endl;
    }
    scalar Qdot_OF = chemistry.Qdot()()[0]/rho[0];
    scalar Qdot_cg = customChemistryModel.heat_release()[0]/rho[0];

    Info << "Hdot OF: " << Qdot_OF << " Hdot CG: " << Qdot_cg <<endl;
    Info << "For Cantera:" << endl;
    Info << T[0]<<", " << p[0]<<", [";
    forAll(Y, specieI)
    {
        Info << Y[specieI][0]<< "," ;
    }
    Info <<endl;
 
    post<< runTime.value() << token::TAB << thermo.T()[0] << token::TAB
        << thermo.p()[0] << endl;

