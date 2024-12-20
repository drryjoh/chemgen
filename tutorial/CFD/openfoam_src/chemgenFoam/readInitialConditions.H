    word constProp(initialConditions.get<word>("constantProperty"));
    if (constProp != "volume")
    {
        FatalError << "in initialConditions, unknown constantProperty type "
            << constProp << nl << " Valid types are: volume."
            << abort(FatalError);
    }

    word fractionBasis(initialConditions.get<word>("fractionBasis"));
    if (fractionBasis != "mass" && fractionBasis != "mole")
    {
        FatalError << "in initialConditions, unknown fractionBasis type " << nl
            << "Valid types are: mass or mole."
            << fractionBasis << abort(FatalError);
    }


    
    scalarList Y0(n_species, Zero);
    scalarList X0(n_species, Zero);
    scalar p0 = initialConditions.get<scalar>("p");
    scalar T0 = initialConditions.get<scalar>("T");
    ChemgenInterface chemgen_interface(mesh);

    Species MW = molecular_weights();

    dictionary fractions(initialConditions.subDict("fractions"));
    if (fractionBasis == "mole")
    {
        forAll(Y0, i)
        {
            fractions.readIfPresent(chemgen_interface.species_name(i), X0[i]);
        }

        scalar mw = 0.0;
        const scalar mTot = sum(X0);
        forAll(Y0, i)
        {
            X0[i] /= mTot;
            mw += MW[i]*X0[i];
        }

        forAll(Y0, i)
        {
            Y0[i] = X0[i]*MW[i]/mw;
        }
    }
    else  // mass fraction
    {
        forAll(Y0, i)
        {
            fractions.readIfPresent(chemgen_interface.species_name(i), Y0[i]);
        }

        scalar invW = 0.0;
        const scalar mTot = sum(Y0);
        forAll(Y0, i)
        {
            Y0[i] /= mTot;
            invW += Y0[i]/MW[i];
        }
        const scalar mw = 1.0/invW;

        forAll(Y0, i)
        {
            X0[i] = Y0[i]*mw/MW[i];
        }
    }

    Species Y0_species  = {scalar(0.0)};
    forAll(Y0, i)
    {
        Y0_species[i] = Y0[i];
    }

    CustomChemistryModel chemgen_chemistry = CustomChemistryModel::TPY(T0, p0, Y0_species, mesh);

    
    volScalarField& T = chemgen_chemistry.T();
    volScalarField& p = chemgen_chemistry.p();
    volScalarField& rho = chemgen_chemistry.rho();
    volScalarField& internal_energy = chemgen_chemistry.internal_energy();
    PtrList<volScalarField>& concentrations = chemgen_chemistry.C();

    PtrList<volScalarField> Y = chemgen_chemistry.Y();

    scalar integratedHeat = 0.0;

    Info << constProp << " will be held constant." << nl
        << " p   = " << p[0] << " [Pa]" << nl
        << " T   = " << T[0] << " [K] " << nl
        << " rho   = " << rho[0] << " [kg/m^3] " << nl
        << " rhoU   = " << internal_energy[0] << " [J/m^3] " << nl
        << endl;

    forAll(Y, i)
    {
        Info<< Y[i].name() <<": "<<Y[i][0]<<endl;
    }

    forAll(concentrations, i)
    {
        Info<< concentrations[i].name() <<": "<<concentrations[i][0]<<endl;
    }
    Info << "T from rhoU" <<""<<endl;
    chemgen_chemistry.temperature_from_internal_energy();
    Info << " T   = " << T[0] << " [K] " << nl;

