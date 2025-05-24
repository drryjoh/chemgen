#include "chemgen.H"

namespace Foam
{

// Private constructor
CustomChemistryModel::CustomChemistryModel(
    const volScalarField& T,
    const volScalarField& p,
    PtrList<volScalarField>& Y,
    const fvMesh& mesh
)
: ChemgenInterface_(mesh),    // 1. Initialize ChemgenInterface
  T_(T),                      // 2. Initialize T_
  p_(p),                      // 3. Initialize p_
  Y_(Y),                      // 4. Bind reference Y_
  rho_(                       // 5. Initialize rho_
      IOobject(
          "rho",
          mesh.time().timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar("zero", dimDensity, 0.0)
  ),
  internal_energy_(           // 6. Initialize internal_energy_
      IOobject(
          "internal_energy",
          mesh.time().timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar("zero", dimEnergy / dimVolume, 0.0)
  ),
  C_(Y.size())                // 7. Initialize C_ with the same size as Y_
{
    rho_ = ChemgenInterface_.density(Y_, p_, T_);
    C_ = ChemgenInterface_.concentrations(Y_, rho_, T_);
    internal_energy_ = ChemgenInterface_.internal_energy(Y_, rho_, T_);
}


// Destructor
CustomChemistryModel::~CustomChemistryModel() {}

// Factory method using TPY fields
CustomChemistryModel 
CustomChemistryModel::TPY(
    const volScalarField& T,
    const volScalarField& p,
    PtrList<volScalarField>& Y,
    const fvMesh& mesh
)
{
    
    return CustomChemistryModel(T, p, Y, mesh);
}


// Factory method using TPY fields
CustomChemistryModel 
CustomChemistryModel::TPY(
    const scalar T0,
    const scalar p0,
    const Species& Y0,
    const fvMesh& mesh
)
{
    volScalarField p
    (
        IOobject
        (
            "p",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar("p", dimPressure, 0)
    );

    volScalarField T
    (
        IOobject
        (
            "T",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar("T", dimTemperature, T0)
    );

    
    PtrList<volScalarField> Y(n_species); // Create PtrList with the same size as Y

    forAll(Y, i)
    {
        Y.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Y" + ChemgenInterface::species_name(i),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimless, Y0[i]) // Initialize with dimensions
            )
        );
    }

    return TPY(T, p, Y, mesh);
}

PtrList<volScalarField> CustomChemistryModel::sourceField() const
{
    return ChemgenInterface_.source_field(Y_, rho_, T_);
}

scalar CustomChemistryModel::solve(const scalar deltaT)
{
    // Generate source terms
    PtrList<volScalarField> sourceField = this->sourceField();

    // Update species fields
    forAll(C_, specieI)
    {
        volScalarField& Ci = C_[specieI];
        forAll(Ci, celli)
        {
            Ci += sourceField[specieI][celli] * deltaT;
        }
    }
    temperature_from_internal_energy();

    return deltaT;
}

volScalarField CustomChemistryModel::heat_release(
    const PtrList<volScalarField>& sourceField
) const
{
    volScalarField Qdot(
        IOobject(
            "Qdot",
            T_.time().timeName(),
            T_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        T_.mesh(),
        dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
    );
    
    forAll(Qdot, celli)
    {
        auto hc = species_enthalpy_mass_specific(298.0);
        forAll(sourceField, i)    
        {
            Qdot[celli] -= hc[i] * sourceField[i][celli]; // Adjust multiplier as needed
        }
    }

    Qdot.correctBoundaryConditions();
    return Qdot;
}

volScalarField CustomChemistryModel::heat_release() const
{
    return heat_release(this->sourceField());
}
// Return species mass fractions (Y_)
PtrList<volScalarField>& CustomChemistryModel::Y()
{
    return Y_;
}

// Return concentrations 
PtrList<volScalarField>& CustomChemistryModel::C()
{
    return C_;
}

// Return temperature field (T_)
volScalarField& CustomChemistryModel::T()
{
    return T_;
}

// Return pressure field (p_)
volScalarField& CustomChemistryModel::p()
{
    return p_;
}

// Return density field (rho_)
volScalarField& CustomChemistryModel::rho()
{
    return rho_;
}

// Return internal energy field (internal_energy_)
volScalarField& CustomChemistryModel::internal_energy()
{
    return internal_energy_;
}

void CustomChemistryModel::temperature_from_internal_energy()
{
    forAll(C_, specieI)
    {
        Info << C_[specieI][0] << endl;
    }
    ChemgenInterface_.temperature_from_internal_energy(C_, internal_energy_, T_);
}

void CustomChemistryModel::update_chemical_state_from_concentrations_and_temperature()
{
    ChemgenInterface_.density_from_concentrations(rho_, C_);   
    ChemgenInterface_.pressure_from_concentrations_temperature(p_, C_, T_);
    ChemgenInterface_.mass_fractions_from_concentrations_density(Y_, C_, rho_);
}


} // namespace Foam
