#ifndef CUSTOMCHEMISTRYMODEL_H
#define CUSTOMCHEMISTRYMODEL_H

#include "chemgenInterface.H" // Include your SourceChemgen struct

namespace Foam
{

    class CustomChemistryModel
    {
        private:
            mutable ChemgenInterface ChemgenInterface_; // Instance of chemgen interface
            volScalarField T_;             // Temperature field
            volScalarField p_;             // Temperature field
            PtrList<volScalarField> Y_;
            volScalarField rho_;           // Density field
            volScalarField internal_energy_;           // Density field
            PtrList<volScalarField> C_;          // Concentrations mass fractions


        public:
        CustomChemistryModel(
            const volScalarField& T,
            const volScalarField& p,
            PtrList<volScalarField>& Y,
            const fvMesh& mesh
        );
        // Destructor
        ~CustomChemistryModel();
        // Factory methods
        static CustomChemistryModel TPY(
            const volScalarField& T,
            const volScalarField& p,
            PtrList<volScalarField>& Y,
            const fvMesh& mesh
        );

        static CustomChemistryModel TPY(
            scalar T0,
            scalar p0,
            const Species& Y0,
            const fvMesh& mesh
        );



            // Methods
            scalar solve(const scalar deltaT);
            volScalarField heat_release(const PtrList<volScalarField>& sourceField) const;
            volScalarField heat_release() const;
            PtrList<volScalarField> sourceField() const;
            PtrList<volScalarField>& Y();
            PtrList<volScalarField>& C();
            volScalarField& T();
            volScalarField& p();
            volScalarField& internal_energy();
            volScalarField& rho();
            void temperature_from_internal_energy();
            void update_chemical_state_from_concentrations_and_temperature();
    };

} // namespace Foam

#include "chemgen.C"

#endif
