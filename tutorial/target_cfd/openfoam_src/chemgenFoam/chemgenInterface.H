#include "src/types_inl.h"
#include "src/multiply_divide.h"
#include "src/pow_gen.h"
#include "src/exp_gen.h"
#include "src/array_handling.h"
#include "src/constants.h"
#include "src/thermally_perfect.h"
#include "src/arrhenius.h"
#include "src/third_body.h"
#include "src/falloff_troe.h"
#include "src/falloff_lindemann.h"
#include "src/falloff_sri.h"
#include "src/pressure_dependent_arrhenius.h"
#include "src/reactions.h"
#include "src/source.h"

struct ChemgenInterface
{
    const fvMesh& mesh_;              // Reference to the mesh for creating volScalarFields

    // Constructor: Initialize members
    ChemgenInterface(const fvMesh& mesh)
        : mesh_(mesh) // Initialize with the provided mesh reference
    {}

    static Foam::stringList species_names() //fix later
    {
        // Hard-coded list of species names
        Foam::stringList species_names_(n_species);
        forAll(species_names_, i)
        {
            species_names_[i] = species_name_gen(i);
        }
        return species_names_;
    }

    static Foam::string species_name(const Foam::label specieI)
    {
        Foam::stringList species_names_ = species_names(); // Get the list of species names
        return species_names_[specieI];                   // Return the name of the requested species
    }

    // Overloaded function call operator to compute source terms
    PtrList<volScalarField> source_field
    (
        const PtrList<volScalarField>& Y,
        const volScalarField& rho,
        const volScalarField& T
    ) const
    {

        PtrList<volScalarField> sourceField_(Y.size()); // Create PtrList with the same size as Y

        forAll(Y, i)
        {
            sourceField_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "source_chemgen_" + Y[i].name(),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0) // Initialize with dimensions
                )
            );
        }


        Species inv_molecular_weights_ = inv_molecular_weights();
        Species concentrations = {scalar(0)};
        Species net_production_rates = {scalar(0)};

        forAll(T, cell_index)
        {
            for(label i=0; i<n_species; i++)
            {
                concentrations[i] = rho[cell_index] * Y[i][cell_index] * inv_molecular_weights_[i];
            }

            net_production_rates = source(concentrations, T[cell_index]) * molecular_weights();

            for(label i=0; i<n_species; i++)
            {
                sourceField_[i][cell_index] = net_production_rates[i];
            }
        }

        return sourceField_;
    }
    
    PtrList<volScalarField> concentrations
    (
        const PtrList<volScalarField>& Y,
        const volScalarField& rho,
        const volScalarField& T
    ) const
    {  
        Species inv_molecular_weights_ = inv_molecular_weights();
        
        PtrList<volScalarField> concentrations(Y.size()); // Create PtrList with the same size as Y

        forAll(Y, i)
        {
            concentrations.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "C" + Y[i].name(),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar("zero", dimMoles/dimVolume, 0.0) // Initialize with dimensions
                )
            );
        }


        forAll(T, cell_index)
        {
            for(label i=0; i<n_species; i++)
            {
                concentrations[i][cell_index] = rho[cell_index] * Y[i][cell_index] * inv_molecular_weights_[i];
            }
        }

        return concentrations;
    }


    PtrList<volScalarField> source_field
    (
        const PtrList<volScalarField>& C,
        const volScalarField& T
    ) const
    {

        PtrList<volScalarField> sourceField_(C.size()); // Create PtrList with the same size as Y

        forAll(C, i)
        {
            sourceField_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "source_chemgen_" + C[i].name(),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0) // Initialize with dimensions
                )
            );
        }


        Species inv_molecular_weights_ = inv_molecular_weights();
        Species concentrations = {scalar(0)};
        Species net_production_rates = {scalar(0)};

        forAll(T, cell_index)
        {
            for(label i=0; i<n_species; i++)
            {
                concentrations[i] = C[i][cell_index];
            }

            net_production_rates = source(concentrations, T[cell_index]) * molecular_weights();

            for(label i=0; i<n_species; i++)
            {
                sourceField_[i][cell_index] = net_production_rates[i];
            }
        }

        return sourceField_;
    }
        
    volScalarField density
    (
        const PtrList<volScalarField>& Y,
        const volScalarField& p,
        const volScalarField& T
    ) const
    {
        Species inv_MW = inv_molecular_weights();

        volScalarField rho_
        (
            IOobject
            (
                "rho",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            p.mesh(),
            dimensionedScalar("zero", dimMass/dimVolume, 0.0)
        );

        forAll(p, cellI)
        {
            scalar inv_MW_avg = scalar(0);
            forAll(Y, specieI)
            {
                inv_MW_avg+= Y[specieI][cellI] * inv_MW[specieI];
            }
            rho_[cellI] = p[cellI]/(T[cellI] * universal_gas_constant() * inv_MW_avg);
        }
        
        return rho_;
    }

    void density_from_concentrations
    (
        volScalarField& rho,
        const PtrList<volScalarField>& C
    ) const
    {

        Species MW = molecular_weights();

        forAll(rho, cellI)
        {
            scalar sum_rho = scalar(0);
            forAll(C, specieI)
            {
                sum_rho+= C[specieI][cellI] * MW[specieI];
            }

            rho[cellI] = sum_rho;
        } 
    }

    void mass_fractions_from_concentrations_density
    (
        PtrList<volScalarField>& Y,
        const PtrList<volScalarField>& C,
        const volScalarField& rho
    ) const
    {

        Species MW = molecular_weights();

        forAll(rho, cellI)
        {
            forAll(C, specieI)
            {
                Y[specieI][cellI] = C[specieI][cellI] * MW[specieI]/rho[cellI];
            }
        } 
    }

    void pressure_from_concentrations_temperature
    (
        volScalarField& p,
        const PtrList<volScalarField>& C,
        const volScalarField T
    ) const
    {
        forAll(p, cellI)
        {
            scalar sum_C = scalar(0);
            forAll(C, specieI)
            {
                sum_C+= C[specieI][cellI];
            }

            p[cellI] = T[cellI] * universal_gas_constant() * sum_C;
        }
    }


    volScalarField internal_energy
    (
        const PtrList<volScalarField>& Y,
        const volScalarField& rho,
        const volScalarField& T
    ) const
    {

        volScalarField internal_energy_
        (
            IOobject
            (
                "internal_energy",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rho.mesh(),
            dimensionedScalar("zero", dimEnergy/dimVolume, 0.0)
        );

        Species inv_molecular_weights_ = inv_molecular_weights();
        Species concentrations = {scalar(0)};
        Species net_production_rates = {scalar(0)};

        forAll(T, cell_index)
        {
            for(label i=0; i<n_species; i++)
            {
                concentrations[i] = rho[cell_index] * Y[i][cell_index] * inv_molecular_weights_[i];
            }
            internal_energy_[cell_index] = internal_energy_volume_specific(concentrations, T[cell_index]);
        }

        return internal_energy_;
    }
    
    void temperature_from_internal_energy
    (
        const PtrList<volScalarField>& C,
        const volScalarField& internal_energy,
        volScalarField& T
    ) const
    {
        Species concentrations = {scalar(0)};
        Species net_production_rates = {scalar(0)};

        forAll(T, cell_index)
        {
            for(label i=0; i<n_species; i++)
            {
                concentrations[i] = C[i][cell_index];
            }

            T[cell_index] = temperature(internal_energy[cell_index], concentrations, int(3)); 
        }
    }

};
