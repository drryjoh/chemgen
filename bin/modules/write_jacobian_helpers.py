import numpy as np

def jacobian_output_text(configuration):
    if not configuration.eigen:
        return configuration.jacobian_function
    elif configuration.eigen_sparse:
        return "SparseMatrix<{scalar}>".format(**vars(configuration))
    else:
        return configuration.jacobian_function_eigen

def jacobian_initialize_text(configuration, sparsity_pattern):
    if not configuration.eigen:
        return "{jacobian} jacobian_net_production_rates = {{{scalar_cast}(0)}};".format(**vars(configuration))
    elif configuration.eigen_sparse:
        return """std::vector<Triplet_> jacobian_triplets;
        jacobian_triplets.reserve({nonzeros});""".format(**vars(configuration), nonzeros=np.sum(sparsity_pattern)+sparsity_pattern.shape[0]) # potentially over-allocate for diagonal fill
    else:
        return "{jacobian_eigen} jacobian_net_production_rates = {jacobian_eigen}::Zero();".format(**vars(configuration))

def update_jacobian_input(configuration):
    if not configuration.eigen:
        return "{jacobian}& jacobian_net_production_rates".format(**vars(configuration))
    elif configuration.eigen_sparse:
        return "std::vector<Triplet_>& jacobian_triplets".format(**vars(configuration))
    else:
        return "{jacobian_eigen}& jacobian_net_production_rates".format(**vars(configuration))

def modify_jacobian_text_eigen(configuration, formatted_text):
    if not configuration.eigen:
        return formatted_text
    elif configuration.eigen_sparse:
        text_split = formatted_text.replace('+=', '=').split('=')
        assert(len(text_split) == 2)
        lhs_text = text_split[0].strip()
        rhs_text = text_split[1].strip()
        ij_text = lhs_text.replace('jacobian_net_production_rates', '').replace(' ', '')
        ij_split = ij_text.strip('][').split('][')
        i_string = ij_split[0]
        j_string = ij_split[1]
        return "        jacobian_triplets.push_back(Triplet_({i}, {j}, {value}));\n".format(i=i_string, j=j_string, value=rhs_text.replace(';', ''))
    else:
        # only modify LHS
        text_split = formatted_text.split('=')
        assert(len(text_split) == 2)
        lhs_text = text_split[0]
        rhs_text = text_split[1]
        lhs_text = lhs_text.replace('][', ',',).replace('[', '(',).replace(']', ')',)
        return lhs_text + "=" + rhs_text

def jacobian_all_text_eigen(configuration, row_index, species_row_to_be_added_text):
    if not configuration.eigen:
        # TODO: modify jacobian_net_production_rates[{row_index}] in place?
        return f"        jacobian_net_production_rates[{row_index}] = add_species_to_chemical_state(jacobian_net_production_rates[{row_index}], {species_row_to_be_added_text});\n"
    elif configuration.eigen_sparse:
        return f"        add_species_to_jacobian_row_eigen_sparse(jacobian_triplets, {row_index}, {species_row_to_be_added_text});\n"
    else:
        return f"        add_species_to_jacobian_row_eigen(jacobian_net_production_rates, {row_index}, {species_row_to_be_added_text});\n"

def dtemperature_source_dspecies_and_dsource_species_dtemperature_text(configuration):
    if configuration.eigen_sparse:
        return """
        {scalar} specific_heat_constant_volume_volume_specific_ = specific_heat_constant_volume_volume_specific(species, temperature);
        Species dtemperature_source_dspecies_1 = scale_gen(divide(species_internal_energy_mole_source_sum(species, temperature), pow2(specific_heat_constant_volume_volume_specific_)), dspecific_heat_constant_volume_volume_specific_dspecies(species, temperature));

        for ({index} i = 0; i < n_species; i++)
        {{
            // Derivative of species source terms with respect to temperature
            jacobian_triplets.push_back(Triplet_(i+1, 0, dsource_species_dtemperature_[i]));

            // Only one term (out of two terms) of dtemperature_source_dspecies so the first row is fully allocated after compression
            // Will update with second term below
            jacobian_triplets.push_back(Triplet_(0, i+1, dtemperature_source_dspecies_1[i]));
        }}
        {diagonal_add}
        {set_jacobian_from_triplets}
        Species dspecies_internal_energy_mole_source_sum_dspecies_ = {{{scalar_cast}(0)}};
        Species species_internal_energy_mole_ = molecular_weights() * species_internal_energy_mass_specific(temperature);

        for ({index} i = 0; i < n_species; ++i)
        {{
            for (SparseMatrix<{scalar}>::InnerIterator it(jacobian_net_production_rates, i+1); it; ++it) // skip first column (temperature)
            {{
                {index} j = it.row();

                // skip first row (energy source)
                if (j == 0) continue;

                {scalar} value = it.value();

                // Undo diagonal_add
                if (i+1 == j)
                {{
                    value -= diagonal_add_;
                }}

                dspecies_internal_energy_mole_source_sum_dspecies_[i] += species_internal_energy_mole_[j-1] * value;
                // it.value(); // value
                // it.row();   // row index
                // it.col();   // col index - here, it is equal to i+1
                // it.index(); // inner index - here, it is equal to it.row()
            }}
        }}

        for ({index} i = 0; i < n_species; i++)
        {{
            // Derivative of temperature source term with respect to concentrations
            // Second part
            for (SparseMatrix<{scalar}>::InnerIterator it(jacobian_net_production_rates, i+1); it; ++it) // skip first column (temperature)
            {{
                if (it.row() != 0) std::cerr << "it.row() != 0" << std::endl; // first row is dense

                it.valueRef() = it.value() - divide(dspecies_internal_energy_mole_source_sum_dspecies_[i], specific_heat_constant_volume_volume_specific_);

                break; // only first row
            }}
        }}
""".format(**vars(configuration), diagonal_add=add_diagonal_text(configuration), set_jacobian_from_triplets=set_jacobian_from_triplets_text(configuration))

    else:
        return """
        {species} dspecies_internal_energy_mole_source_sum_dspecies_ = {{{scalar_cast}(0)}};
        {species} species_internal_energy_mole_ = molecular_weights() * species_internal_energy_mass_specific(temperature);

        for ({index} i = 0; i < n_species; i++)
        {{
            {species} jacobian_column;
            for ({index} j = 0; j < n_species; j++)
            {{
                jacobian_column[j] = jacobian_net_production_rates{left}j+1{mid}i+1{right};
            }}
            dspecies_internal_energy_mole_source_sum_dspecies_[i] = dot(species_internal_energy_mole_, jacobian_column);
        }}

        {species} dtemperature_source_dspecies_ = dtemperature_source_dspecies(temperature, species, dspecies_internal_energy_mole_source_sum_dspecies_);

        for ({index} i = 0; i < n_species; i++)
        {{
            // Derivative of temperature source term with respect to concentrations
            {store_dtemperature_source_dspecies}
            // Derivative of species source terms with respect to temperature
            {store_dsource_species_dtemperature}
        }}
""".format(**vars(configuration),
           store_dtemperature_source_dspecies=modify_jacobian_text_eigen(configuration, "jacobian_net_production_rates[0][i+1] = dtemperature_source_dspecies_[i];"),
           store_dsource_species_dtemperature=modify_jacobian_text_eigen(configuration, "jacobian_net_production_rates[i+1][0] = dsource_species_dtemperature_[i];"),
           scale_jacobian="= scale_gen(scaling_factor, jacobian_net_production_rates)" if not configuration.eigen else "*= scaling_factor")

def scale_jacobian_text(configuration):
    return """
        // Scale Jacobian
        if (scaling_factor != 1)
        {{
            jacobian_net_production_rates {scale_jacobian};
        }}
""".format(**vars(configuration), scale_jacobian="= scale_gen(scaling_factor, jacobian_net_production_rates)" if not configuration.eigen else "*= scaling_factor")

def add_diagonal_text(configuration):
    if configuration.eigen_sparse:
        return """
        // Divide diagonal_add by scaling_factor here since we will do J = scaling_factor*J at the end
        {scalar} diagonal_add_ = divide(diagonal_add, scaling_factor);

        // Add to diagonal - append to triplets (before scaling by scaling_factor) for two reasons:
        // (1) Ensure memory is allocated for the diagonal
        // (2) So we don't have to access solely the diagonal later
        if (diagonal_add != 0)
        {{
            for ({index} i = 0; i < n_variables; i++)
            {{
                jacobian_triplets.push_back(Triplet_(i, i, diagonal_add_));
            }}
        }}
""".format(**vars(configuration))

    else:
        return """
        // Add to diagonal
        if (diagonal_add != 0)
        {{
            for ({index} i = 0; i < n_variables; i++)
            {{
                jacobian_net_production_rates{left}i{mid}i{right} += diagonal_add;
            }}
        }}
""".format(**vars(configuration))

def write_scale_jacobian(file, configuration):
    file.write(scale_jacobian_text(configuration))

def write_add_diagonal(file, configuration):
    if configuration.eigen_sparse:
        # diagonal add dealt with separately (either in write_helpers_eigen_sparse or dtemperature_source_dspecies_and_dsource_species_dtemperature_text)
        return

    file.write(add_diagonal_text(configuration))

def set_jacobian_from_triplets_text(configuration):
    return """
        SparseMatrix<{scalar}> jacobian_net_production_rates(n_variables, n_variables);
        jacobian_net_production_rates.setFromTriplets(jacobian_triplets.begin(), jacobian_triplets.end());
""".format(**vars(configuration))

def write_helpers_eigen_sparse(file, temperature_equation, configuration):
    assert(configuration.eigen_sparse)
    
    if not temperature_equation:
        file.write(add_diagonal_text(configuration))
        file.write(set_jacobian_from_triplets_text(configuration))
