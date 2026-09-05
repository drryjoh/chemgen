 constexpr Species molecular_weights()  {return {double(31.998),double(54.092),double(28.009999999999998),double(44.009),double(18.015)};}
 constexpr Species inv_molecular_weights()  {return {double(0.03125195324707794),double(0.018487022110478445),double(0.03570153516601214),double(0.022722624917630486),double(0.055509297807382736)};}
    // Define the species names as a fixed-size array
    #pragma once
    #include <string>

    static constexpr std::array<const char*, 5> species_names_gen()
    {
        return {"O2", "C4H6", "CO", "CO2", "H2O"};
    }

    // Return the species name for a given index
    static std::string species_name_gen(int index)
    {
        constexpr auto names = species_names_gen(); // Get the list of species names use auto for now
        return names[index]; // Return the name of the requested species
    }
    // Return the species name for a given index
    int species_index_gen(const char* name)
    {
        constexpr auto names = species_names_gen(); // Get the list of species names use auto for now
        for(int i = 0; i<n_species; i++)
        {
            if (std::strcmp(names[i], name) == 0)
            {
                return i;
            }
        }
        return -1;
    }
    